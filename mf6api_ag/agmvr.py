from modflowapi import Callbacks, run_simulation
import numpy as np
import flopy
import os


class ModflowAgmvr(object):
    """
    Modflow6 API AG pacakge that uses the MVR to facilitate irrigation

    Package is based on the ETDEMAND AG method in MODFLOW-NWT

    Parameters
    ----------
    sim : flopy.mf6.MFSimulation object
    ag_type : str
        "etdemand" is the only supported method currently
    mvr_name : str
        short name for the Modflow6 MVR package. Ex. "mvr_0"
    """

    def __init__(self, sim, ag_type="etdemand", mvr_name="mvr"):
        self.sim = sim
        name = list(sim.model_names)[0]
        self.name = name.upper()
        self.gwf = sim.get_model(name)
        self.mvr_name = mvr_name
        self.mvr = self.gwf.get_package(self.mvr_name)
        if self.mvr is None:
            raise AssertionError(
                "MVR or AGMVR package does not exist for this model"
            )
        if not isinstance(self.mvr, flopy.mf6.ModflowGwfmvr):
            # assume this is an Agmvr package!
            # need to pop agmvr, replace with mvr package and re-write
            # to file
            perioddata = {}
            for i, recarray in self.mvr.perioddata.data.items():
                record = []
                for rec in recarray:
                    record.append(
                        (
                            rec["pname1"],
                            rec["id1"],
                            rec["pname2"],
                            rec["id2"],
                            "UPTO",
                            rec["value"],
                        )
                    )
                perioddata[i] = record

            self.gwf.remove_package(self.mvr_name)
            mvr = flopy.mf6.ModflowMvr(
                self.gwf,
                maxmvr=self.mvr.maxmvr.array,
                maxpackages=self.mvr.maxpackages.array,
                packages=self.mvr.packages.array,
                perioddata=perioddata,
                budget_filerecord=self.mvr.budget_filerecord.array,
                budgetcsv_filerecord=self.mvr.budgetcsv_filerecord.array,
                pname=self.mvr.package_name,
            )

            self.sim.write_simulation()
            self.mvr_name = mvr.path[-1]

        if ag_type.lower() in ("specified", "trigger", "etdemand"):
            self.ag_type = ag_type.lower()

        self.sim_wells = False
        self.sim_maw = False
        self.sim_diversions = False
        self.sim_lak = False

        self.sfr_evap_var = "EVAP"

        self.lak_outrate_var = "OUTRATE"
        self.lak_sarea_var = "SAREA"
        self.lak_evap_var = "EVAPORATION"

    def _set_package_stress_period_data(self, ml, pkg):
        """
        General method to set stress period data from MVR
        (MAXQ, Delivery cells, etc)

        Parameters
        ----------
        mf6 : ModflowApi object
        kper : int
            zero based stress period number
        pkg : str
            package abbreviation (ex. well, maw, sfr)

        Return
        ------
        tuple : shape of max_q for dimensionalizing sup, supold, and aetold
        """
        heads = ml.X
        botm = ml.dis.botm.values
        pet = mf6.get_value(self.pet_addr)
        if pkg in ("well", "maw"):
            addr = self.__dict__[f"{pkg}_addr"]
            data = mf6.get_value(addr)
            if pkg == "well":
                max_q = np.abs(np.copy(data.T[0]))
                data[:, 0] = 0
            else:
                max_q = np.abs(np.copy(data))
                data[:] = 0

            mf6.set_value(addr, data)
        else:
            max_addr = self.__dict__[f"{pkg}q_for_mvr_addr"]
            qformvr = mf6.get_value(max_addr)
            max_q = np.zeros(np.abs(qformvr.shape))

        mvr = mf6.get_value(self.mvr_value_addr)
        pkg_name = self.__dict__[f"{pkg}_name"]
        if kper in self.mvr.perioddata.data:
            recarray = self.mvr.perioddata.data[kper]
            irridx = np.where((recarray["pname1"] == pkg_name) & (recarray["pname2"] == self.uzf_name.lower()))[0]
            if len(irridx) > 0:
                if pkg_name in ("sfr", "lak"):
                    mvr[irridx] = 0

                irrids = sorted(np.unique(recarray[irridx]["id1"]))
            else:
                irrids = []

            mf6.set_value(self.mvr_value_addr, mvr)
        else:
            irrids = []
            recarray = np.recarray((0,), dtype=[("id1", int),
                                                ("pname1", object),
                                                ("id2", int),
                                                ("pname2", object),
                                                ("value", float)])

        active = np.zeros(max_q.shape, dtype=int)
        active[irrids] = 1

        irrigated_cells = []
        irrigated_proportion = []
        irrigation_efficiency = []
        application_fraction = []
        mvr_index = []

        for irrid in range(max_q.size):
            icells = []
            iprop = []
            irreff = []
            appfrac = []
            idx = np.where(
                (recarray["id1"] == irrid) & (recarray["pname1"] == pkg_name)
            )[0]
            if len(idx) > 0:
                icells = recarray[idx]["id2"]
                iprop = recarray[idx]["value"] / np.sum(recarray[idx]["value"])
                if pkg in ("sfr", "lak"):
                    for mvr_rec in idx:
                        max_q[irrid] += recarray[mvr_rec]["value"]

                if "irr_eff" in recarray.dtype.names:
                    irreff = recarray[idx]["irr_eff"]
                    appfrac = recarray[idx]["app_frac"]
                else:
                    irreff = np.ones((len(icells),))
                    appfrac = np.ones((len(icells),))

            irrigation_efficiency.append(irreff)
            application_fraction.append(appfrac)
            irrigated_cells.append(icells)
            irrigated_proportion.append(iprop)
            mvr_index.append(idx)

        # adjustment for multiple providers in same pkg irrigating a node
        if pkg == "well":
            nodelist = mf6.get_value(self.well_nodelist)
            hk = mf6.get_value(self.k11_addr)

        wf_adj = np.ones((max_q.shape[0], pet.shape[0]))
        cnode_provider = {}
        cnode_weights = {}
        max_nodes = 1
        for ix, crop_nodes in enumerate(irrigated_cells):
            if len(crop_nodes) > 0:
                for node in crop_nodes:
                    if pkg == "well":
                        gw_node = nodelist[ix] - 1
                        if heads[gw_node] <= -1e+30 or (heads[gw_node] - botm[gw_node]) < 0:
                            weight = 0
                        else:
                            weight = (heads[gw_node] - botm[gw_node]) * hk[gw_node]
                    else:
                        weight = max_q[ix]
                    if node in cnode_provider:
                        cnode_provider[node].append(ix)
                        cnode_weights[node].append(weight)
                    else:
                        cnode_provider[node] = [ix]
                        cnode_weights[node] = [weight]

        for node, provider in cnode_provider.items():
            weight = np.array(cnode_weights[node]) / np.sum(
                cnode_weights[node])
            wf_adj[provider, node] = weight

        setattr(self, f"{pkg}_active", active)
        setattr(self, f"{pkg}_maxq", max_q)
        setattr(self, f"{pkg}_irrigated_cells", irrigated_cells)
        setattr(self, f"{pkg}_irrigated_proportion", irrigated_proportion)
        setattr(self, f"{pkg}_irrigation_efficiency", irrigation_efficiency)
        setattr(self, f"{pkg}_application_fraction", application_fraction)
        setattr(self, f"{pkg}_mvr_index", mvr_index)
        setattr(self, f"{pkg}sup", np.zeros(max_q.shape))
        setattr(self, f"{pkg}supold", np.zeros(max_q.shape))
        setattr(self, f"{pkg}aetold", np.zeros(max_q.shape))
        setattr(self, f"{pkg}mp_weights", wf_adj)

    def set_stress_period_data(self, ml):
        """
        Method to set stress period data from MVR (MAXQ, Delivery cells, etc)

        Parameters
        ----------
        ml : ModflowApi.ApiModel object

        """
        if self.sim_wells:
            self._set_package_stress_period_data(ml, "well")

        if self.sim_maw:
            self._set_package_stress_period_data(ml, "maw")

        if self.sim_diversions:
            self._set_package_stress_period_data(ml, "sfr")
            self.sfr_evap_old = ml.sfr.get_advanced_var(self.sfr_evap_var)

        if self.sim_lak:
            self._set_package_stress_period_data(ml, "lak")
            self.lak_evap_old = ml.lak.get_advanced_var(self.lak_evap_var)

        if ml.kper in self.mvr.perioddata.data:
            t = self.mvr.perioddata.data[ml.kper]
            self.ag_active = np.where(
                t["pname2"] == self.uzf_name.lower(), True, False
            )
        else:
            pass

    def callback(self, sim, cstep):
        """
        Callback function for running the agmvr model
        """
        # todo: could support multimodel by checking if UZF is active
        #   and if the MVR package is also active
        # todo: could support ET package and simple gwet calculation
        #   if UZF is not active
        ml = sim.get_model()
        if cstep == Callbacks.stress_period_start:
            self.set_stress_period_data(ml)