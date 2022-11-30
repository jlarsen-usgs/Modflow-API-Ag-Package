from modflowapi import ModflowApi
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

        self.totim = np.add.accumulate(self.gwf.modeltime.perlen)

        self.sim_wells = False
        self.sim_maw = False
        self.sim_diversions = False
        self.sim_lak = False

        self.well_name = None
        self.maw_name = None
        self.sfr_name = None
        self.lak_name = None
        self.uzf_name = None
        self.dis_name = None
        self.output_filename = os.path.join(
            self.gwf.model_ws, f"{self.name.lower()}_ag.out"
        )

        self.get_pkgs()

    def get_pkgs(self):
        """
        Method to get AG package names

        """
        pkg_names = self.mvr.packages

        well = None
        sfr = None
        uzf = None
        maw = None
        lak = None
        for ix, name in enumerate(pkg_names.array):
            name = name[-1]
            pkg = self.gwf.get_package(name)
            if pkg is None:
                pkg = self.gwf.get_package(f"{name[:3]}_0")
            if isinstance(pkg, flopy.mf6.ModflowGwfwel):
                well = name
                self.sim_wells = True
            elif isinstance(pkg, flopy.mf6.ModflowGwfmaw):
                maw = name
                self.sim_maw = True
            elif isinstance(pkg, flopy.mf6.ModflowGwfsfr):
                sfr = name
                self.sim_diversions = True
            elif isinstance(pkg, flopy.mf6.ModflowGwflak):
                lak = name
                self.sim_lak = True
            elif isinstance(pkg, flopy.mf6.ModflowGwfuzf):
                uzf = name
            else:
                pass

        if self.sim_wells:
            self.well_name = well
            self.well_output = {}
        if self.sim_maw:
            self.maw_name = maw
            self.maw_output = {}
        if self.sim_diversions:
            self.sfr_name = sfr
            self.sfr_output = {}
        if self.sim_lak:
            self.lak_name = lak
            self.lak_output = {}

        if uzf is None:
            raise AssertionError("UZF must be active to simulate AG")

        self.uzf_name = uzf
        self.dis_name = self.gwf.dis.name[0].upper()
        self.npf_name = self.gwf.npf.name[0].upper()
        self.ag_active = np.zeros((1,))

    def create_addresses(self, mf6):
        """
        Method to create MODFLOW-API addresses

        Parameters
        ----------
        mf6 : ModflowApi object

        """
        solutiongroup = self.sim.name_file.solutiongroup.data
        sln_grp = 1
        for sln_grp, recarray in solutiongroup.items():
            if self.name.lower() in recarray["slnmnames"]:
                sln_grp += 1
                break

        self.maxiter_addr = mf6.get_var_address("MXITER", f"SLN_{sln_grp}")
        self.head_addr = mf6.get_var_address("X", self.name)
        self.botm_addr = mf6.get_var_address(
            "BOT", self.name, self.dis_name.upper()
        )
        self.k11_addr = mf6.get_var_address(
            "K11", self.name, self.npf_name.upper()
        )
        self.area_addr = mf6.get_var_address(
            "UZFAREA", self.name, self.uzf_name.upper()
        )
        self.vks_addr = mf6.get_var_address(
            "VKS", self.name, self.uzf_name.upper()
        )
        self.pet_addr = mf6.get_var_address(
            "PET", self.name, self.uzf_name.upper()
        )
        self.aet_addr = mf6.get_var_address(
            "ETACT", self.name, self.uzf_name.upper()
        )
        self.gwet_addr = mf6.get_var_address(
            "GWET", self.name, self.uzf_name.upper()
        )
        self.uzf_mvr_addr = mf6.get_var_address(
            "QFROMMVR", self.name, self.uzf_name.upper()
        )
        self.mvr_id1_addr = mf6.get_var_address(
            "ID1", self.name, self.mvr_name.upper()
        )
        self.mvr_id2_addr = mf6.get_var_address(
            "ID2", self.name, self.mvr_name.upper()
        )
        self.mvr_type_addr = mf6.get_var_address(
            "IMVRTYPE", self.name, self.mvr_name.upper()
        )
        self.mvr_value_addr = mf6.get_var_address(
            "VALUE", self.name, self.mvr_name.upper()
        )

        self.uzf_shape = mf6.get_value(self.uzf_mvr_addr).shape
        self.applied_irrigation = np.zeros(self.uzf_shape)

        if self.sim_wells:
            self.wellq_for_mvr_addr = mf6.get_var_address(
                "QFORMVR", self.name, self.well_name.upper()
            )
            self.well_addr = mf6.get_var_address(
                "BOUND", self.name, self.well_name.upper()
            )

            self.well_nodelist = mf6.get_var_address(
                "NODELIST", self.name, self.well_name.upper()
            )

        if self.sim_maw:
            self.mawq_for_mvr_addr = mf6.get_var_address(
                "QFORMVR", self.name, self.maw_name.upper()
            )
            self.maw_addr = mf6.get_var_address(
                "RATE", self.name, self.maw_name.upper()
            )

        if self.sim_diversions:
            self.sfrq_for_mvr_addr = mf6.get_var_address(
                "QFORMVR", self.name, self.sfr_name.upper()
            )
            self.sfrq_to_mvr_addr = mf6.get_var_address(
                "QTOMVR", self.name, self.sfr_name.upper()
            )
            self.sfr_evap_addr = mf6.get_var_address(
                "EVAP", self.name, self.sfr_name.upper()
            )
            self.sfr_length_addr = mf6.get_var_address(
                "LENGTH", self.name, self.sfr_name.upper()
            )
            self.sfr_width_addr = mf6.get_var_address(
                "WIDTH", self.name, self.sfr_name.upper()
            )
        if self.sim_lak:
            self.lak_addr = mf6.get_var_address(
                "OUTRATE", self.name, self.lak_name.upper()
            )
            self.lakq_for_mvr_addr = mf6.get_var_address(
                "QFORMVR", self.name, self.lak_name.upper()
            )
            self.lakq_to_mvr_addr = mf6.get_var_address(
                "QTOMVR", self.name, self.lak_name.upper()
            )
            self.lak_sarea_addr = mf6.get_var_address(
                "SAREA", self.name, self.lak_name.upper()
            )
            self.lak_evap_addr = mf6.get_var_address(
                "EVAPORATION", self.name, self.lak_name.upper()
            )

    def _set_package_stress_period_data(self, mf6, kper, pkg):
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
        heads = mf6.get_value(self.head_addr)
        botm = mf6.get_value(self.botm_addr)
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

    def set_stress_period_data(self, mf6, kper):
        """
        Method to set stress period data from MVR (MAXQ, Delivery cells, etc)

        Parameters
        ----------
        mf6 : ModflowApi object
        kper : int
            zero based stress period number
        """
        if self.sim_wells:
            self._set_package_stress_period_data(mf6, kper, "well")

        if self.sim_maw:
            self._set_package_stress_period_data(mf6, kper, "maw")

        if self.sim_diversions:
            self._set_package_stress_period_data(mf6, kper, "sfr")
            self.sfr_evap_old = mf6.get_value(self.sfr_evap_addr)

        if self.sim_lak:
            self._set_package_stress_period_data(mf6, kper, "lak")
            self.lak_evap_old = mf6.get_value(self.lak_evap_addr)
        if kper in self.mvr.perioddata.data:
            t = self.mvr.perioddata.data[kper]
            self.ag_active = np.where(
                t["pname2"] == self.uzf_name.lower(), True, False
            )
        else:
            pass

    def zero_mvr(self, mf6):
        """
        Method to zero out MVR values before initial solve
        """
        mvr = mf6.get_value(self.mvr_value_addr)
        mvr = np.where(self.ag_active, 0, mvr)
        mf6.set_value(self.mvr_value_addr, mvr)

    def _set_etdemand_variables(self, mf6, pkg):
        """
        General method to calculate and set common etdemand calculation
        variables

        Parameters
        ----------
        mf6 : ModflowApi object
        pkg : str
            package type abbr ("well", "maw", "sfr")


        Returns
        -------
        """
        pet = mf6.get_value(self.pet_addr)
        vks = mf6.get_value(self.vks_addr)
        area = np.tile(mf6.get_value(self.area_addr), self.gwf.modelgrid.nlay)
        aet = mf6.get_value(self.aet_addr)
        gwet = mf6.get_value(self.gwet_addr)
        maxq = getattr(self, f"{pkg}_maxq")
        application_fraction = getattr(self, f"{pkg}_application_fraction")
        irrigated_cells = getattr(self, f"{pkg}_irrigated_cells")
        mp_adj = getattr(self, f"{pkg}mp_weights")
        nodelist = None
        gw_avail = None
        if pkg == "well":
            nodelist = mf6.get_value(self.well_nodelist)
            head = mf6.get_value(self.head_addr)
            botm = mf6.get_value(self.botm_addr)
            gw_avail = np.zeros(maxq.shape)

        crop_vks = np.zeros(maxq.shape)
        crop_pet = np.zeros(maxq.shape)
        crop_aet = np.zeros(maxq.shape)
        crop_gwet = np.zeros(maxq.shape)
        app_frac = np.zeros(maxq.shape)
        prev_applied = np.zeros(maxq.shape)
        for ix, crop_nodes in enumerate(irrigated_cells):
            if len(crop_nodes) > 0:
                crop_vks[ix] = np.sum(vks[crop_nodes] * area[crop_nodes])
                crop_pet[ix] = np.sum(
                    pet[crop_nodes] * area[crop_nodes] * mp_adj[ix, crop_nodes]
                )
                crop_aet[ix] = np.sum(
                    aet[crop_nodes] * area[crop_nodes] * mp_adj[ix, crop_nodes]
                )
                crop_gwet[ix] = np.sum(gwet[crop_nodes] * mp_adj[ix, crop_nodes])
                app_frac[ix] = np.mean(application_fraction[ix])
                prev_applied[ix] = np.sum(self.applied_irrigation[crop_nodes])
                if pkg == "well":
                    node = nodelist[ix] - 1
                    if head[node] <= -1e+30 or (head[node] - botm[node]) < 0:
                        gw_avail[ix] = 0
                    else:
                        gw_avail[ix] = (head[node] - botm[node]) * area[node]

       # gw_avail = np.where(gw_avail > 0, gw_avail, 0)
        crop_aet = np.where(np.isnan(crop_gwet), crop_aet, crop_aet + crop_gwet)
        if pkg in ("well", "maw"):
            crop_aet = np.where(crop_vks < 1e-30, crop_pet, crop_aet)

        return crop_pet, crop_aet, crop_vks, app_frac, prev_applied, gw_avail

    def _gw_demand_etdemand(self, mf6, kstp, delt, kiter, pkg):
        """
        Generalized method to determine groundwater use demand

        Parameters
        ----------
        mf6 : ModflowApi object
        kstp : int
            modflow6 time step
        delt : float
            length of time step
        kiter : int
            iteration number
        """
        (
            crop_pet,
            crop_aet,
            crop_vks,
            app_frac,
            prev_applied,
            gw_avail
        ) = self._set_etdemand_variables(mf6, pkg)

        maxq = getattr(self, f"{pkg}_maxq")
        pkg_sup = getattr(self, f"{pkg}sup")
        pkg_supold = getattr(self, f"{pkg}supold")
        pkg_aetold = getattr(self, f"{pkg}aetold")
        pkg_active = getattr(self, f"{pkg}_active")
        pkg_addr = getattr(self, f"{pkg}_addr")
        pkg_mvr_index = getattr(self, f"{pkg}_mvr_index")
        application_fraction = getattr(self, f"{pkg}_application_fraction")
        irrigated_proportion = getattr(self, f"{pkg}_irrigated_proportion")
        irrigation_efficiency = getattr(self, f"{pkg}_irrigation_efficiency")
        irrigated_cells = getattr(self, f"{pkg}_irrigated_cells")
        pkg_output = getattr(self, f"{pkg}_output")
        pkg_name = getattr(self, f"{pkg}_name")

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            sup = np.zeros(maxq.shape)
            if kstp == 0:
                supold = np.zeros(maxq.shape)
            else:
                supold = pkg_supold
        else:
            sup = pkg_sup
            supold = pkg_supold

        factor = self._calculate_factor(
            crop_pet, crop_aet, pkg_aetold, sup, supold, kiter
        )
        factor *= app_frac

        qonly = np.where(sup + factor > crop_vks, crop_vks, sup + factor)

        setattr(self, f"{pkg}supold", sup)
        setattr(self, f"{pkg}aetold", crop_aet)

        if self.sim_diversions:
            qonly = qonly - prev_applied
            qonly = np.where(qonly < 0, 0, qonly)

        pumping = np.where(qonly > maxq, -1 * maxq, -1 * qonly)
        if gw_avail is not None:
            # adjustment for well demand based on gw availability
            pumping = np.where(np.abs(pumping) > gw_avail, -1 * gw_avail, pumping)

        pumping = np.where(np.abs(pumping) <= 1e-10, 0, pumping)

        setattr(self, f"{pkg}sup", np.abs(pumping))

        active_ix = np.where(pkg_active)[0]

        if len(active_ix) > 0:
            wells = mf6.get_value(pkg_addr)
            if pkg in ("well",):
                wells[active_ix, 0] = pumping[active_ix]
            else:
                wells[active_ix] = pumping[active_ix]
            mf6.set_value(pkg_addr, wells)

            mvr = mf6.get_value(self.mvr_value_addr)
            for well in active_ix:
                idx = pkg_mvr_index[well]
                app_frac_proportion = (
                    application_fraction[well]
                    / np.sum(application_fraction[well])
                ) / (1 / len(idx))
                mvr[idx] = (
                    (np.abs(pumping[well]) * irrigated_proportion[well])
                    * app_frac_proportion
                    * irrigation_efficiency[well]
                )
                self.applied_irrigation[irrigated_cells[well]] = (
                    (np.abs(pumping[well]) * irrigated_proportion[well])
                    * app_frac_proportion
                    * irrigation_efficiency[well]
                )

            mf6.set_value(self.mvr_value_addr, mvr)

            # store output...
            stp_output = []
            kkstp = mf6.get_value(mf6.get_var_address("KSTP", "TDIS"))[0]
            kper = mf6.get_value(mf6.get_var_address("KPER", "TDIS"))[0]
            for well in active_ix:
                idx = pkg_mvr_index[well]
                rec = (
                    kstp,
                    kper,
                    kkstp,
                    pkg_name,
                    well + 1,
                    crop_pet[well],
                    crop_aet[well],
                    sup[well] + factor[well],
                    np.abs(pumping[well]),
                    np.sum(mvr[idx]),
                )
                stp_output.append(rec)

            pkg_output[kstp] = stp_output
            setattr(self, f"{pkg}_output", pkg_output)

    def well_demand_etdemand(self, mf6, kstp, delt=1, kiter=1):
        """
        Method to determine groundwater use demand for well pkg irrigation

        Parameters
        ----------
        mf6 : ModflowApi object
        kstp : int
            modflow6 time step
        delt : float
            length of time step
        kiter : int
            iteration number
        """
        self._gw_demand_etdemand(mf6, kstp, delt, kiter, "well")

    def maw_demand_etdemand(self, mf6, kstp, delt=1, kiter=1):
        """
        Method to determine groundwater use demand for multi aquifer wells

        Parameters
        ----------
        mf6 : ModflowApi object
        kstp : int
            modflow6 time step
        delt : float
            length of time step
        kiter : int
            iteration number
        """
        self._gw_demand_etdemand(mf6, kstp, delt, kiter, "maw")

    def _sw_demand_etdemand(self, mf6, kstp, delt, kiter, pkg):
        """
        Method to determine surface-water demand

        Parameters
        ----------
        mf6 : ModflowApi object
        kstp : int
            modflow6 time step
        delt : float
            length of time step
        kiter : int
            iteration number
        pkg : str
            package name ("sfr" or "lak")
        """
        (
            crop_pet,
            crop_aet,
            crop_vks,
            app_frac,
            prev_applied,
            gw_avail
        ) = self._set_etdemand_variables(mf6, pkg)

        pkgq_to_mvr_addr = getattr(self, f"{pkg}q_to_mvr_addr")
        pkgq_for_mvr_addr = getattr(self, f"{pkg}q_for_mvr_addr")
        maxq = getattr(self, f"{pkg}_maxq")
        qtomvr = mf6.get_value(pkgq_to_mvr_addr)
        qformvr = mf6.get_value(pkgq_for_mvr_addr)
        pkg_active = getattr(self, f"{pkg}_active")
        pkg_mvr_index = getattr(self, f"{pkg}_mvr_index")
        pkg_evap_addr = getattr(self, f"{pkg}_evap_addr")
        application_fraction = getattr(self, f"{pkg}_application_fraction")
        irrigated_proportion = getattr(self, f"{pkg}_irrigated_proportion")
        irrigation_efficiency = getattr(self, f"{pkg}_irrigation_efficiency")
        irrigated_cells = getattr(self, f"{pkg}_irrigated_cells")
        evap = getattr(self, f"{pkg}_evap_old").copy()
        pkg_output = getattr(self, f"{pkg}_output")
        pkg_name = getattr(self, f"{pkg}_name")

        if pkg == "sfr":
            length = mf6.get_value(self.sfr_length_addr)
            width = mf6.get_value(self.sfr_width_addr)
            sarea = length * width
        else:
            sarea = mf6.get_value(self.lak_sarea_addr)

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            setattr(self, f"{pkg}sup", np.zeros(maxq.shape))
            setattr(self, f"{pkg}supold", np.zeros(maxq.shape))
            setattr(self, f"{pkg}aetold", crop_aet)
            qtomvr[:] = 0
            if kstp == 0:
                setattr(self, f"{pkg}aetold", np.zeros(crop_aet.shape))

        pkg_sup = qtomvr
        pkg_supold = getattr(self, f"{pkg}supold")
        pkg_aetold = getattr(self, f"{pkg}aetold")
        factor = self._calculate_factor(
            crop_pet, crop_aet, pkg_aetold, pkg_sup, pkg_supold, kiter
        )
        factor *= app_frac

        qonly = np.where(
            pkg_sup + factor > crop_vks, crop_vks, pkg_sup + factor
        )
        factor = np.where(factor < 0, 0, factor)

        if self.sim_lak and pkg == "sfr":
            only = qonly - prev_applied
            qonly = np.where(qonly < 0, 0, qonly)

        pkg_supold[:] = qtomvr[:]
        pkg_sup += factor
        pkg_aetold = crop_aet

        setattr(self, f"{pkg}supold", pkg_supold)
        setattr(self, f"{pkg}sup", pkg_sup)
        setattr(self, f"{pkg}aetold", pkg_aetold)

        dvflw = np.where(qonly >= maxq, maxq, qonly)
        dvflw = np.where(dvflw > qformvr, qformvr, dvflw)
        active_ix = np.where(pkg_active)[0]

        diversions = mf6.get_value(self.mvr_value_addr)
        tot_div_requested = np.zeros(maxq.shape)
        for provider in active_ix:
            idx = pkg_mvr_index[provider]
            app_frac_proportion = (
                application_fraction[provider]
                / np.sum(application_fraction[provider])
            ) / (1 / len(idx))
            div_requested = (
                dvflw[provider] * irrigated_proportion[provider]
            ) * app_frac_proportion
            tot_div_requested[provider] = np.sum(div_requested)
            div_inefficient = div_requested * irrigation_efficiency[provider]
            diversions[idx] = div_inefficient
            evap[provider] += (
                np.sum(div_requested - div_inefficient) / sarea[provider]
            )
            self.applied_irrigation[irrigated_cells[provider]] = (
                dvflw[provider] * irrigated_proportion[provider]
            ) * app_frac_proportion
        mf6.set_value(self.mvr_value_addr, diversions)
        mf6.set_value(pkg_evap_addr, evap)

        stp_output = []
        kkstp = mf6.get_value(mf6.get_var_address("KSTP", "TDIS"))[0]
        kper = mf6.get_value(mf6.get_var_address("KPER", "TDIS"))[0]
        for provider in active_ix:
            idx = pkg_mvr_index[provider]
            rec = (
                kstp,
                kper,
                kkstp,
                pkg_name,
                provider + 1,
                crop_pet[provider],
                crop_aet[provider],
                pkg_sup[provider],
                tot_div_requested[provider],
                np.sum(diversions[idx]),
            )
            stp_output.append(rec)

        pkg_output[kstp] = stp_output
        setattr(self, f"{pkg}_output", pkg_output)

    def sfr_demand_etdemand(self, mf6, kstp, delt=1, kiter=1):
        """
        Method to determine surface-water demand from SFR

        Parameters
        ----------
        mf6 : ModflowApi object
        kstp : int
            modflow6 time step
        delt : float
            length of time step
        kiter : int
            iteration number
        """
        self._sw_demand_etdemand(mf6, kstp, delt, kiter, "sfr")

    def lak_demand_etdemand(self, mf6, kstp, delt=1, kiter=1):
        """
        Method to determine surface-water demand from LAK

        Parameters
        ----------
        mf6 : ModflowApi object
        kstp : int
            modflow6 time step
        delt : float
            length of time step
        kiter : int
            iteration number
        """
        self._sw_demand_etdemand(mf6, kstp, delt, kiter, "lak")

    def _calculate_factor(
        self, crop_pet, crop_aet, aetold, sup, supold, kiter, accel=1
    ):
        """
        Method to calculate unmet demand based on PET and AET difference for
        crops

        Parameters
        ----------
        crop_pet : np.ndarray
            array of crop potential ET values
        crop_aet : np.ndarray
            array of crop actual ET values
        aetold : np.ndarray
            array of actual ET from previous iteration
        sup : np.ndarray
            current demand
        supold : np.ndarray
            demand from previous iteration
        kiter : float
            iteration number
        accel : float
        """
        et_diff = crop_pet - crop_aet
        factor = et_diff
        det = crop_aet - aetold
        dq = sup - supold

        if kiter > 0:
            # original equation logic commented out. I don't think
            # that we should continuously apply et_diff if det is small
            # instead there should be zero change, this prevents runaway
            # pumping from the equation...
            # factor = np.where(np.abs(det) > 1e-05,
            #                   dq * et_diff / det,
            #                   factor)
            factor = np.where(np.abs(det) > 1e-6, dq * et_diff / det, 1e-05)

        factor = np.where(factor > et_diff * accel, et_diff * accel, factor)
        # original equation logic commented out. I believe this to be one
        # of the culprits for causing runaway pumping values!
        # factor = np.where(factor < et_diff, et_diff, factor)
        factor = np.where(factor < 1e-05, 1e-05, factor)
        return factor

    def _write_output(self):
        """
        Method for writing AG output to file
        """
        with open(self.output_filename, "w") as foo:
            header = [
                "kstp",
                "kper",
                "stp",
                "pkg",
                "pid",
                "pet",
                "aet",
                "etdemand",
                "q_from_provider",
                "q_to_receiver",
            ]
            hdr_str = "{:>8} {:>8} {:>8} {:>10} {:>8} {:>15} {:>15} {:>15} {:>15} {:>15}\n"
            foo.write(hdr_str.format(*header))
            if self.sim_wells:
                self.__write_pkg_output(foo, "well")
            if self.sim_maw:
                self.__write_pkg_output(foo, "maw")
            if self.sim_diversions:
                self.__write_pkg_output(foo, "sfr")
            if self.sim_lak:
                self.__write_pkg_output(foo, "lak")

    def __write_pkg_output(self, fobj, pkg):
        """
        Generalized method to write all package outputs

        Parameters
        ----------
        fobj : FileObject
        pkg : str
            package string ("well", "sfr", "maw", "lak")
        """
        fmt_str = "{:>8} {:>8} {:>8} {:>10} {:8d} {:15f} {:15f} {:15f} {:15f} {:15f}\n"
        for _, output in getattr(self, f"{pkg}_output").items():
            for rec in output:
                fobj.write(fmt_str.format(*rec))

    @classmethod
    def load_output(cls, filename):
        """
        Method to load Agmvr output into a pandas dataframe

        Parameters:
        ----------
        filename : str
            Agmvr output file path
        """
        import pandas as pd

        df = pd.read_csv(filename, delim_whitespace=True)
        df["pid"] -= 1
        df["kstp"] -= 1
        return df

    def _save_variable_addresses(self, mf6):
        """
        Method that allows developers to dump out BMI variable addresses
        to a text file. Do not call directly, use the develop=True kwarg
        on run_moel() to dump variable addresses

        Parameters
        ----------
        mf6 : ModflowApi object

        """
        input_vars = mf6.get_input_var_names()
        output_vars = mf6.get_output_var_names()

        with open("input_vars.txt", "w") as foo:
            for i in input_vars:
                foo.write(f"{i}\n")

        with open("output_vars.txt", "w") as foo:
            for i in output_vars:
                foo.write(f"{i}\n")

    def run_model(self, dll=None, **kwargs):
        """
        Method to run MODFLOW6 with the MF6 API AG package

        Parameters
        ----------
        dll : str
            path to modflow API ".dll" or ".so"
        **kwargs : keyword arguments
            "develop" : bool
                flag to print out BMI input and output variable addresses for
                developer information about model data available through the
                BMI for a given model
        """
        develop = kwargs.pop("develop", False)

        if dll is None:
            import platform
            if platform.system().lower() == "linux":
                dll_name = "libmf6.so"
            elif platform.system().lower() == "darwin":
                dll_name = "libmf6.dylib"
            else:
                dll_name = "libmf6.dll"
            sws = os.path.abspath(os.path.dirname(__file__))
            dll = os.path.join(sws, "..", "bin", dll_name)

        mf6 = ModflowApi(
            dll,
            working_directory=self.sim.simulation_data.mfpath.get_sim_path(),
        )

        mf6.initialize()

        if develop:
            self._save_variable_addresses(mf6)

        self.create_addresses(mf6)
        prev_time = 0
        current_time = mf6.get_current_time()
        end_time = mf6.get_end_time()
        max_iter = mf6.get_value(self.maxiter_addr)

        kperold = 0
        kstp = 0
        while current_time < end_time:
            delt = current_time - prev_time
            dt = mf6.get_time_step()
            mf6.prepare_time_step(dt)
            kiter = 0

            kper = mf6.get_value(mf6.get_var_address("KPER", "TDIS"))[0] - 1
            if kper != kperold or current_time == 0.0:
                if current_time == 0:
                    pass
                else:
                    kperold += 1

                self.set_stress_period_data(mf6, kper)

            n_solutions = mf6.get_subcomponent_count()
            for sol_id in range(1, n_solutions + 1):
                mf6.prepare_solve(sol_id)

                if kiter == 0:
                    # need to zero out mvr values
                    self.zero_mvr(mf6)
                    mf6.solve(sol_id)

                while kiter < max_iter:
                    self.applied_irrigation = np.zeros(self.uzf_shape)

                    if self.sim_lak:
                        self.lak_demand_etdemand(mf6, kstp, delt, kiter)

                    if self.sim_diversions:
                        self.sfr_demand_etdemand(mf6, kstp, delt, kiter)

                    if self.sim_maw:
                        self.maw_demand_etdemand(mf6, kstp, delt, kiter)

                    if self.sim_wells:
                        self.well_demand_etdemand(mf6, kstp, delt, kiter)

                    has_converged = mf6.solve(sol_id)
                    kiter += 1
                    if has_converged:
                        break

                mf6.finalize_solve(sol_id)

            mf6.finalize_time_step()
            prev_time = current_time
            current_time = mf6.get_current_time()
            kstp += 1

            if not has_converged:
                print("model did not converge")

        try:
            mf6.finalize()
            success = True
        except:
            raise RuntimeError()

        self._write_output()
        return success


def plugin():
    return ModflowAgmvr


def dfn():
    dfn_name = "gwf-agmvr.dfn"
    return (
        dfn_name,
        os.path.join(os.path.abspath(os.path.dirname(__file__)), dfn_name)
    )