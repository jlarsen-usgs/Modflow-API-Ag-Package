from modflowapi import ModflowApi
import numpy as np
import flopy


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
    def __init__(self, sim, ag_type, mvr_name):
        self.sim = sim
        name = list(sim.model_names)[0]
        self.name = name.upper()
        self.gwf = sim.get_model(name)
        self.mvr_name = mvr_name
        self.mvr = self.gwf.get_package(self.mvr_name)
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
                            rec["value"]
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

        if ag_type.lower() in ('specified', 'trigger', 'etdemand'):
            self.ag_type = ag_type.lower()

        self.totim = np.add.accumulate(self.gwf.modeltime.perlen)

        self.sim_wells = False
        self.sim_maw = False
        self.sim_diversions = False

        self.well_name = None
        self.maw_name = None
        self.sfr_name = None
        self.uzf_name = None

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
        for ix, name in enumerate(pkg_names.array):
            name = name[-1]
            pkg = self.gwf.get_package(name)
            if isinstance(pkg, flopy.mf6.ModflowGwfwel):
                well = name
                self.sim_wells = True
            elif isinstance(pkg, flopy.mf6.ModflowGwfmaw):
                maw = name
                self.sim_maw = True
            elif isinstance(pkg, flopy.mf6.ModflowGwfsfr):
                sfr = name
                self.sim_diversions = True
            elif isinstance(pkg, flopy.mf6.ModflowGwfuzf):
                uzf = name
            else:
                pass

        if self.sim_wells:
            self.well_name = well
        if self.sim_maw:
            self.maw_name = maw
        if self.sim_diversions:
            self.sfr_name = sfr
        if uzf is None:
            raise AssertionError("UZF must be active to simulate AG")

        self.uzf_name = uzf
        self.dis_name = self.gwf.dis.name[0].upper()

    def create_addresses(self, mf6):
        """
        Method to create MODFLOW-API addresses

        Parameters
        ----------
        mf6 : ModflowApi object

        """
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
            self.well_addr = mf6.get_var_address(
                "BOUND", self.name, self.well_name.upper()
            )

        if self.sim_maw:
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
            max_q = np.zeros(qformvr.shape)

        mvr = mf6.get_value(self.mvr_value_addr)
        pkg_name = self.__dict__[f"{pkg}_name"]
        recarray = self.mvr.perioddata.data[kper]
        irridx = np.where(recarray['pname1'] == pkg_name)[0]
        if len(irridx) > 0:
            if pkg_name in ("sfr", ):
                mvr[irridx] = 0

            irrids = sorted(np.unique(recarray[irridx]["id1"]))
        else:
            irrids = []

        mf6.set_value(self.mvr_value_addr, mvr)

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
            idx = np.where((recarray["id1"] == irrid) & (recarray["pname1"] == pkg_name))[0]
            if len(idx) > 0:
                icells = recarray[idx]["id2"]
                iprop = recarray[idx]["value"] / np.sum(recarray[idx]["value"])
                if pkg in ("sfr",):
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

        setattr(self, f"{pkg}_active", active)
        setattr(self, f"{pkg}_maxq", max_q)
        setattr(self, f"{pkg}_irrigated_cells", irrigated_cells)
        setattr(self, f"{pkg}_irrigated_proportion", irrigated_proportion)
        setattr(self, f"{pkg}_irrigation_efficiency", irrigation_efficiency)
        setattr(self, f"{pkg}_application_fraction", application_fraction)
        setattr(self, f"{pkg}_mvr_index", mvr_index)

        return max_q.shape

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
            mq_shape = self._set_package_stress_period_data(mf6, kper, 'well')
            self.wellsup = np.zeros(mq_shape)
            self.wellsupold = np.zeros(mq_shape)
            self.aetold = np.zeros(mq_shape)

        if self.sim_maw:
            mq_shape = self._set_package_stress_period_data(mf6, kper, "maw")
            self.mawsup = np.zeros(mq_shape)
            self.mawsupold = np.zeros(mq_shape)
            self.mawaetold = np.zeros(mq_shape)

        if self.sim_diversions:
            mq_shape = self._set_package_stress_period_data(mf6, kper, "sfr")
            self.sfrsup = np.zeros(mq_shape)
            self.sfrsupold = np.zeros(mq_shape)
            self.idaetold = np.zeros(mq_shape)
            self.evap_old = mf6.get_value(self.sfr_evap_addr)

    def zero_mvr(self, mf6):
        """
        Method to zero out MVR values before initial solve
        """
        mvr = mf6.get_value(self.mvr_value_addr)
        mvr[:] = 0
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
        area = mf6.get_value(self.area_addr)
        aet = mf6.get_value(self.aet_addr)
        gwet = mf6.get_value(self.gwet_addr)
        maxq = getattr(self, f"{pkg}_maxq")
        application_fraction = getattr(self, f"{pkg}_application_fraction")
        irrigated_cells = getattr(self, f"{pkg}_irrigated_cells")

        crop_vks = np.zeros(maxq.shape)
        crop_pet = np.zeros(maxq.shape)
        crop_aet = np.zeros(maxq.shape)
        crop_gwet = np.zeros(maxq.shape)
        app_frac = np.zeros(maxq.shape)
        prev_applied = np.zeros(maxq.shape)
        for ix, crop_nodes in enumerate(irrigated_cells):
            if len(crop_nodes) > 0:
                crop_vks[ix] = np.sum(vks[crop_nodes] * area[crop_nodes])
                crop_pet[ix] = np.sum(pet[crop_nodes] * area[crop_nodes])
                crop_aet[ix] = np.sum(aet[crop_nodes] * area[crop_nodes])
                crop_gwet[ix] = np.sum(gwet[crop_nodes])
                app_frac[ix] = np.mean(application_fraction[ix])
                prev_applied[ix] = np.sum(self.applied_irrigation[crop_nodes])

        crop_aet += crop_gwet
        if pkg in ("well", "maw"):
            crop_aet = np.where(crop_vks < 1e-30, crop_pet, crop_aet)

        return crop_pet, crop_aet, crop_vks, app_frac, prev_applied

    def gw_demand_etdemand(self, mf6, kstp, delt=1, kiter=1):
        """
        Method to determine groundwater use demand

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
        crop_pet, crop_aet, crop_vks, app_frac, sfr_applied = self._set_etdemand_variables(mf6, "well")

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            sup = np.zeros(self.well_maxq.shape)
            if kstp == 0:
                supold = np.zeros(self.well_maxq.shape)
            else:
                supold = self.wellsupold
        else:
            sup = self.wellsup
            supold = self.wellsupold

        factor = self._calculate_factor(crop_pet, crop_aet, self.aetold, sup, supold, kiter)
        factor *= app_frac

        qonly = np.where(sup + factor > crop_vks, crop_vks, sup + factor)

        self.wellsupold = sup
        self.aetold = crop_aet

        if self.sim_diversions:
            qonly = qonly - sfr_applied
            qonly = np.where(qonly < 0, 0, qonly)

        pumping = np.where(qonly > self.well_maxq, -1 * self.well_maxq, -1 * qonly)
        pumping = np.where(np.abs(pumping) <= 1e-10, 0, pumping)

        self.wellsup = np.abs(pumping)

        active_ix = np.where(self.well_active)[0]

        if len(active_ix) > 0:
            wells = mf6.get_value(self.well_addr)
            wells[active_ix, 0] = pumping[active_ix]
            mf6.set_value(self.well_addr, wells)

            mvr = mf6.get_value(self.mvr_value_addr)
            for well in active_ix:
                idx = self.well_mvr_index[well]
                app_frac_proportion = (self.well_application_fraction[well] / np.sum(self.well_application_fraction[well])) / (1 / len(idx))
                mvr[idx] = (np.abs(pumping[well]) * self.well_irrigated_proportion[well]) * app_frac_proportion * self.well_irrigation_efficiency[well]
                self.applied_irrigation[self.well_irrigated_cells[well]] = \
                    (np.abs(pumping[well]) * self.well_irrigated_proportion[well]) * app_frac_proportion * self.well_irrigation_efficiency[well]

            mf6.set_value(self.mvr_value_addr, mvr)

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
        crop_pet, crop_aet, crop_vks, app_frac, sfr_applied = self._set_etdemand_variables(mf6, "maw")

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            sup = np.zeros(self.maw_maxq.shape)
            if kstp == 0:
                supold = np.zeros(self.maw_maxq.shape)
            else:
                supold = self.mawsupold
        else:
            sup = self.mawsup
            supold = self.mawsupold

        factor = self._calculate_factor(crop_pet, crop_aet, self.mawaetold, sup, supold, kiter)
        factor *= app_frac

        qonly = np.where(sup + factor > crop_vks, crop_vks, sup + factor)

        self.mawsupold = sup
        self.mawaetold = crop_aet

        if self.sim_diversions:
            qonly = qonly - sfr_applied
            qonly = np.where(qonly < 0, 0, qonly)

        pumping = np.where(qonly > self.maw_maxq, -1 * self.maw_maxq, -1 * qonly)
        pumping = np.where(np.abs(pumping) <= 1e-10, 0, pumping)

        self.mawsup = np.abs(pumping)

        active_ix = np.where(self.maw_active)[0]

        if len(active_ix) > 0:
            maws = mf6.get_value(self.maw_addr)
            maws[active_ix] = pumping[active_ix]
            mf6.set_value(self.maw_addr, maws)

            mvr = mf6.get_value(self.mvr_value_addr)
            for maw in active_ix:
                idx = self.maw_mvr_index[maw]
                app_frac_proportion = (self.maw_application_fraction[maw] / np.sum(self.maw_application_fraction[maw])) / (1 / len(idx))
                mvr[idx] = (np.abs(pumping[maw]) * self.maw_irrigated_proportion[maw]) * app_frac_proportion * self.maw_irrigation_efficiency[maw]
                self.applied_irrigation[self.maw_irrigated_cells[maw]] = \
                   (np.abs(pumping[maw]) * self.maw_irrigated_proportion[maw]) * app_frac_proportion * self.maw_irrigation_efficiency[maw]

            mf6.set_value(self.mvr_value_addr, mvr)

    def sw_demand_etdemand(self, mf6, kstp, delt=1, kiter=1):
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
        """
        qtomvr = mf6.get_value(self.sfrq_to_mvr_addr)
        crop_pet, crop_aet, crop_vks, app_frac, prev_applied = self._set_etdemand_variables(mf6, "sfr")

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            self.sfrsupold = np.zeros(self.sfr_maxq.shape)
            self.sfrsup = np.zeros(self.sfr_maxq.shape)
            qtomvr[:] = 0
            self.idaetold[:] = crop_aet[:]
            if kstp == 0:
                self.idaetold[:] = 0

        sup = qtomvr
        supold = self.sfrsupold
        factor = self._calculate_factor(crop_pet, crop_aet, self.idaetold, sup, supold, kiter)
        factor *= app_frac

        qonly = np.where(sup + factor > crop_vks, crop_vks, sup + factor)
        factor = np.where(factor < 0, 0, factor)

        self.sfrsupold[:] = qtomvr[:]
        self.sfrsup += factor
        self.idaetold = crop_aet

        dvflw = np.where(qonly >= self.sfr_maxq, self.sfr_maxq, qonly)

        active_ix = np.where(self.sfr_active)[0]

        if len(active_ix) > 0:
            diversions = mf6.get_value(self.mvr_value_addr)
            length = mf6.get_value(self.sfr_length_addr)
            width = mf6.get_value(self.sfr_width_addr)
            evap = self.evap_old.copy()
            for seg in active_ix:
                idx = self.sfr_mvr_index[seg]
                app_frac_proportion = (self.sfr_application_fraction[seg] / np.sum(self.sfr_application_fraction[seg])) / (1 / len(idx))
                div_requested = (dvflw[seg] * self.sfr_irrigated_proportion[seg]) * app_frac_proportion
                div_inefficient = div_requested * self.sfr_irrigation_efficiency[seg]
                diversions[idx] = div_inefficient
                evap[seg] += np.sum(div_requested - div_inefficient) / (length[seg] * width[seg])
                self.applied_irrigation[self.sfr_irrigated_cells[seg]] = \
                    (dvflw[seg] * self.sfr_irrigated_proportion[seg]) * app_frac_proportion
            mf6.set_value(self.mvr_value_addr, diversions)
            mf6.set_value(self.sfr_evap_addr, evap)

    def _calculate_factor(self, crop_pet, crop_aet, aetold, sup, supold, kiter, accel=1):
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
            # original nwt equation logic commented out. I don't think
            # that we should continuously apply et_diff if det is small
            # instead there should be zero change, this prevents runaway
            # pumping from the equation...
            # factor = np.where(np.abs(det) > 1e-05,
            #                   dq * et_diff / det,
            #                   factor)
            factor = np.where(np.abs(det) > 1e-6,
                              dq * et_diff / det,
                              0)

        factor = np.where(factor > et_diff * accel,
                          et_diff * accel,
                          factor)
        # original equation logic commented out. I believe this to be one
        # of the culprits for causing runaway pumping values!
        # factor = np.where(factor < et_diff, et_diff, factor)
        factor = np.where(factor < 1e-05, 1e-05, factor)
        return factor

    def run_model(self, dll):
        """
        Method to run MODFLOW6 with the MF6 API AG package

        Parameters
        ----------
        dll : str
            path to modflow API ".dll" or ".so"
        """
        mf6 = ModflowApi(
            dll,
            working_directory=self.sim.simulation_data.mfpath.get_sim_path()
        )

        mf6.initialize()

        self.create_addresses(mf6)
        prev_time = 0
        current_time = mf6.get_current_time()
        end_time = mf6.get_end_time()

        input_vars = mf6.get_input_var_names()
        output_vars = mf6.get_output_var_names()

        max_iter = mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))

        with open("input_vars.txt", "w") as foo:
            for i in input_vars:
                foo.write(f"{i}\n")

        with open("output_vars.txt", "w") as foo:
            for i in output_vars:
                foo.write(f"{i}\n")

        kper = 0
        kstp = 0
        while current_time < end_time:
            delt = current_time - prev_time
            dt = mf6.get_time_step()
            mf6.prepare_time_step(dt)
            kiter = 0

            if current_time in self.totim or current_time == 0.:
                if current_time == 0:
                    pass
                else:
                    kper += 1

                self.set_stress_period_data(mf6, kper)

            n_solutions = mf6.get_subcomponent_count()
            for sol_id in range(1, n_solutions + 1):
                mf6.prepare_solve(sol_id)

                if kiter == 0:
                    # need to zero out mvr values
                    self.zero_mvr(mf6)
                    mf6.solve(sol_id)

                self.applied_irrigation = np.zeros(self.uzf_shape)

                while kiter < max_iter:
                    if self.sim_diversions:
                        self.sw_demand_etdemand(mf6, kstp, delt, kiter)

                    if self.sim_maw:
                        self.maw_demand_etdemand(mf6, kstp, delt, kiter)

                    if self.sim_wells:
                        self.gw_demand_etdemand(mf6, kstp, delt, kiter)

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

        return success
