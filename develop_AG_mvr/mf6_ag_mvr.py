from modflowapi import ModflowApi
import numpy as np
import flopy


class ModflowAgmvr(object):
    """

    """
    def __init__(self, sim, ag_type, mvr_name):
        self.sim = sim
        name = list(sim.model_names)[0]
        self.name = name.upper()
        self.gwf = sim.get_model(name)
        self.mvr_name = mvr_name
        self.mvr = self.gwf.get_package(self.mvr_name)

        if ag_type.lower() in ('specified', 'trigger', 'etdemand'):
            self.ag_type = ag_type.lower()

        self.totim = np.add.accumulate(self.gwf.modeltime.perlen)

        self.sim_wells = False
        self.sim_diversions = False

        self.well_name = None
        self.sfr_name = None
        self.uzf_name = None

        self.get_pkgs()
        self.create_arrays()


    def get_pkgs(self):
        """
        Method to get AG package names

        """
        pkg_names = self.mvr.packages

        well = None
        sfr = None
        uzf = None
        for ix, name in enumerate(pkg_names.array):
            name = name[-1]
            pkg = self.gwf.get_package(name)
            if isinstance(pkg, flopy.mf6.ModflowGwfwel):
                well = name
                self.sim_wells = True
            elif isinstance(pkg, flopy.mf6.ModflowGwfsfr):
                sfr = name
                self.sim_diversions = True
            elif isinstance(pkg, flopy.mf6.ModflowGwfuzf):
                uzf = name
            else:
                pass

        if self.sim_wells:
            self.well_name = well
        if self.sim_diversions:
            self.sfr_name = sfr
        if uzf is None:
            raise AssertionError("UZF must be active to simulate AG")

        self.uzf_name = uzf
        self.dis_name = self.gwf.dis.name[0].upper()

    def create_arrays(self):
        """
        Method to dimensionalize arrays for AG calculations

        """
        # create well arrays
        period_data = self.mvr.perioddata.array
        if self.sim_wells:
            mxwlnum = 0
            for recarray in period_data:
                wlidx = np.where(recarray['pname1'] == self.well_name)[0]
                if len(wlidx) > 0:
                    wlids = recarray[wlidx]["id1"]
                    mxwls = np.max(wlids) + 1
                    if mxwls > mxwlnum:
                        mxwlnum = mxwls


            self.well_mq = np.zeros((mxwlnum,), dtype=int)
            self.well_active = np.zeros((mxwlnum,), dtype=int)
            self.sup = np.zeros((mxwlnum,))
            self.supold = np.zeros((mxwlnum,))
            self.well_num = np.arange(mxwlnum, dtype=int)

        # create diversion arrays
        if self.sim_diversions:
            mxdivnum = 0
            for recarray in period_data:
                dividx = np.where(recarray["pname1"] == self.sfr_name)[0]
                if len(dividx) > 0:
                    # todo: this is zero based....
                    divids = recarray[dividx]["id1"]
                    mxdiv = np.max(divids) + 1
                    if mxdiv > mxdivnum:
                        mxdivnum = mxdiv

            self.diversion_mq = np.zeros((mxdivnum,), dtype=int)
            self.diversion_active = np.zeros((mxdivnum,), dtype=int)
            self.idsup = np.zeros((mxdivnum,))
            self.idsupold = np.zeros((mxdivnum,))
            self.diversion_num = np.arange(mxdivnum, dtype=int)

    def create_addresses(self, mf6):
        """
        Method to create MODFLOW-API addresses

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

            self.wellq_for_mvr_addr = mf6.get_var_address(
                "QFORMVR", self.name, self.well_name.upper()
            )
            self.wellq_to_mvr_addr = mf6.get_var_address(
                "QTOMVR", self.name, self.well_name.upper()
            )

        if self.sim_diversions:
            self.sfrq_for_mvr_addr = mf6.get_var_address(
                "QFORMVR", self.name, self.sfr_name.upper()
            )
            self.sfrq_to_mvr_addr = mf6.get_var_address(
                "QTOMVR", self.name, self.sfr_name.upper()
            )
            self.sfrqt_for_mvr_addr = mf6.get_var_address(
                "QTFORMVR", self.name, self.sfr_name.upper()
            )
            self.sfrq_sim_to_mvr_addr = mf6.get_var_address(
                "SIMTOMVR", self.name, self.sfr_name.upper()
            )
            self.sfr_qextoutflow_addr = mf6.get_var_address(
                "QEXTOUTFLOW", self.name, self.sfr_name.upper()
            )

    def set_stress_period_data(self, mf6, kper):
        """
        Method to set stress period data from MVR (MAXQ, Delivery cells, etc)

        """
        if self.sim_wells:
            qformvr = mf6.get_value(self.wellq_for_mvr_addr)
            qtomvr = mf6.get_value(self.wellq_to_mvr_addr)
            max_q = np.min(np.vstack((qformvr, qtomvr)), axis=0)
            recarray = self.mvr.perioddata.data[kper]
            wlidx = np.where(recarray['pname1'] == self.well_name)[0]
            if len(wlidx) > 0:
                wellids = sorted(np.unique(recarray[wlidx]["id1"]))
            else:
                wellids = []

            active = np.zeros(max_q.shape, dtype=int)
            active[wellids] = 1

            irrigated_cells = []
            irrigated_proportion = []
            mvr_index = []
            for wlid in range(max_q.size):
                icells = []
                iprop = []
                idx = np.where(recarray["id1"] == wlid)[0]
                if len(idx) > 0:
                    icells = recarray[idx]["id2"]
                    iprop = recarray[idx]["value"] / np.sum(recarray[idx]["value"])


                irrigated_cells.append(icells)
                irrigated_proportion.append(iprop)
                mvr_index.append(idx)

            self.well_active = active
            self.well_maxq = max_q
            self.well_irrigated_cells = irrigated_cells
            self.well_irrigated_proportion = irrigated_proportion
            self.well_mvr_index = mvr_index
            self.sup = np.zeros(max_q.shape)
            self.supold = np.zeros(max_q.shape)
            self.aetold = np.zeros(max_q.shape)

        if self.sim_diversions:
            qformvr = mf6.get_value(self.sfrq_for_mvr_addr)
            qtomvr = mf6.get_value(self.sfrq_to_mvr_addr)
            max_q = np.min(np.vstack((qformvr, qtomvr)), axis=0)
            recarray = self.mvr.perioddata.data[kper]
            sfridx = np.where(recarray["pname1"] == self.sfr_name)[0]
            if len(sfridx) > 0:
                sfrids = sorted(np.unique(recarray[sfridx]["id1"]))
            else:
                sfrids = []

            active = np.zeros(max_q.shape, dtype=int)
            active[sfrids] = 1

            irrigated_cells = []
            irrigated_proportion = []
            mvr_index = []
            for sfrid in range(max_q.size):
                icells = []
                iprop = []
                idx = np.where(recarray["id1"] == sfrid)[0]
                if len(idx) > 0:
                    icells = recarray[idx]["id2"]
                    iprop = recarray[idx]["value"] / np.sum(recarray[idx]["value"])

                irrigated_cells.append(icells)
                irrigated_proportion.append(iprop)
                mvr_index.append(idx)

            self.sfr_active = active
            self.sfr_maxq = max_q
            self.sfr_irrigated_cells = irrigated_cells
            self.sfr_irrigated_proportion = irrigated_proportion
            self.sfr_mvr_index = mvr_index
            self.idsup = np.zeros(max_q.shape)
            self.idsupold = np.zeros(max_q.shape)
            self.idaetold = np.zeros(max_q.shape)

    def gw_demand_etdemand(self, mf6, kstp, delt=1, kiter=1):
        """
        Method to determine groundwater use demand

        Parameters
        ----------
        mf6 : ModflowApi object
        delt : float
        kiter : int
        """
        pet = mf6.get_value(self.pet_addr)
        vks = mf6.get_value(self.vks_addr)
        area = mf6.get_value(self.area_addr)
        aet = mf6.get_value(self.aet_addr)

        crop_vks = np.zeros(self.well_maxq.shape)
        crop_pet = np.zeros(self.well_maxq.shape)
        crop_aet = np.zeros(self.well_maxq.shape)
        sfr_applied = np.zeros(self.well_maxq.shape)
        for ix, crop_nodes in enumerate(self.well_irrigated_cells):
            if len(crop_nodes) > 0:
                crop_vks[ix] = np.sum(vks[crop_nodes] * area[crop_nodes])
                crop_pet[ix] = np.sum(pet[crop_nodes] * area[crop_nodes])
                crop_aet[ix] = np.sum(aet[crop_nodes] * area[crop_nodes])
                sfr_applied[ix] = np.sum(self.applied_irrigation[crop_nodes])

        crop_aet = np.where(crop_vks < 1e-30, crop_pet, crop_aet)

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            sup = np.zeros(self.well_maxq.shape)
            if kstp == 0:
                supold = np.zeros(self.well_maxq.shape)
            else:
                supold = self.supold
        else:
            sup = self.sup
            supold = self.supold

        factor = self._calculate_factor(crop_pet, crop_aet, self.aetold, sup, supold, kiter)

        qonly = np.where(sup + factor > crop_vks, crop_vks, sup + factor)

        self.supold = sup
        self.aetold = crop_aet

        if self.sim_diversions:
            qonly = qonly - sfr_applied
            qonly = np.where(qonly < 0, 0, qonly)

        pumping = np.where(qonly > self.well_maxq, self.well_maxq, -1 * qonly)
        pumping = np.where(np.abs(pumping) < 1e-10, 0, pumping)

        self.sup = np.abs(pumping)

        active_ix = np.where(self.well_active)[0]

        if len(active_ix) > 0:
            wells = mf6.get_value(self.well_addr)
            wells[active_ix, 0] = pumping[active_ix]
            mf6.set_value(self.well_addr, wells)

            mvr = mf6.get_value(self.mvr_value_addr)
            for well in active_ix:
                idx = self.well_mvr_index[well]
                mvr[idx] = pumping[well] * self.well_irrigated_proportion[well]
                self.applied_irrigation[self.well_irrigated_cells[well]] = \
                    pumping[well] * self.well_irrigated_proportion[well]

            mf6.set_value(self.mvr_value_addr, mvr)

    def sw_demand_etdemand(self, mf6, kstp, delt=1, kiter=1):
        """
        Method to determine surface-water demand

        Parameters
        ----------
        mf6 : ModflowApi object
        delt : float
        kiter : int

        """
        qtomvr = mf6.get_value(self.sfrq_to_mvr_addr)
        pet = mf6.get_value(self.pet_addr)
        vks = mf6.get_value(self.vks_addr)
        area = mf6.get_value(self.area_addr)
        aet = mf6.get_value(self.aet_addr)

        crop_vks = np.zeros(self.sfr_maxq.shape)
        crop_pet = np.zeros(self.sfr_maxq.shape)
        crop_aet = np.zeros(self.sfr_maxq.shape)
        for ix, crop_nodes in enumerate(self.sfr_irrigated_cells):
            if len(crop_nodes) > 0:
                crop_vks[ix] = np.sum(vks[crop_nodes] * area[crop_nodes])
                crop_pet[ix] = np.sum(pet[crop_nodes] * area[crop_nodes])
                crop_aet[ix] = np.sum(aet[crop_nodes] * area[crop_nodes])

        crop_aet = np.where(crop_vks < 1e-30, crop_pet, crop_aet)

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            self.idsupold = np.zeros(self.sfr_maxq.shape)
            qtomvr[:] = 0
            self.idaetold = crop_aet
            if kstp == 0:
                self.idaetold[:] = 0

        sup = qtomvr
        supold = self.idsupold
        factor = self._calculate_factor(crop_pet, crop_aet, self.idaetold, sup, supold, kiter)

        qonly = np.where(sup + factor > crop_vks, crop_vks, sup + factor)
        factor = np.where(factor < 0, 0, factor)

        self.idsupold = qtomvr
        self.idsup += factor
        self.idaetold = crop_aet

        dvflw = np.where(qonly >= self.sfr_maxq, self.sfr_maxq, qonly)

        active_ix = np.where(self.sfr_active)[0]

        if len(active_ix) > 0:
            diversions = mf6.get_value(self.mvr_value_addr)
            for seg in active_ix:
                idx = self.sfr_mvr_index[seg]
                diversions[idx] = dvflw[seg] * self.sfr_irrigated_proportion[seg]
                self.applied_irrigation[self.sfr_irrigated_cells[seg]] = \
                    dvflw[seg] * self.sfr_irrigated_proportion[seg]
            mf6.set_value(self.mvr_value_addr, diversions)

    def _calculate_factor(self, crop_pet, crop_aet, aetold, sup, supold, kiter, accel=1):
        """
        Method to calculate unmet demand for crops

        Parameters
        ----------
        crop_pet :
        crop_aet :
        aetold :
        sup :
        supold :
        kiter : float
        accel : float
        """
        et_diff = crop_pet - crop_aet
        factor = et_diff
        det = crop_aet - aetold
        dq = sup - supold

        if kiter > 0:
            factor = np.where(np.abs(det) > 1e-05,
                              dq * et_diff / det,
                              factor)

        factor = np.where(factor > et_diff * accel,
                          et_diff * accel,
                          factor)
        factor = np.where(factor < et_diff, et_diff, factor)
        factor = np.where(factor < 1e-05, 1e-05, factor)
        return factor

    def run_model(self, dll):
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

            n_solutions = mf6.get_subcomponent_count()
            for sol_id in range(1, n_solutions + 1):
                mf6.prepare_solve(sol_id)
                mf6.solve(sol_id)

                self.applied_irrigation = np.zeros(self.uzf_shape)
                # updated stress period information
                if current_time in self.totim or current_time == 0.:
                    if current_time == 0:
                        pass
                    else:
                        kper += 1

                    self.set_stress_period_data(mf6, kper)

                while kiter < max_iter:
                    if self.sim_diversions:
                        self.sw_demand_etdemand(mf6, kstp, delt, kiter)

                    if self.sim_wells:
                        self.gw_demand_etdemand(mf6, kstp, delt, kiter)

                    has_converged = mf6.solve(sol_id)
                    uzf2 = mf6.get_value(self.uzf_mvr_addr)
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