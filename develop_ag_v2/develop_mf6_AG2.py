from modflowapi import ModflowApi
import numpy as np


class Modflow6Ag(object):
    """

    """
    def __init__(self, sim, ag):

        self.sim = sim
        self.ag = ag
        name = list(sim.model_names)[0]
        self.gwf = sim.get_model(name)
        self.name = name.upper()
        self.ncpl = self.gwf.modelgrid.ncpl
        self.nlay = self.gwf.modelgrid.nlay
        self.totim = np.add.accumulate(self.gwf.modeltime.perlen)
        self.area = None
        self.top = None
        self.top_tile = None
        self.botm = None
        self.nsteps = np.sum(self.gwf.modeltime.nstp)

        self.etdemand = self.ag.etdemand
        self.trigger = self.ag.trigger

        self.sim_wells = self.ag.irrigation_well
        if self.ag.irrwell is None:
            self.sim_wells = False

        self.nummaxwell = self.ag.nummaxwell

        self.sim_diversions = self.ag.irrigation_diversion
        self.numirrdiversion = self.ag.numirrdiversions

        self.sim_supplemental = self.ag.supplemental_well
        self.numsupwell = self.ag.numsupwells
        self.irrigation_return = []

        if self.sim_wells or self.sim_supplemental:
            self.well_max_q = np.zeros((self.nummaxwell,))
            self.well_out = np.zeros((self.nsteps, self.nummaxwell))
            self.well_demand_out = np.zeros((self.nsteps, self.nummaxwell))
            self._sup = np.zeros((self.nummaxwell,))
            self._supold = np.zeros((self.nummaxwell,))
            self.irrwell_num = None
            self.irrwell = None
            self.well_triggerfact = None

        if self.sim_diversions:
            self.div_info = self.gwf.sfr.diversions.array
            self.div_info.sort(order='idv', axis=0)
            self.sfr_max_q = np.zeros((self.numirrdiversion,))
            self.sfr_out = np.zeros((self.nsteps, self.numirrdiversion))
            self.sfr_demand_out = np.zeros((self.nsteps, self.numirrdiversion))
            self._idsup = np.zeros((self.numirrdiversion,))
            self._idsupold = np.zeros((self.numirrdiversion,))

            self.sfr_irrtrigger = None
            self.sfr_irrtriggerflg = None
            self.irrdiversion_num = None
            self.irrdiversion = None

        if self.sim_supplemental:
            self.supwell_num = None
            self.supwell = None

        # set these here as placeholders to avoid if statements in run loop.
        # todo: these need to be changed to numirrwell size/numirrdiversion size
        self.well_irrperiod = np.zeros((1,))
        self.well_timeinperiod = np.ones((1,)) * 1e+30

        self.sfr_irrperiod = np.zeros((1,), )
        self.sfr_timeinperiod = np.ones((1,)) * 1e+30


    def create_addresses(self, mf6):
        sto_name = self.gwf.sto.name[0].upper()
        dis_name = self.gwf.dis.name[0].upper()
        # try:
        #     rch_name = self.gwf.rcha.name[0].upper()
        # except AttributeError:
        #     rch_name = self.gwf.rch.name[0].upper()

        # try:
        #     evt_name = self.gwf.evta.name[0].upper()
        # except AttributeError:
        #     evt_name = self.gwf.evt.name[0].upper()

        try:
            uzf_name = self.gwf.uzf.name[0].upper()
        except AttributeError:
            uzf_name = None

        try:
            sfr_name = self.gwf.sfr.name[0].upper()
        except AttributeError:
            sfr_name = None

        self.area = mf6.get_value(
            mf6.get_var_address("AREA", self.name, dis_name)
        )

        self.sy = mf6.get_value(
            mf6.get_var_address("SY", self.name, sto_name)
        )
        self.ss = mf6.get_value(
            mf6.get_var_address("SS", self.name, sto_name)
        )
        self.botm = mf6.get_value(
            mf6.get_var_address("BOT", self.name, dis_name)
        )

        self.top = mf6.get_value(
            mf6.get_var_address("TOP", self.name, dis_name)
        )

        self.top_tile = np.tile(self.top, self.nlay)

        # todo: set the well name!!!!
        self.iusesy_addr = mf6.get_var_address("IUSESY", self.name, sto_name)
        self.well_addr = mf6.get_var_address("BOUND", self.name, "WEL_0")
        self.well_node_addr = mf6.get_var_address(
            "NODELIST", self.name, "WEL_0"
        )
        # self.rch_addr = mf6.get_var_address("BOUND", self.name, rch_name)
        # self.evt_addr = mf6.get_var_address("BOUND", self.name, evt_name)
        self.head_addr = mf6.get_var_address("X", self.name)
        if uzf_name is not None:
            self.area_addr = mf6.get_var_address("UZFAREA", self.name, uzf_name)
            self.vks_addr = mf6.get_var_address("VKS", self.name, uzf_name)
            self.pet_addr = mf6.get_var_address("PET", self.name, uzf_name)
            self.uzet_addr = mf6.get_var_address("UZET", self.name, uzf_name)
            self.uz_gwet_addr = mf6.get_var_address("GWET", self.name, uzf_name)
            self.aet_addr = mf6.get_var_address("ETACT", self.name, uzf_name)
            self.uz_rhs_addr = mf6.get_var_address("RHS", self.name, uzf_name)
            self.uz_finf_addr = mf6.get_var_address("FINF", self.name, uzf_name)
            self.uz_sinf_addr = mf6.get_var_address("SINF", self.name, uzf_name)
            self.uz_hcof_addr = mf6.get_var_address("HCOF", self.name, uzf_name)
        if sfr_name is not None:
            self.sfr_divflow_addr = mf6.get_var_address("DIVFLOW", self.name, sfr_name)
            self.sfr_dsflow_addr = mf6.get_var_address("DSFLOW", self.name, sfr_name)

    def set_max_q_well(self, mf6):
        """
        Method to set the max-q of AG wells

        Parameters
        ----------
        mf6 : ModflowApi object
        """
        well = mf6.get_value(self.well_addr)
        self.well_max_q = np.copy(well.T[0])
        # change max_q to zero and reset based on AG-Demand
        well[:, 0] = 0
        mf6.set_value(self.well_addr, well)

    def set_max_q_diversion(self, mf6):
        """
        Method to set the max-q of AG diversions

        Parameters
        ----------
        mf6 : ModflowApi object
        """
        sfr = mf6.get_value(self.sfr_divflow_addr)
        self.sfr_max_q = np.copy(sfr)

    def __fmt_irr_spd(self, spd, kper):
        """
        Method to get irrwell and irrdiversion data into an easy
        to use format

        Parameters
        ----------
        spd : dict
            dict of irrwell or irrdiversion np.record arrays
        kper : int
            stress period number (zero based)

        """
        irrwell_info = []
        for rec in spd[kper]:
            num = rec['numcell']
            t = [
                (self.gwf.modelgrid.get_node(
                    (0, rec[f"i{ix}"], rec[f"j{ix}"])
                )[0],
                 rec[f"eff_fact{ix}"],
                 rec[f"field_fact{ix}"]
                 )
                for ix in range(num)
            ]
            irrwell_info.append(np.array(t, dtype=[("node", int),
                                                   ("eff_fact", float),
                                                   ("field_fact", float)]))
        return irrwell_info

    def __fmt_supwell_spd(self, spd, kper):
        """
        Method reformat supwell data into an easy to use format

        Parameters
        ----------
        spd : dict
            dict of supwell np.record arrays
        kper : int
            stress period number (zero based)

        """
        supwell_info = []
        for rec in spd[kper]:
            num = rec['numcell']
            t = [(rec[f"segid{ix}"] - 1,
                  rec[f"fracsup{ix}"],
                  rec[f"fracsupmax{ix}"])
                 for ix in range(num)
            ]
            supwell_info.append(np.array(t, dtype=[("segid", int),
                                                   ("fracsup", float),
                                                   ("fracsupmax", float)]))
        return supwell_info

    def set_irrwell_stress_period_data(self, kper):
        """
        Method to set stress period data from wells

        Parameters:
        ----------
        kper : int
            stress period number (zero based)
        """
        self.irrwell_num = self.ag.irrwell[kper]['wellid']
        self.well_irrperiod = self.ag.irrwell[kper]['period']
        self.well_triggerfact = self.ag.irrwell[kper]['triggerfact']
        if kper == 0:
            # todo: change this to length of well list!
            self.well_timeinperiod = np.ones(self.well_irrperiod.size) * 1e+30
        self.irrwell = self.__fmt_irr_spd(self.ag.irrwell, kper)

    def set_irrdiversion_stress_period_data(self, kper):
        """
        Method to set stress period data from diversions

        Parameters:
        ----------
        kper : int
            stress period number (zero based)
        """
        self.irrdiversion_num = self.ag.irrdiversion[kper]['segid'] - 1
        self.sfr_irrperiod = self.ag.irrdiversion[kper]['period']
        self.sfr_triggerfact = self.ag.irrdiversion[kper]['triggerfact']
        if kper == 0:
            # todo: this needs to be changed to length of diversion list!
            self.sfr_timeinperiod = np.ones(self.sfr_irrperiod.size) * 1e+30
        self.irrdiversion = self.__fmt_irr_spd(self.ag.irrdiversion, kper)

    def set_supwell_stress_period_data(self, kper):
        """
        Method to set stress period data from supplemental wells

        Parameters:
        ----------
        kper : int
            stress period number (zero based)
        """
        self.supwell_num = self.ag.supwell[kper]['wellid']
        self.supwell = self.__fmt_supwell_spd(self.ag.supwell, kper)

    def conjuctive_demand_etdemand(self, mf6, kstp, delt=1, kiter=1, kper=0):
        """
        Method to determine conjunctive use demand

        Parameters
        ----------
        mf6 : ModflowApi object
        delt : float
        kiter : int
        """
        pet = mf6.get_value(self.pet_addr)
        vks = mf6.get_value(self.vks_addr)
        area = mf6.get_value(self.area_addr)
        divflow = mf6.get_value(self.sfr_divflow_addr)
        dsflow = mf6.get_value(self.sfr_dsflow_addr)
        aet = mf6.get_value(self.aet_addr)

        crop_vks = np.zeros((len(self.irrdiversion_num,)))
        crop_pet = np.zeros((len(self.irrdiversion_num,)))
        crop_aet = np.zeros((len(self.irrdiversion_num,)))
        out_factor = np.zeros(self.sfr_max_q.shape)
        for ix, record in enumerate(self.irrdiversion):
            crop_nodes = record["node"]
            crop_vks[ix] = np.sum(vks[crop_nodes] * area[crop_nodes])
            crop_pet[ix] = np.sum(pet[crop_nodes] * area[crop_nodes])
            crop_aet[ix] = np.sum(aet[crop_nodes] * area[crop_nodes])

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            self._supact = np.zeros(self.sfr_max_q.shape)
            self._idsupold = np.zeros(self.sfr_max_q.shape)
            self._aetoldsw = crop_aet
            divflow[self.irrdiversion_num] = 0
            if kstp == 0:
                self._aetoldsw = crop_aet * 0

        sup = divflow[self.irrdiversion_num] # + actual???? must be supplementary pumping.....
        supold = self._idsupold[self.irrdiversion_num] # + actualold??? supplementary pumping....

        factor = self.__calculate_factor(
            crop_pet,
            crop_aet,
            self._aetoldsw,
            sup,
            supold,
            kiter,
            1
        )

        qonly = np.where(sup + factor > crop_vks, crop_vks, sup + factor)
        factor = np.where(factor < 0, 0, factor)

        self._idsupold[self.irrdiversion_num] = divflow[self.irrdiversion_num]
        # self._idsupold[self.irrdiversion_num] = idsflow[self.irrdiversion_num]
        self._supact[self.irrdiversion_num] += factor
        self._aetold = crop_aet

        with open('etdiv_factor.txt', "a") as foo:
            foo.write(f"{kper},{kstp},{crop_pet[0]},{crop_aet[0]},{self._supact[0]},{self._idsupold[0]},{factor[0]}\n")

        dvflw = np.where(qonly >= self.sfr_max_q,
                         self.sfr_max_q,
                         qonly)

        div_info = self.div_info[self.irrdiversion_num]
        fmaxflow = dsflow[div_info["rno"]]
        # fmaxflow = divflow[self.irrdiversion_num]

        dvflw = np.where(dvflw > fmaxflow, fmaxflow, dvflw)
        divflow[self.irrdiversion_num] = dvflw

        mf6.set_value(self.sfr_divflow_addr, divflow)

        # for supplemental pumping this needs to be dimensioned as sfr ndiv
        out_factor[self.irrdiversion_num] = factor

        return out_factor, divflow

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

        crop_vks = np.zeros((len(self.irrwell_num,)))
        crop_pet = np.zeros((len(self.irrwell_num,)))
        crop_aet = np.zeros((len(self.irrwell_num)))
        for ix, record in enumerate(self.irrwell):
            crop_nodes = record["node"]
            crop_vks[ix] = np.sum(vks[crop_nodes] * area[crop_nodes])
            crop_pet[ix] = np.sum(pet[crop_nodes] * area[crop_nodes])
            crop_aet[ix] = np.sum(aet[crop_nodes] * area[crop_nodes])

        crop_aet = np.where(crop_vks < 1e-30, crop_pet, crop_aet)

        if delt < 1:
            crop_aet = crop_pet


        # todo: Change the dimensionalization of self._aetold, this needs to
        #   be of numirrwells dimension! Change the dimension of self._sup
        #   and self._supold. This can be set to irrwell_num.size
        if kiter == 0:
            sup = np.zeros(self.well_max_q.shape)[self.irrwell_num]
            if kstp == 0:
                supold = np.zeros(self.well_max_q.shape)[self.irrwell_num]
                self._aetold = crop_aet * 0
            else:
                supold = self._supold[self.irrwell_num]
        else:
            sup = self._sup[self.irrwell_num]
            supold = self._supold[self.irrwell_num]

        factor = self.__calculate_factor(crop_pet, crop_aet, self._aetold, sup, supold, kiter, kstp)

        qonly = np.where(sup + factor > crop_vks, crop_vks, sup + factor)

        self._supold[self.irrwell_num] = sup
        self._aetold = crop_aet

        pumping = np.where(qonly > np.abs(self.well_max_q[self.irrwell_num]),
                           self.well_max_q[self.irrwell_num],
                           -1 * qonly)

        pumping = np.where(np.abs(pumping) < 1e-10,
                           0,
                           pumping)

        self._sup[self.irrwell_num] = np.abs(pumping)

        wells = mf6.get_value(self.well_addr)
        wells[self.irrwell_num] = pumping
        self.well_out[kstp, self.irrwell_num] = pumping
        self.well_demand_out[kstp, self.irrwell_num] = pumping
        mf6.set_value(self.well_addr, wells)
        return pumping

    def supplemental_pumping(self, mf6, conj_demand, divflow, kstp, delt=1.):
        """
        Method to calculate suplemental pumping in a conjuctive use
        scenario

        Parameters
        ----------
        mf6 : ModflowApi object
        conj_demand : np.ndarray
            crop demand for conjunctive use
        divflow : np.ndarray
            diversion flow volume
        kstp : int
            time step
        """
        # todo: check that this completely correct. not sure that the dimesnsionalization is...
        sup_demand = np.zeros((len(self.supwell_num),))
        div_list = []
        for ix, well in enumerate(self.supwell):
            segid = well['segid']
            div_list += [i for i in segid]

            sup_demand[ix] = np.sum(
                well['fracsup'] * (well['fracsupmax'] * (conj_demand[segid] - divflow[segid]))
            )

        sup_demand = np.where(sup_demand > 0, sup_demand, 0)

        if self.trigger:
            sup_demand = np.where(self.sfr_timeinperiod - delt < self.sfr_irrperiod,
                                  sup_demand,
                                  0)

        wells = mf6.get_value(self.well_addr)
        sup_pump = np.where(sup_demand + wells[self.supwell_num] > np.abs(self.well_max_q[self.supwell_num]),
                            -1 * np.abs(self.well_max_q[self.supwell_num]) - np.abs(wells[self.supwell_num]),
                            -1 * sup_demand)

        wells[self.supwell_num] += sup_pump
        wells = np.where(wells < self.well_max_q,
                         self.well_max_q,
                         wells)

        # if self.etdemand:
        #     self._idsupold[div_list] += -1 * wells[self.supwell_num][:, 0]

        self.well_demand_out[kstp, self.supwell_num] = sup_pump
        self.well_out[kstp, self.supwell_num] = sup_pump
        mf6.set_value(self.well_addr, wells)
        return sup_pump

    def __calculate_factor(self, crop_pet, crop_aet, aetold, sup, supold, kiter, kstp, accel=1):
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

        with open("mf6_debug.out", "a") as foo:
            s = f"{kstp + 1},{kiter + 1},{crop_pet[0]},{dq[0]},{det[0]},{crop_aet[0]},{aetold[0]},{sup[0]},{supold[0]},{factor[0]}\n"
            foo.write(s)

        return factor

    def conjunctive_demand(self, mf6, kstp, delt=1, kiter=1):
        """
        Method to calculate triggered surface water demand.

        Parameters
        ----------
        mf6 : ModflowApi object
        delt : float
        """
        pet = mf6.get_value(self.pet_addr)
        dsflow = mf6.get_value(self.sfr_dsflow_addr)
        aet = mf6.get_value(self.aet_addr)
        gwet = mf6.get_value(self.uz_gwet_addr)
        area = mf6.get_value(self.area_addr)
        uzet = mf6.get_value(self.uzet_addr)

        crop_pet = np.zeros((len(self.irrdiversion_num),))
        crop_aet = np.zeros((len(self.irrdiversion_num),))
        for ix, record in enumerate(self.irrdiversion):
            crop_nodes = record["node"]
            crop_pet[ix] = np.sum(pet[crop_nodes] * area[crop_nodes])
            crop_uzet = uzet[crop_nodes] / delt
            crop_gwet = gwet[crop_nodes]
            crop_aet[ix] = np.sum((crop_gwet + crop_uzet))
            # crop_aet[ix] = np.sum(aet[crop_nodes] * area[crop_nodes])

        demand = self.sfr_max_q

        if self.trigger:
            factor = np.where(crop_pet > 1e-30,
                              crop_aet / crop_pet,
                              1)

            with open("div_factor.txt", "a") as foo:
                foo.write(f"{factor[0]}, {kstp + 1}, {kiter}, {crop_aet[0]}, {crop_pet[0]}\n")

            if (kstp, kiter) == (0, 0):
                dvflw = demand * 0
            else:
                mask = np.where((self.sfr_timeinperiod > self.sfr_irrperiod) &
                                (factor <= self.sfr_triggerfact),
                                True,
                                False)

                self.sfr_timeinperiod[mask] = 0

                dvflw = np.where(self.sfr_timeinperiod - delt < self.sfr_irrperiod,
                                 demand[self.irrdiversion_num],
                                 0)
        else:
            dvflw = demand[self.irrdiversion_num]

        div_info = self.div_info[self.irrdiversion_num]
        fmaxflow = dsflow[div_info["rno"]]

        dvflw = np.where(dvflw > fmaxflow, fmaxflow, dvflw)

        divflow = mf6.get_value(self.sfr_divflow_addr)
        divflow[self.irrdiversion_num] = dvflw
        mf6.set_value(self.sfr_divflow_addr, divflow)

        return demand, divflow

    def gw_demand(self, mf6, kstp, delt, kiter=1):
        """
        Method to calculate pumping demand from groundwater using
        the trigger option

        Parameters
        ----------
        mf6 : modflowApi object
        kstp : int
        delt : float
            timestep length
        """
        pet = mf6.get_value(self.pet_addr)
        area = mf6.get_value(self.area_addr)
        aet = mf6.get_value(self.aet_addr)

        crop_pet = np.zeros((len(self.irrwell_num),))
        crop_aet = np.zeros((len(self.irrwell_num),))
        for ix, record in enumerate(self.irrwell):
            crop_nodes = record["node"]
            crop_pet[ix] = np.sum(pet[crop_nodes] * area[crop_nodes])
            crop_aet[ix] = np.sum(aet[crop_nodes] * area[crop_nodes])

        demand = self.well_max_q[self.irrwell_num]

        if self.trigger:
            factor = np.where(crop_pet > 1e-30,
                              crop_aet/crop_pet,
                              1)

            with open("factor.txt", "a") as foo:
                foo.write(f"{factor[0]}, {kstp}, {kiter}\n")

            if (kstp, kiter) == (0, 0):
                pumpage = demand * 0
            else:
                mask = np.where((self.well_timeinperiod > self.well_irrperiod) &
                                (factor < self.well_triggerfact),
                                True,
                                False)

                self.well_timeinperiod[mask] = 0

                pumpage = np.where(self.well_timeinperiod - delt < self.well_irrperiod,
                                   demand,
                                   0)

        else:
            pumpage = demand

        wells = mf6.get_value(self.well_addr)
        wells[self.irrwell_num] = pumpage
        self.well_out[kstp, self.irrwell_num] = pumpage
        self.well_demand_out[kstp, self.irrwell_num] = demand
        mf6.set_value(self.well_addr, wells)
        return pumpage

    def apply_efficiency_factors(self, mf6, finfold, pumpage=None, divflow=None, sup_pump=None):
        """
        Method to apply efficiency factors to diversion calculations

        Parameters
        ----------
        mf6 : ModflowApi
        pumpage : np.ndarray, None
            numpy array of pumping values from wells
        divflow : np.ndarray, None
            numpy array of diversion values from diversion segments
        sup_pump : np.ndarray, None
            numpy array of supplemental pumping values from wells
        finfold : np.ndarray
        """
        area = mf6.get_value(self.area_addr)
        finf = mf6.get_value(self.uz_finf_addr)
        # todo: work on putting together output for ag_package

        return_rates = np.zeros((self.gwf.modelgrid.ncpl,))
        if pumpage is not None and pumpage:
            for ix, record in enumerate(self.irrwell):
                crop_nodes = record["node"]
                eff_fact = record['eff_fact']
                field_fact = record['field_fact']
                crop_area = area[crop_nodes]
                vol = np.ones((len(crop_nodes),)) * pumpage[ix]
                subvol = ((1 - eff_fact) * vol * field_fact)
                subrate = -1 * subvol / crop_area
                return_rates[crop_nodes] += subrate

        if divflow is not None and divflow:
            for ix, record in enumerate(self.irrdiversion):
                crop_nodes = record["node"]
                eff_fact = record['eff_fact']
                field_fact = record['field_fact']
                crop_area = area[crop_nodes]
                vol = np.ones((len(crop_nodes),)) * divflow[ix]
                subvol = ((1 - eff_fact) * vol * field_fact)
                subrate = subvol / crop_area
                return_rates[crop_nodes] += subrate

            if sup_pump is not None and sup_pump:
                for ix, record in enumerate(self.supwell):
                    segid = record['segid']
                    for seg in segid:
                        idx = np.where(self.irrdiversion_num == seg)[0]
                        for idxx in idx:
                            irr_record = self.irrdiversion[idxx]
                            crop_nodes = irr_record["node"]
                            eff_fact = irr_record['eff_fact']
                            field_fact = irr_record['field_fact']
                            crop_area = area[crop_nodes]

                            vol = np.ones((len(crop_nodes),)) * sup_pump[ix]
                            subvol = ((1 - eff_fact) * vol * field_fact)
                            subrate = -1 * subvol / crop_area
                            # note: not applying return from SUP allows for trend match in via triggered irrigation....
                            #  this is most likely applied in NWT/GSFLOW during the IRRWELL LOOP
                            #  and then in effect not properly applied in the case of no IRRWELLS!
                            # return_rates[crop_nodes] += subrate

        # finf = finfold + return_rates
        finf = finf + return_rates
        mf6.set_value(self.uz_finf_addr, finf)
        return return_rates

    def run_model(self, dll,):
        mf6 = ModflowApi(dll, working_directory=self.sim.simulation_data.mfpath.get_sim_path())
        mf6.initialize()
        prev_time = 0
        current_time = mf6.get_current_time()
        end_time = mf6.get_end_time()

        input_vars = mf6.get_input_var_names()
        output_vars = mf6.get_output_var_names()

        with open("input_vars.txt", "w") as foo:
            for i in input_vars:
                foo.write(f"{i}\n")

        with open("output_vars.txt", "w") as foo:
            for i in output_vars:
                foo.write(f"{i}\n")

        self.create_addresses(mf6)
        # todo: maybe add this to the create_addresses routine
        max_iter = mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))

        # prepare the iteration loops
        kper = 0
        kstp = 0
        pumping = []
        pumping2 = []
        div = []
        sup_p = []
        while current_time < end_time:
            delt = current_time - prev_time
            dt = mf6.get_time_step()
            mf6.prepare_time_step(dt)
            kiter = 0

            if current_time in self.totim or current_time == 0.:
                if current_time == 0:
                    if self.sim_wells:
                        self.set_max_q_well(mf6)
                        self.set_irrwell_stress_period_data(kper)
                    if self.sim_diversions:
                        self.set_max_q_diversion(mf6)
                        self.set_irrdiversion_stress_period_data(kper)
                        if self.sim_supplemental:
                            if not self.sim_wells:
                                self.set_max_q_well(mf6)
                            self.set_supwell_stress_period_data(kper)

                if current_time in self.totim:
                    kper += 1
                    if self.sim_wells:
                        self.set_irrwell_stress_period_data(kper)
                    if self.sim_diversions:
                        self.set_max_q_diversion(mf6)
                        self.set_irrdiversion_stress_period_data(kper)
                        if self.sim_supplemental:
                            if not self.sim_wells:
                                self.set_max_q_well(mf6)
                            self.set_supwell_stress_period_data(kper)

            n_solutions = mf6.get_subcomponent_count()
            for sol_id in range(1, n_solutions + 1):

                mf6.prepare_solve(sol_id)
                finfold = mf6.get_value(self.uz_finf_addr).copy()
                while kiter < max_iter:
                    if self.etdemand:
                        if kiter == 0:
                            has_converged = mf6.solve(sol_id)
                        well_demand, divflow, sup_demand = None, None, None
                        if self.sim_wells:
                            well_demand = self.gw_demand_etdemand(mf6, kstp, delt, kiter)
                        if self.sim_diversions:
                            conj_demand, divflow = self.conjuctive_demand_etdemand(mf6, kstp, delt, kiter)
                            if self.sim_supplemental:
                                sup_demand = self.supplemental_pumping(mf6, conj_demand, divflow, kstp, delt)

                    else:
                        if (kstp, kiter) == (0, 0):
                            has_converged = mf6.solve(sol_id)
                        well_demand, divflow, sup_demand = None, None, None
                        if self.sim_wells:
                            well_demand = self.gw_demand(mf6, kstp, delt, kiter)
                        if self.sim_diversions:
                            conj_demand, divflow = self.conjunctive_demand(mf6, kstp, delt, kiter)
                            if self.sim_supplemental:
                                sup_demand = self.supplemental_pumping(mf6, conj_demand, divflow, kstp, delt)

                    self.apply_efficiency_factors(mf6, finfold, well_demand, divflow, sup_demand)

                    has_converged = mf6.solve(sol_id)
                    kiter += 1
                    if has_converged:
                        break

                mf6.finalize_solve(sol_id)

            self.sfr_timeinperiod = np.where(self.sfr_timeinperiod - delt < self.sfr_irrperiod,
                                             self.sfr_timeinperiod + delt,
                                             self.sfr_timeinperiod)

            self.well_timeinperiod = np.where(self.well_timeinperiod - delt < self.well_irrperiod,
                                              self.well_timeinperiod + delt,
                                              self.well_timeinperiod)

            # print(well_demand)
            # sup_p.append(sup_demand)
            # div.append(divflow)
            # pumping.append(conj_demand)
            # pumping2.append(well_demand)
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

        # print(f"conjunctive_demand 2-ac almonds {np.nansum(pumping) * 0.000810714 :.2f}")
        # print(f"sfr use 2-ac almonds {np.sum(div) * 0.000810714 :.2f}")
        # print(f"supplemental pumping use 2-ac almonds {np.sum(sup_p) * 0.000810714 :.2f}")
        # print(f"gw pumping 2-ac almonds {np.sum(pumping2) * 0.000810714 :.2f}")
        return success
