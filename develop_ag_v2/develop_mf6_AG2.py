from modflowapi import ModflowApi
import numpy as np
import pandas as pd
import flopy
import os


def create_test_model(name):
    import flopy as fp

    sim_ws = os.path.join(".",)
    sim = fp.mf6.MFSimulation(name, sim_ws=sim_ws)

    # create TDIS with monthly stress periods and daily time steps
    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    period_data = [(i, i, 1.0) for i in perlen]
    tdis = fp.mf6.ModflowTdis(
        sim,
        nper=12,
        perioddata=tuple(period_data),
        time_units="days"
    )

    # create IMS
    ims = fp.mf6.ModflowIms(sim, complexity="MODERATE")

    # create model!
    gwf = fp.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        print_input=True,
        print_flows=True
    )

    # define delc and delr to equal approximately 1 acre
    dis = fp.mf6.ModflowGwfdis(
        gwf,
        nrow=10,
        ncol=10,
        delr=63.6,
        delc=63.6,
        top=100,
        length_units='meters'
    )

    ic = fp.mf6.ModflowGwfic(gwf, strt=95)
    npf = fp.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True)
    sto = fp.mf6.ModflowGwfsto(gwf, iconvert=1)

    stress_period_data = {
        i: [[(0, 4, 4), -100.], [(0, 9, 9), -100.]] for i in range(12)
    }
    wel = fp.mf6.ModflowGwfwel(gwf, stress_period_data=stress_period_data)

    # create RCH and EVT packages from DAVIS monthly average CIMIS data
    cimis_data = os.path.join("..", "data", "davis_monthly_ppt_eto.txt")
    df = pd.read_csv(cimis_data)

    recharge = {i: v / perlen[i] for i, v in enumerate(df.ppt_avg_m.values)}
    rch = fp.mf6.ModflowGwfrcha(gwf, recharge=recharge)

    surface = {i: 100 for i in range(12)}
    eto = {i: v / perlen[i] for i, v in enumerate(df.eto_avg_m.values)}
    depth = {i: 3 for i in range(12)}
    evt = fp.mf6.ModflowGwfevta(gwf, surface=surface, rate=eto, depth=depth)

    # build a UZF package
    nuzfcells = 100
    ntrailwaves = 7
    nwavesets = 40
    package_data = []
    cnt = 0
    for i in range(10):
        for j in range(10):
            rec = (cnt, (0, i, j), 1, 0, 0.33, 8.64, 0.05, 0.35, 0.08, 5)
            package_data.append(rec)
            cnt += 1

    period_data = {}
    for i in range(12):
        cnt = 0
        spd = []
        for _ in range(10):
            for _ in range(10):
                rec = (
                    cnt,
                    df.ppt_avg_m.values[i]/perlen[i],
                    df.eto_avg_m.values[i]/perlen[i],
                    4,
                    0.06,
                    100,
                    100,
                    3
                )
                spd.append(rec)
                cnt += 1
        period_data[i] = spd

    uzf = fp.mf6.ModflowGwfuzf(
        gwf,
        simulate_et=True,
        nuzfcells=nuzfcells,
        ntrailwaves=ntrailwaves,
        nwavesets=nwavesets,
        packagedata=package_data,
        perioddata=period_data
    )

    budget_file = f"{name}.cbc"
    head_file = f"{name}.hds"
    saverecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(10)}
    printrecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(10)}
    oc = fp.mf6.ModflowGwfoc(gwf,
                             budget_filerecord=budget_file,
                             head_filerecord=head_file,
                             saverecord=saverecord,
                             printrecord=printrecord)

    sim.write_simulation()

    return sim, gwf


def create_ag_package_etdemand():
    """
    Method to create an ETDEMAND AG package
    """
    ml = flopy.modflow.Modflow("etdemand", version='mfnwt')
    dis = flopy.modflow.ModflowDis(
        ml,
        nlay=1,
        nrow=10,
        ncol=10,
        nper=12,
        delr=63.6,
        delc=63.6
    )

    options = flopy.utils.OptionBlock(
        "ETDEMAND IRRIGATION_WELL 1 2 MAXWELLS 2".lower(),
        flopy.modflow.ModflowAg
    )

    well_list = flopy.modflow.ModflowAg.get_empty(2, block="well")

    x = [[0, 4, 4, -100.], [0, 9, 9, -100.]]
    for ix, rec in enumerate(well_list):
        well_list[ix] = tuple(x[ix])

    irrwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 2, "irrwell")
        spd[0] = (0, 2, 0, 1, 0, 0, 0, 0.5, 0, 1, 0, 0.5)
        irrwell[i] = spd

    ag = flopy.modflow.ModflowAg(
        ml,
        options=options,
        well_list=well_list,
        irrwell=irrwell,
        nper=12
    )
    ag.write_file()
    return ag


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

        self.etdemand = self.ag.etdemand
        self.trigger = self.ag.trigger
        self.sim_wells = self.ag.maxwells
        self.nummaxwell = self.ag.nummaxwell
        self.irrwell = self.ag.irrwell
        self.nsteps = np.sum(self.gwf.modeltime.nstp)

        self.sim_diversions = False

        if self.sim_wells:
            self.well_max_q = np.zeros((self.nummaxwell,))
            self.well_q = np.zeros((self.nsteps, self.nummaxwell))
            if self.trigger:
                self.well_timeinperiod = np.zeros(
                    (self.nsteps, self.nummaxwell)
                )
                self.well_irrperiod = np.zeros(
                    (self.nsteps, self.nummaxwell)
                )
            self._sup = np.zeros((self.nummaxwell))
            self._supold = np.zeros((self.nummaxwell))
            self.irrwell_num = None
            self.irrwell = None

        if self.sim_diversions:
            self.sfr_max_q = None
            self.sfr_timeinperiod = None
            self.sfr_irrperiod = None
            self.sfr_q = None

    def create_addresses(self, mf6):
        sto_name = self.gwf.sto.name[0].upper()
        dis_name = self.gwf.dis.name[0].upper()
        try:
            rch_name = self.gwf.rcha.name[0].upper()
        except AttributeError:
            rch_name = self.gwf.rch.name[0].upper()

        try:
            evt_name = self.gwf.evta.name[0].upper()
        except AttributeError:
            evt_name = self.gwf.evt.name[0].upper()

        try:
            uzf_name = self.gwf.uzf.name[0].upper()
        except AttributeError:
            uzf_name = None

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
        self.rch_addr = mf6.get_var_address("BOUND", self.name, rch_name)
        self.evt_addr = mf6.get_var_address("BOUND", self.name, evt_name)
        self.head_addr = mf6.get_var_address("X", self.name)
        self.area_addr = mf6.get_var_address("UZFAREA", self.name, uzf_name)
        self.vks_addr = mf6.get_var_address("VKS", self.name, uzf_name)
        self.pet_addr = mf6.get_var_address("PET", self.name, uzf_name)
        self.uzet_addr = mf6.get_var_address("UZET", self.name, uzf_name)
        self.uz_gwet_addr = mf6.get_var_address("GWET", self.name, uzf_name)

    def set_max_q_well(self, mf6):
        """
        Method to set the max-q of AG wells

        Parameters
        ----------
        mf6 : ModflowApi object
        """
        well = mf6.get_value(self.well_addr)
        self.well_max_q = well.T[0]
        # change max_q to zero and reset based on AG-Demand
        well[0, :] = 0
        mf6.set_value(self.well_addr, well)

    def set_irrwell_stress_period_data(self, kper):
        """
        Method to set stress period data from wells

        Parameters:
        ----------
        kper : int
            stress period number (zero based)
        """
        self.irrwell_num = self.ag.irrwell[kper]['wellid']
        irrwell_info = []
        for rec in self.ag.irrwell[kper]:
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
                                                   ("eff_fact", int),
                                                   ("field_fact", int)]))
        self.irrwell = irrwell_info

    def conjuctive_demand_uzf(self, mf6, crop_nodes, aetold, supold, delt=1, kiter=1):
        """
        Method to determine conjunctive use demand

        Parameters
        ----------
        mf6 : ModflowApi object
        crop_nodes : np.ndarray
        aetold : float
        supold : float
        delt : float
        timestep : float
        """
        pet = mf6.get_value(self.pet_addr)
        vks = mf6.get_value(self.vks_addr)
        area = mf6.get_value(self.area_addr)
        gwet = mf6.get_value(self.uz_gwet_addr)
        uzet = mf6.get_value(self.uzet_addr)

        sup = 1
        # need to get number of ag-nodes, then calculate area
        crop_vks = vks[crop_nodes] * area[crop_nodes]
        crop_pet = np.sum(pet[crop_nodes] * area[crop_nodes])
        crop_uzet = uzet[crop_nodes]
        crop_gwet = gwet[crop_nodes]
        crop_aet = np.sum(crop_gwet + crop_uzet / delt)

        factor = self.__calculate_factor(
            crop_pet,
            crop_aet,
            aetold,
            sup,
            supold,
            kiter
        )

        return factor, crop_aet, sup

    def gw_demand_uzf(self, mf6, delt=1, kiter=1):
        """
        Method to determine groundwater use demand

        Parameters
        ----------
        mf6 : ModflowApi object
        crop_nodes : np.ndarray
        aetold : float
        supold : float
        delt : float
        timestep : float
        """
        pet = mf6.get_value(self.pet_addr)
        vks = mf6.get_value(self.vks_addr)
        area = mf6.get_value(self.area_addr)
        gwet = mf6.get_value(self.uz_gwet_addr)
        uzet = mf6.get_value(self.uzet_addr)

        crop_vks = np.zeros((len(self.irrwell_num,)))
        crop_pet = np.zeros((len(self.irrwell_num,)))
        crop_aet = np.zeros((len(self.irrwell_num)))
        for ix, record in enumerate(self.irrwell):
            crop_nodes = record["node"]
            crop_vks[ix] = np.sum(vks[crop_nodes] * area[crop_nodes])
            crop_pet[ix] = np.sum(pet[crop_nodes] * area[crop_nodes])
            crop_uzet = uzet[crop_nodes] / delt
            crop_gwet = gwet[crop_nodes]
            crop_aet[ix] = np.nansum(crop_gwet + crop_uzet)

        crop_aet = np.where(crop_vks < 1e-30, crop_pet, crop_aet)

        # todo: maybe change this to current time
        if delt < 1:
            crop_aet = crop_pet
            if kiter == 0:
                sup = self.well_max_q[self.irrwell_num]
                supold = self.well_max_q[self.irrwell_num]
                self._aetold = crop_aet
            else:
                sup = self._sup[self.irrwell_num]
                supold = self._supold[self.irrwell_num]
        else:
            sup = self._sup[self.irrwell_num]
            supold = self._supold[self.irrwell_num]

        factor = self.__calculate_factor(crop_pet, crop_aet, self._aetold, sup, supold, kiter)

        # todo: need to add the factor = supold + factor term to finish eq. 17
        self._supold[self.irrwell_num] = sup
        self._sup[self.irrwell_num] = sup + factor
        self._aetold = crop_aet

        # todo: need to reset factor using a numpy where statement.
        factor = np.where(factor > crop_vks, crop_vks, factor)

        return factor

    def __calculate_factor(self, crop_pet, crop_aet, aetold, sup, supold, kiter, accel=1):
        """
        Method to calculate unmet demand for crops

        Parameters
        ----------
        crop_pet :
        crop_aet :
        aetold :
        sup :
        supold :
        timestep : float

        """
        et_diff = crop_pet - crop_aet
        factor = et_diff
        det = crop_aet - aetold
        dq = sup - supold

        if kiter > 0:
            factor = np.where(np.abs(det) > 1e-05,
                              dq * et_diff / det,
                              factor)

        factor = np.where(factor < et_diff, et_diff, factor)
        factor = np.where(factor < 1e-05, 1e-05, factor)
        factor = np.where(factor > et_diff * accel,
                          et_diff * accel,
                          factor)

        return factor

    def sw_demand_trigger(self, mf6, crop_nodes, delt):
        """
        Method to calculate triggered surface water demand.

        Parameters
        ----------
        mf6 : ModflowApi object
        crop_nodes : np.ndarray
        delt : float
        timestep : float
        """
        # todo: will need to get iseg from SFR or wrap and pass in iseg at
        #   some point!
        if delt < 1:
            return 0, 0, 0

        pet = mf6.get_value(self.pet_addr)
        area = mf6.get_value(self.area_addr)
        gwet = mf6.get_value(self.uz_gwet_addr)
        uzet = mf6.get_value(self.uzet_addr)

        crop_pet = np.sum(pet[crop_nodes] * area[crop_nodes])
        crop_uzet = uzet[crop_nodes]
        crop_gwet = gwet[crop_nodes]
        crop_aet = np.sum((crop_gwet + crop_uzet / delt) / area)

        # todo: need to set pet/aet to an array based on SFR iseg!
        #   then we add each crop's demand to the SFR iseg array before
        #   performing calculations.
        return 0, 0, 0

    def gw_demand_trigger(self,  mf6, crop_nodes, delt):
        if delt < 1:
            return 0, 0, 0

        pet = mf6.get_value(self.pet_addr)
        area = mf6.get_value(self.area_addr)
        gwet = mf6.get_value(self.uz_gwet_addr)
        uzet = mf6.get_value(self.uzet_addr)

        crop_pet = np.sum(pet[crop_nodes] * area[crop_nodes])
        crop_uzet = uzet[crop_nodes]
        crop_gwet = gwet[crop_nodes]
        crop_aet = np.sum((crop_gwet + crop_uzet / delt) / area)

        factor = 1
        if crop_pet > 1e-30:
            factor = crop_aet/crop_pet

        # todo: need to figure out clafication for the remainder of this
        #   equation.
        return 0, 0, 0

    def run_model(self, dll,):
        # todo: or alternatively, we can the well_list from a modflow WEL package
        mf6 = ModflowApi(dll)
        mf6.initialize()
        prev_time = 0
        current_time = mf6.get_current_time()
        end_time = mf6.get_end_time()

        """
        input_vars = mf6.get_input_var_names()
        output_vars = mf6.get_output_var_names()

        with open("input_vars.txt", "w") as foo:
            for i in input_vars:
                foo.write(f"{i}\n")

        with open("output_vars.txt", "w") as foo:
            for i in output_vars:
                foo.write(f"{i}\n")
        """
        self.create_addresses(mf6)
        # todo: add this to the create_addresses routine
        max_iter = mf6.get_value(mf6.get_var_address("MXITER", "SLN_1"))

        # prepare the iteration loops
        kper = 0
        kstp = 0
        aetold = 0
        supold = 0
        pumping = []
        pumping2 = []
        while current_time < end_time:
            delt = current_time - prev_time
            dt = mf6.get_time_step()
            mf6.prepare_time_step(dt)
            kiter = 0

            if current_time in self.totim or current_time == 0.:
                if current_time == 0 and self.sim_wells:
                    self.set_max_q_well(mf6)
                    self.set_irrwell_stress_period_data(kper)

                if current_time in self.totim:
                    kper += 1
                    self.set_irrwell_stress_period_data(kper)
            # conjunctive_demand, aetold, supold = self.conjuctive_demand_uzf(
            #     mf6, [0, 1], aetold, supold, delt=delt, timestep=current_time
            # )

            n_solutions = mf6.get_subcomponent_count()
            for sol_id in range(1, n_solutions + 1):

                mf6.prepare_solve(sol_id)
                while kiter < max_iter:

                    if self.etdemand:
                        if self.sim_wells:
                            well_demand = self.gw_demand_uzf(mf6, delt, kiter)
                        if self.sim_diversions:
                            conj_demand = self.conjuctive_demand_uzf(mf6, delt, kiter)
                    # todo: need to solve ETDEMAND in here!!!!!!!!!
                    #   and set up the newtonian iteration in python!!!!!
                    has_converged = mf6.solve(sol_id)
                    kiter += 1
                    if has_converged:
                        break

                print(well_demand)
                pumping2.append(well_demand)
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

        print(f"conjunctive_demand 2-ac almonds {np.sum(pumping) * 0.000810714 :.2f}")
        print(f"gw_demand 2-ac almonds {np.sum(pumping2) * 0.000810714 :.2f}")
        return success


if __name__ == "__main__":
    dll = os.path.join("..", "modflow-bmi", "libmf6.dll")
    sim, gwf = create_test_model("GWF")
    ag = create_ag_package_etdemand()

    mf6ag = Modflow6Ag(sim, ag)
    mf6ag.run_model(dll, )

