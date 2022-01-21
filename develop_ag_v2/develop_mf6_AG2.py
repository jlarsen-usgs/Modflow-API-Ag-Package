from modflowapi import ModflowApi
import numpy as np
import pandas as pd
import flopy
import os


def create_test_model(name):

    sim_ws = os.path.join(".",)
    sim = flopy.mf6.MFSimulation(name, sim_ws=sim_ws)

    # create TDIS with monthly stress periods and daily time steps
    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    period_data = [(i, i, 1.0) for i in perlen]
    tdis = flopy.mf6.ModflowTdis(
        sim,
        nper=12,
        perioddata=tuple(period_data),
        time_units="days"
    )

    # create IMS
    ims = flopy.mf6.ModflowIms(sim, complexity="COMPLEX")

    # create model!
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        print_input=True,
        print_flows=True
    )

    # define delc and delr to equal approximately 1 acre
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nrow=10,
        ncol=10,
        delr=63.6,
        delc=63.6,
        top=100,
        length_units='meters'
    )

    ic = flopy.mf6.ModflowGwfic(gwf, strt=95)
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True)
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1)

    stress_period_data = {
        i: [[(0, 4, 4), -100.], [(0, 9, 9), -100.], [(0, 6, 6), -50.]] for i in range(12)
    }
    wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data=stress_period_data)

    # create RCH and EVT packages from DAVIS monthly average CIMIS data
    cimis_data = os.path.join("..", "data", "davis_monthly_ppt_eto.txt")
    df = pd.read_csv(cimis_data)

    # recharge = {i: v / perlen[i] for i, v in enumerate(df.ppt_avg_m.values)}
    # rch = flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)

    # surface = {i: 100 for i in range(12)}
    # eto = {i: v / perlen[i] for i, v in enumerate(df.eto_avg_m.values)}
    # depth = {i: 3 for i in range(12)}
    # evt = flopy.mf6.ModflowGwfevta(gwf, surface=surface, rate=eto, depth=depth)

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

    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        simulate_et=True,
        nuzfcells=nuzfcells,
        ntrailwaves=ntrailwaves,
        nwavesets=nwavesets,
        packagedata=package_data,
        perioddata=period_data
    )

    # create SFR package
    nreaches = 11
    package_data = []
    for i in range(nreaches):
        ustrf = 1.0
        if i in (0, 9, 10):
            ncon = 1
        elif i == 6:
            ncon = 3
        else:
            ncon = 2

        if i == 6:
            ndiv = 1
        else:
            ndiv = 0

        cellid = (0, i, 7)
        kh = 0.000015
        if i == 10:
            kh = 0.0000015
            cellid = (0, 6, 6)
            ustrf = 0.0

        rch_data = (i, cellid, 100, 5, 0.02, 99, 0.5, kh, 0.03, ncon, ustrf, ndiv)
        package_data.append(rch_data)

    connection_data = []
    for i in range(nreaches):
        if i == 0:
            cd = [i, -1 * (i + 1)]
        elif i == 9:
            cd = [i, i - 1]
        elif i == 10:
            cd = [i, 6]
        elif i == 6:
            cd = [i, i - 1, -1 * (i + 1), -10]
        else:
            cd = [i, i - 1, -1 * (i + 1)]
        connection_data.append(cd)

    diversions = [[6, 0, 10, "UPTO"]]

    period_data = {
        i: [(0, "INFLOW", 8),  (6, "DIVERSION", 0, 10)] for i in range(12)
    }

    sfr = flopy.mf6.ModflowGwfsfr(gwf,
                               nreaches=nreaches,
                               packagedata=package_data,
                               connectiondata=connection_data,
                               diversions=diversions,
                               perioddata=period_data)

    budget_file = f"{name}.cbc"
    head_file = f"{name}.hds"
    saverecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(10)}
    printrecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(10)}
    oc = flopy.mf6.ModflowGwfoc(gwf,
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
        "ETDEMAND IRRIGATION_WELL 1 2 SUPPLEMENTAL_WELL 1 1 "
        "MAXWELLS 3 IRRIGATION_DIVERSION 1 2".lower(),
        flopy.modflow.ModflowAg
    )

    well_list = flopy.modflow.ModflowAg.get_empty(2, block="well")

    x = [[0, 4, 4, -100.], [0, 9, 9, -100.], [0, 6, 6, -50]]
    for ix, rec in enumerate(well_list):
        well_list[ix] = tuple(x[ix])

    irrwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 2, "irrwell")
        spd[0] = (0, 2, 10, 0.0, 4, 4, 0, 0.5, 5, 4, 0, 0.5)
        irrwell[i] = spd

    irrdiversion = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 2, "irrdiversion")
        # todo: update triggerfact and period parameters
        spd[0] = (1, 2, 10, 0.0, 6, 6, 0, 0.5, 7, 6, 0, 0.5)
        irrdiversion[i] = spd

    supwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 1, "supwell")
        spd[0] = (2, 1, 1, 0.9, 1)
        supwell[i] = spd

    ag = flopy.modflow.ModflowAg(
        ml,
        options=options,
        well_list=well_list,
        irrwell=irrwell,
        irrdiversion=irrdiversion,
        supwell=supwell,
        nper=12
    )
    ag.write_file()
    return ag


def create_ag_package_trigger():
    """
    Method to create an TRIGGER AG package
    """
    ml = flopy.modflow.Modflow("trigger", version='mfnwt')
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
        "TRIGGER IRRIGATION_WELL 1 2 SUPPLEMENTAL_WELL 1 1 "
        "MAXWELLS 3 IRRIGATION_DIVERSION 1 2".lower(),
        flopy.modflow.ModflowAg
    )

    well_list = flopy.modflow.ModflowAg.get_empty(2, block="well")

    x = [[0, 4, 4, -100.], [0, 9, 9, -100.], [0, 6, 6, -50]]
    for ix, rec in enumerate(well_list):
        well_list[ix] = tuple(x[ix])

    irrwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 2, "irrwell")
        spd[0] = (0, 2, 10, 0.1, 4, 4, 0, 0.5, 5, 4, 0, 0.5)
        irrwell[i] = spd

    irrdiversion = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 2, "irrdiversion")
        spd[0] = (1, 2, 10, 0.1, 6, 6, 0, 0.5, 7, 6, 0, 0.5)
        irrdiversion[i] = spd

    supwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 1, "supwell")
        spd[0] = (2, 1, 1, 1, 1)
        supwell[i] = spd

    ag = flopy.modflow.ModflowAg(
        ml,
        options=options,
        well_list=well_list,
        irrwell=irrwell,
        irrdiversion=irrdiversion,
        supwell=supwell,
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
        self.nsteps = np.sum(self.gwf.modeltime.nstp)

        self.etdemand = self.ag.etdemand
        self.trigger = self.ag.trigger

        self.sim_wells = self.ag.irrigation_well
        self.nummaxwell = self.ag.nummaxwell

        self.sim_diversions = self.ag.irrigation_diversion
        self.numirrdiversion = self.ag.numirrdiversions

        self.sim_supplemental = self.ag.supplemental_well
        self.numsupwell = self.ag.numsupwells
        self.irrigation_return = []

        if self.sim_wells:
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

    def conjuctive_demand_etdemand(self, mf6, delt=1, kiter=1):
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

        if delt < 1:
            crop_aet = crop_pet

        if kiter == 0:
            sup = divflow[self.irrdiversion_num]
            supold = np.zeros(self.sfr_max_q.shape)[self.irrdiversion_num]
            self._supact = np.zeros(self.sfr_max_q.shape)
            self._idsupold = np.zeros(self.sfr_max_q.shape)
            self._aetoldsw = crop_aet
        else:
            sup = divflow[self.irrdiversion_num]
            supold = self._idsupold[self.irrdiversion_num]

        factor = self.__calculate_factor(
            crop_pet,
            crop_aet,
            self._aetoldsw,
            sup,
            supold,
            kiter,
            1
        )

        factor = np.where(factor < 0, 0, factor)

        self._aetoldsw = crop_aet
        self._idsupold[self.irrdiversion_num] = divflow[self.irrdiversion_num]
        self._supact[self.irrdiversion_num] += factor

        dvflw = np.where(self._supact >= self.sfr_max_q,
                         self.sfr_max_q,
                         self._supact)

        div_info = self.div_info[self.irrdiversion_num]
        fmaxflow = dsflow[div_info["rno"]]

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

        if kstp == 243:
            print('break')

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

        with open("factor.txt", "a") as foo:
            foo.write(f"{qonly[0]}, {kstp}, {kiter}\n")

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
        for ix, well in enumerate(self.supwell):
            segid = well['segid']
            sup_demand[ix] = np.sum(
                well['fracsup'] * (well['fracsupmax'] * (conj_demand[segid] - divflow[segid]))
            )

        sup_demand = np.where(sup_demand > 0, sup_demand, 0)

        if self.trigger:
            sup_demand = np.where(self.sfr_timeinperiod - delt < self.sfr_irrperiod,
                                  sup_demand,
                                  0)

        sup_pump = np.where(sup_demand > np.abs(self.well_max_q[self.supwell_num]),
                            np.abs(self.well_max_q[self.supwell_num]),
                            -1 * sup_demand)

        # todo: need a smart way to apply sup!
        wells = mf6.get_value(self.well_addr)
        wells[self.supwell_num] = sup_pump
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

    def conjunctive_demand(self, mf6, delt=1):
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

        crop_pet = np.zeros((len(self.irrdiversion_num),))
        crop_aet = np.zeros((len(self.irrdiversion_num),))
        for ix, record in enumerate(self.irrdiversion):
            crop_nodes = record["node"]
            crop_pet[ix] = np.sum(pet[crop_nodes])
            crop_aet[ix] = np.sum(aet[crop_nodes])

        demand = self.sfr_max_q

        if self.trigger:
            factor = np.where(crop_pet > 1e-30,
                              crop_aet / crop_pet,
                              1)

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
        # finf = mf6.get_value(self.uz_finf_addr)

        return_rates = np.zeros((self.gwf.modelgrid.ncpl,))
        if pumpage is not None and pumpage:
            for ix, record in enumerate(self.irrwell):
                crop_nodes = record["node"]
                eff_fact = record['eff_fact']
                field_fact = record['field_fact']
                crop_area = area[crop_nodes]
                vol = np.ones((len(crop_nodes),)) * pumpage[ix] / len(crop_nodes)
                subvol = ((1 - eff_fact) * vol * field_fact)
                subrate = -1 * subvol / crop_area
                return_rates[crop_nodes] += subrate

        if divflow is not None and divflow:
            for ix, record in enumerate(self.irrdiversion):
                crop_nodes = record["node"]
                eff_fact = record['eff_fact']
                field_fact = record['field_fact']
                crop_area = area[crop_nodes]
                vol = np.ones((len(crop_nodes),)) * divflow[ix] / len(crop_nodes)
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

                            vol = np.ones((len(crop_nodes),)) * sup_pump[ix] / len(crop_nodes)
                            subvol = ((1 - eff_fact) * vol * field_fact)
                            subrate = -1 * subvol / crop_area
                            return_rates[crop_nodes] += subrate

        finf = finfold + return_rates
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
                            conj_demand, divflow = self.conjuctive_demand_etdemand(mf6, delt, kiter)
                            if self.sim_supplemental:
                                sup_demand = self.supplemental_pumping(mf6, conj_demand, divflow, kstp, delt)

                    else:
                        well_demand, divflow, sup_demand = None, None, None
                        if self.sim_wells:
                            well_demand = self.gw_demand(mf6, kstp, delt, kiter)
                        if self.sim_diversions:
                            conj_demand, divflow = self.conjunctive_demand(mf6, delt)
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


if __name__ == "__main__":
    dll = os.path.join("..", "modflow-bmi", "libmf6.dll")
    sim, gwf = create_test_model("GWF")
    ag = create_ag_package_etdemand()
    ag2 = create_ag_package_trigger()

    mf6ag = Modflow6Ag(sim, ag2)
    mf6ag.run_model(dll, )

    # todo: create methods to store and write
    #  demands, pumpage, sfr, and diversions

    # todo: setup tests for MF6 AG package (compare to mf-nwt AG1).
    # todo: create methods that allow this class to be used by other
    #   driver classes (ie. let it be called by another class, but
    #   still be able to perform the setup etc....
