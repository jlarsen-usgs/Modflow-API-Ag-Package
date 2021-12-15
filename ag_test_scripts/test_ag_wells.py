import flopy
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sws = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(sws, "..", "develop_ag_v2"))
from develop_mf6_AG2 import Modflow6Ag



def build_mf6(name):
    sim_ws = os.path.join(sws, "..", "data", "mf6_test_ag_wells")
    sim = flopy.mf6.MFSimulation(name, sim_ws=sim_ws)

    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    period_data = [(i, i, 1.0) for i in perlen]
    tdis = flopy.mf6.ModflowTdis(
        sim,
        nper=12,
        perioddata=tuple(period_data),
        time_units="days"
    )

    ims = flopy.mf6.ModflowIms(sim, complexity="COMPLEX")

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
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True, icelltype=1)
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1)

    stress_period_data = {
        i: [[(0, 4, 4), -10.], [(0, 9, 9), -10.], [(0, 6, 6), -50.]] for i in
        range(12)
    }
    wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data=stress_period_data)

    cimis_data = os.path.join("..", "data", "davis_monthly_ppt_eto.txt")
    df = pd.read_csv(cimis_data)

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

    budget_file = f"{name}.cbc"
    head_file = f"{name}.hds"
    saverecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(10)}
    printrecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(10)}
    oc = flopy.mf6.ModflowGwfoc(gwf,
                                budget_filerecord=budget_file,
                                head_filerecord=head_file,
                                saverecord=saverecord,
                                printrecord=printrecord)

    ag = build_ag_package(name, model_ws=sim_ws)
    sim.write_simulation()
    ag.write_file()
    return sim, gwf, ag


def build_nwt(name):
    model_ws = os.path.join(sws, "..", "data", "nwt_test_ag_wells")

    ml = flopy.modflow.Modflow(
        name,
        version="mfnwt",
        exe_name="mfnwt",
        model_ws=model_ws
    )

    nrow = 10
    ncol = 10

    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    # period_data = [(i, i, 1.0) for i in perlen]

    dis = flopy.modflow.ModflowDis(
        ml,
        nlay=1,
        nrow=10,
        ncol=10,
        nper=12,
        delr=63.6,
        delc=63.6,
        top=100,
        perlen=list(perlen),
        nstp=list(perlen),
        steady=False,
        itmuni=4,
        lenuni=2
    )

    bas = flopy.modflow.ModflowBas(
        ml,
        ibound=np.ones((nrow, ncol), dtype=int),
        strt=np.ones((nrow, ncol)) * 95,
    )

    upw = flopy.modflow.ModflowUpw(
        ml,
        ipakcb=52,
        laytyp=1,
        hk=1,
        ss=1e-05,
        sy=0.15,
    )

    stress_period_data = {
        i: [[(0, 4, 4), -10.], [(0, 9, 9), -10.], [(0, 6, 6), -50.]] for i in
        range(12)
    }

    # build a UZF package
    cimis_data = os.path.join("..", "data", "davis_monthly_ppt_eto.txt")
    df = pd.read_csv(cimis_data)

    finf = {}
    pet = {}
    extdp = {}
    extwc = {}
    for i in range(12):
        finf[i] = np.ones((nrow, ncol)) * df.ppt_avg_m.values[i] / perlen[i]
        pet[i] = np.ones((nrow, ncol)) * df.eto_avg_m.values[i] / perlen[i]
        extdp[i] = np.ones((nrow, ncol)) * 4
        extwc[i] = np.ones((nrow, ncol)) * 0.06

    uzf = flopy.modflow.ModflowUzf1(
        ml,
        nuztop=1,
        iuzfopt=1,
        irunflg=0,
        ipakcb=52,
        ietflg=1,
        ntrail2=7,
        nsets=1000,
        surfdep=0.33,
        iuzfbnd=np.ones((nrow, ncol)),
        vks=8.64,
        eps=5,
        thts=0.35,
        thtr=0.05,
        thti=0.08,
        specifythtr=True,
        specifythti=True,
        finf=finf,
        pet=pet,
        extdp=extdp,
        extwc=extwc

    )

    stress_period_data = {}
    for kper, nts in enumerate(perlen):
        for ts in range(nts):
            stress_period_data[(kper, ts)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(
        ml,
        stress_period_data=stress_period_data
    )

    nwt = flopy.modflow.ModflowNwt(
        ml
    )

    ag = build_ag_package(name, ml, model_ws)
    ml.write_input()
    nam = os.path.join(model_ws, f"{name}.nam")
    with open(nam, "a") as foo:
        foo.write(f"DATA              53  {name}.well1.txt")

    return ml, ag


def build_ag_package(name, ml=None, model_ws=None):

    if ml is None:
        # this is mf6
        ml = flopy.modflow.Modflow(name, model_ws=model_ws, version="mfnwt")
        dis = flopy.modflow.ModflowDis(
            ml,
            nlay=1,
            nrow=10,
            ncol=10,
            nper=12,
            delr=63.6,
            delc=63.6
        )
    else:
        dis = ml.dis

    if name == "specified":
        options = flopy.utils.OptionBlock(
            "IRRIGATION_WELL 1 2 "
            "MAXWELLS 1 WELLCBC 52 TIMESERIES_WELL".lower(),
            flopy.modflow.ModflowAg
        )
    else:
        options = flopy.utils.OptionBlock(
            f"{name} IRRIGATION_WELL 1 2 "
            "MAXWELLS 1 WELLCBC 52 TIMESERIES_WELL".lower(),
            flopy.modflow.ModflowAg
        )

    timeseries_well = flopy.modflow.ModflowAg.get_empty(1, block='time series')
    timeseries_well[0] = ("well", 1, 53)

    well_list = flopy.modflow.ModflowAg.get_empty(1, block="well")
    x = [[0, 4, 4, -10.],]  # [0, 9, 9, -100.], [0, 6, 6, -50]]
    for ix, rec in enumerate(well_list):
        well_list[ix] = tuple(x[ix])

    irrwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 2, "irrwell")
        spd[0] = (0, 2, 31, 0.2, 4, 4, 0.9, 0.5, 5, 4, 0.9, 0.5)
        irrwell[i] = spd

    ag = flopy.modflow.ModflowAg(
        ml,
        options=options,
        time_series=timeseries_well,
        well_list=well_list,
        irrwell=irrwell,
        nper=dis.nper
    )

    return ag


def run_mfnwt(name):
    ml, ag = build_nwt(name)
    success, buff = ml.run_model()
    if not success:
        raise AssertionError


def run_mf6(name):
    from modflowapi import ModflowApi
    dll = os.path.join("..", "modflow-bmi", "libmf6.dll")
    sim, gwf, ag = build_mf6(name)
    # mf6 = ModflowApi(dll, working_directory=os.path.join(sws, "..", "data", "mf6_test_ag_wells"))
    # mf6.initialize()
    mf6ag = Modflow6Ag(sim, ag)
    mf6ag.run_model(dll)


def compare_output(name):
    cbc_name = "{}.cbc".format(name)
    nwt_ws = os.path.join(sws, "..", "data", "nwt_test_ag_wells")
    mf6_ws = os.path.join(sws, "..", "data", "mf6_test_ag_wells")

    nwt_cbc = os.path.join(nwt_ws, cbc_name)
    mf6_cbc = os.path.join(mf6_ws, cbc_name)

    nwt_cbc = flopy.utils.CellBudgetFile(nwt_cbc)
    print(nwt_cbc.get_unique_record_names())
    nwt_ag_well = nwt_cbc.get_data(text="AG WE")
    # print(nwt_ag_well)
    l = []
    for rec in nwt_ag_well:
        l.append(rec["q"][0])

    nwt_pump = np.sum(l)
    print(nwt_pump)

    mf6_cbc = flopy.utils.CellBudgetFile(mf6_cbc)
    print(mf6_cbc.get_unique_record_names())
    mf6_pump = mf6_cbc.get_data(text="WEL")
    l2 = []
    for rec in mf6_pump:
        l2.append(rec['q'][0])

    mf6_pump = np.sum(l2)
    print(mf6_pump)

    plt.plot(range(1, len(l) + 1), l2, "r-")
    plt.plot(range(1, len(l) + 1), l, "b--")
    plt.show()


if __name__ == "__main__":
    model_names = ("etdemand", "trigger", "specified")
    # todo: look at how PET and AET are being calculated,
    #   trigger and etdemand are showing different results from
    #   mf6AG

    run_mfnwt(model_names[1])
    run_mf6(model_names[1])
    compare_output(model_names[1])