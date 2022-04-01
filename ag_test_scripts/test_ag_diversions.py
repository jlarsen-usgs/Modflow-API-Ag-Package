import flopy
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sws = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(sws, "..", "develop_ag_v2"))
from develop_mf6_AG2 import Modflow6Ag

from math import log10, floor

def read_sfr_list_data(f):
    d = {"kper": [], "kstp": [], "iseg": [], "inflow": [], "outflow": []}
    read_header = False
    read_one_more = False
    read_data = False
    kper, kstp = None, None
    with open(f) as foo:
        for line in foo:
            if read_data:
                if "--------------" in line:
                    read_data = False
                else:
                    t = line.split()
                    d['kper'].append(kper)
                    d['kstp'].append(kstp)
                    d["iseg"].append(int(t[0]))
                    d["inflow"].append(float(t[2]))
                    d["outflow"].append(float(t[3]))
            elif read_one_more:
                read_one_more = False
                read_data = True
            elif read_header:
                if "NUMBER" in line:
                    read_one_more = True
                    read_header = False
            if "SFR_0 PACKAGE - SUMMARY OF FLOWS FOR EACH CONTROL VOLUME" in line:
                read_header = True
                t = line.split()
                kper = int(t[-3])
                kstp = int(t[-1])

    df = pd.DataFrame.from_dict(d)
    return df


def round_to_n(x, n):
    if x == 0:
        return 0
    t = round(x, -int(floor(log10(abs(x))) - (n - 1)))
    return t


def build_mf6(name):
    sim_ws = os.path.join(sws, "..", "data", "mf6_test_ag_diversions")
    sim = flopy.mf6.MFSimulation(name, sim_ws=sim_ws)

    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    period_data = [(i, i, 1.0) for i in perlen]
    tdis = flopy.mf6.ModflowTdis(
        sim,
        nper=12,
        perioddata=tuple(period_data),
        time_units="days"
    )
    #  headtol=2e-03,
    #  fluxtol=50,
    #  thickfact=1e-08,
    n = 2e-03
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        complexity="SIMPLE",
        no_ptcrecord=["ALL"],
        outer_dvclose=2e-03,
        outer_maximum=50,
        # inner_dvclose=1e-06,
        # inner_maximum=50,
        rcloserecord=[1e-10, "L2NORM_RCLOSE"],
        scaling_method="L2NORM", # "L2NORM",
        linear_acceleration="BICGSTAB",
        under_relaxation="DBD",
        under_relaxation_gamma=0.0,
        under_relaxation_theta=0.97,
        under_relaxation_kappa=0.0001
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        print_input=True,
        print_flows=True,
        newtonoptions="NEWTON UNDER_RELAXATION",
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

    # stress_period_data = {i: [[(0, 1, 1), -50.], ] for i in range(12)}
    # wel = flopy.mf6.ModflowGwfwel(gwf, stress_period_data=stress_period_data)

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

        rch_data = (
        i, cellid, 100, 5, 0.02, 99, 0.5, kh, 0.03, ncon, ustrf, ndiv)
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
        i: [(0, "INFLOW", 50), (6, "DIVERSION", 0, 10)] for i in range(12)
    }

    sfr = flopy.mf6.ModflowGwfsfr(gwf,
                                  nreaches=nreaches,
                                  packagedata=package_data,
                                  connectiondata=connection_data,
                                  diversions=diversions,
                                  perioddata=period_data,
                                  unit_conversion=86400,
                                  save_flows=True)

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
                    round_to_n(df.ppt_avg_m.values[i]/perlen[i], 5),
                    round_to_n(df.eto_avg_m.values[i]/perlen[i], 5),
                    5,
                    0.06,
                    -1.1,
                    -75.0,
                    1.0
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
        perioddata=period_data,
        unsat_etwc=True,
        linear_gwet=True,
        simulate_gwseep=True
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
    model_ws = gwf.model_ws
    sfr_name = gwf.sfr.filename
    uzf_name = gwf.uzf.filename
    mf6_dev_no_final_check(model_ws, sfr_name)
    mf6_dev_no_final_check(model_ws, uzf_name)
    return sim, gwf, ag


def mf6_dev_no_final_check(model_ws, fname):
    contents = []
    with open(os.path.join(model_ws, fname)) as foo:
        for line in foo:
            if "options" in line.lower():
                contents.append(line)
                contents.append("  DEV_NO_FINAL_CHECK\n")
            else:
                contents.append(line)

    with open(os.path.join(model_ws, fname), "w") as foo:
        for line in contents:
            foo.write(line)


def build_nwt(name):
    model_ws = os.path.join(sws, "..", "data", "nwt_test_ag_diversions")

    ml = flopy.modflow.Modflow(
        name,
        version="mfnwt",
        exe_name="mfnwt.exe",
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
    # build a SFR package
    reach_data = flopy.modflow.ModflowSfr2.get_empty_reach_data(11)
    for i in range(0, 10):
        outreach = i + 2
        if outreach == 11:
            outreach = 0
        reach_data[i] = (
        10 * i + 7, 0, i, 7, i + 1, 1, 100, 99, 0.02, 0.5, 0.000015, 0.35, 0.08, 5, 8.64, i + 1, outreach
        )
    reach_data[10] = (10 * 6 + 6, 0, 6, 6, 11, 1, 100, 99, 0.02, 0.5, 0.0000015, 0.35, 0.08, 5, 8.64, 11, 0)

    segdata = {}
    for per in range(12):
        segment_data = flopy.modflow.ModflowSfr2.get_empty_segment_data(11)
        for i in range(1, 11):
            outseg = i + 1
            if outseg == 11:
                outseg = 0
            flow = 0
            if i == 1:
                flow = 50
            segment_data[i - 1] = (
            i, 1, outseg, 0, 0, 0, flow, 0, 0, 0, 0.03, 0, 0, 0, 0, 0, 0.000015,
            0.5, 99 - ((i -1) * 0.02), 5, 0, 0, 0, 0, 0, 0.000015, 0.5,  99 - (i * 0.02), 5, 0, 0, 0, 0, 0)

        divflow = 10
        segment_data[10] = (
            11, 1, 0, 6, 0, 0, divflow, 0, 0, 0, 0.03, 0, 0, 0, 0, 0, 0.000015,
            0.5, 98.88, 5, 0, 0, 0, 0, 0, 0.000015, 0.5, 98.87, 5, 0, 0, 0, 0, 0)

        segdata[per] = segment_data

    sfr = flopy.modflow.ModflowSfr2(
        ml,
        nstrm=11,
        nss=11,
        const=86400,
        dleak=0.0001,
        ipakcb=52,
        reach_data=reach_data,
        segment_data=segdata,

    )

    # build a UZF package
    cimis_data = os.path.join("..", "data", "davis_monthly_ppt_eto.txt")
    df = pd.read_csv(cimis_data)

    finf = {}
    pet = {}
    extdp = {}
    extwc = {}
    for i in range(12):
        finf[i] = np.ones((nrow, ncol)) * round_to_n(
            df.ppt_avg_m.values[i] / perlen[i], 5)
        pet[i] = np.ones((nrow, ncol)) * round_to_n(
            df.eto_avg_m.values[i] / perlen[i], 5)
        extdp[i] = np.ones((nrow, ncol)) * 5
        extwc[i] = np.ones((nrow, ncol)) * 0.06

    uzf = flopy.modflow.ModflowUzf1(
        ml,
        nuztop=1,
        iuzfopt=1,
        irunflg=0,
        ipakcb=52,
        ietflg=1,
        ntrail2=7,
        nsets=40,
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
        ml,
        headtol=2e-03,
        fluxtol=50,
        thickfact=1e-08,
        options="SIMPLE"
    )

    gage = flopy.modflow.ModflowGage(
        ml,
        numgage=2,
        gage_data=[(6, 1, 55, 0), (11, 1, 56, 5)],
        files=["gage6.txt", "gage11.txt"]
    )

    ag = build_ag_package(name, ml, model_ws)
    ml.write_input()
    nam = os.path.join(model_ws, f"{name}.nam")
    with open(nam, "a") as foo:
        foo.write(f"DATA              53  {name}.diversion.txt")

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

        # if this is modflow 6 set segid to diversion id
        segid = 0
    else:
        dis = ml.dis
        segid = 11

    if name == "specified":
        options = flopy.utils.OptionBlock(
            "IRRIGATION_DIVERSION 1 2 "
            "TIMESERIES_DIVERSION".lower(),
            flopy.modflow.ModflowAg
        )
    elif name == "trigger":
        options = flopy.utils.OptionBlock(
            f"TRIGGER IRRIGATION_DIVERSION 1 2 "
            "TIMESERIES_DIVERSION".lower(),
            flopy.modflow.ModflowAg
        )
    else:
        options = flopy.utils.OptionBlock(
            f"{name} 1.0 IRRIGATION_DIVERSION 1 2 "
            "TIMESERIES_DIVERSION".lower(),
            flopy.modflow.ModflowAg
        )

    timeseries_diversion = flopy.modflow.ModflowAg.get_empty(1, block='time series')
    timeseries_diversion[0] = ("diversion", segid, 53)

    irrdiversion = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(1, 2, "irrdiversion")
        spd[0] = (segid, 2, 7, 0.2, 4, 4, 0, 0.5, 5, 4, 0, 0.5)
        irrdiversion[i] = spd

    # well_list = flopy.modflow.ModflowAg.get_empty(1, block="well")
    # x = [[0, 1, 1, -50.], ]
    # for ix, rec in enumerate(well_list):
    #     well_list[ix] = tuple(x[ix])

    # irrwell = {}
    # for i in range(12):
    #     spd = flopy.modflow.ModflowAg.get_empty(1, 2, "irrwell")
    #     spd[0] = (0, 2, 31, 0.2, 2, 1, 1.0, 0.5, 2, 2, 1.0, 0.5)
    #     irrwell[i] = spd

    ag = flopy.modflow.ModflowAg(
        ml,
        options=options,
        time_series=timeseries_diversion,
        # well_list=well_list,
        # irrwell=irrwell,
        irrdiversion=irrdiversion,
        nper=dis.nper
    )

    return ag


def run_mfnwt(name):
    ml, ag = build_nwt(name)
    success, buff = ml.run_model()
    if not success:
        raise AssertionError


def run_mf6(name):
    dll = os.path.join("..", "modflow-bmi", "libmf6.dll")
    sim, gwf, ag = build_mf6(name)
    mf6ag = Modflow6Ag(sim, ag)
    mf6ag.run_model(dll)


def compare_output(name):
    lst_name = "{}.lst".format(name)
    csv_name = "{}.diversion.txt".format(name)
    nwt_ws = os.path.join(sws, "..", "data", "nwt_test_ag_diversions")
    mf6_ws = os.path.join(sws, "..", "data", "mf6_test_ag_diversions")

    nwt_div = os.path.join(nwt_ws, csv_name)
    mf6_lst = os.path.join(mf6_ws, lst_name)

    nwt_div = pd.read_csv(nwt_div, delim_whitespace=True)
    nwt_dis = nwt_div["SW-DIVERSION"].sum()

    # mf6_uzf_dis = mf6_cbc.get_data(text="UZF-GWET")
    # create utility to look at sfr diversions from
    mf6_div = read_sfr_list_data(mf6_lst)
    mf6_div = mf6_div[mf6_div.iseg == 11]
    mf6_dis = mf6_div["inflow"].sum()

    print(nwt_dis, mf6_dis)

    # x = pd.read_csv('factor.txt')
    # x = x.groupby(by="kstp", as_index=False)["factor"].max()
    # x['mf6'] = l2
    # x['mfnwt'] = l
    # x = x[x["mf6"] < 0]

    plt.plot(range(1, len(mf6_div) + 1), mf6_div['inflow'].values, "r-", label='mf6 ag')
    plt.plot(range(1, len(nwt_div) + 1), nwt_div["SW-DIVERSION"].values, "b--", label="mfnwt ag")
    # plt.plot(range(1, len(d) + 1), d2, "g-", label='mf6 uzf discharge')
    # plt.plot(range(1, len(d) + 1), d, "k--", label="mfnwt uzf discharge")
    # plt.plot(range(1, 366), [0.2] * 365, 'c--')
    # plt.plot(x.kstp.values, x.factor.values * -1, "k-", label="aet/pet")
    plt.legend(loc=0)
    #
    plt.show()


def heads(name):
    name = f"{name}.hds"
    head_file = os.path.join("..", "data", "nwt_test_ag_diversions", name)

    hds = flopy.utils.HeadFile(head_file)
    data = hds.get_alldata()


def compare_etdemand():
    nwt_debug = os.path.join("..", "data", "nwt_test_ag_diversions", "debug_etdemand.out")
    mf6_debug = os.path.join("etdiv_factor.txt")

    nwt_header = ["kper", "kstp", "iseg", "crop_pet", "crop_aet", "sup", "supold", "factor"]

    dfnwt = pd.read_csv(nwt_debug, names=nwt_header, delim_whitespace=True)
    dfmf6 = pd.read_csv(mf6_debug)
    dfmf6.kper += 1
    dfmf6.kstp += 1
    dfnwt.to_csv("etdiv_factor_nwt.txt", index=False)
    print('break')


if __name__ == "__main__":
    debug = False
    if not debug:
        try:
            os.remove(os.path.join("..", "data", "nwt_test_ag_diversions", "debug.out"))
        except:
            pass
        try:
            os.remove(os.path.join("..", "data", "nwt_test_ag_diversions", "debug_trigger.out"))
        except:
            pass
        try:
            os.remove(os.path.join("..", "data", "nwt_test_ag_diversions", "debug_etdemand.out"))
        except:
            pass
        try:
            os.remove(os.path.join("div_factor.txt"))
        except:
            pass
        with open(os.path.join("div_factor.txt"), "w") as foo:
            foo.write("factor,kstp,kiter,crop_aet,crop_pet\n")
        try:
            os.remove(os.path.join("etdiv_factor.txt"))
        except:
            pass
        with open(os.path.join("etdiv_factor.txt"), "w") as foo:
            foo.write("kper,kstp,crop_pet,crop_aet,sup,supold,factor\n")

        model_names = ("etdemand", "trigger", "specified")
        model = 0
        run_mfnwt(model_names[model])
        run_mf6(model_names[model])
        compare_output(model_names[model])
        heads(model_names[model])
    compare_etdemand()