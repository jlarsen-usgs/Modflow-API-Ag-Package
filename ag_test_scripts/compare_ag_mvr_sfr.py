import flopy
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sws = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(sws, "..", "develop_AG_mvr"))
from mf6_ag_mvr import ModflowAgmvr

from math import log10, floor


def round_to_n(x, n):
    if x == 0:
        return 0
    t = round(x, -int(floor(log10(abs(x))) - (n - 1)))
    return t


def build_mf6(name, headtol=None, fluxtol=None):
    sim_ws = os.path.join(sws, "..", "data", "mf6_etdemand_test_problems")
    sim = flopy.mf6.MFSimulation(name, sim_ws=sim_ws)

    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    period_data = [(i, i, 1.0) for i in perlen]
    tdis = flopy.mf6.ModflowTdis(
        sim,
        nper=12,
        perioddata=tuple(period_data),
        time_units="days"
    )

    if headtol is None:
        if name == "etdemand":
            headtol = 0.0570641530019691
        elif name == "trigger":
            headtol = 0.0570641530019691

    if fluxtol is None:
        if name == "etdemand":
            fluxtol = 213.1677138100136
        elif name == "trigger":
            fluxtol = 213.1677138100136

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        complexity="COMPLEX",
        no_ptcrecord=["ALL"],
        outer_dvclose=headtol,
        outer_maximum=fluxtol,
        rcloserecord=[1e-10, "L2NORM_RCLOSE"],
        scaling_method="L2NORM",
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

    # create SFR package
    nreaches = 10
    package_data = []
    for i in range(nreaches):
        ustrf = 1.0
        if i in (0, 9):
            ncon = 1
        else:
            ncon = 2
        ndiv = 0
        cellid = (0, i, 7)
        kh = 0.000015
        rch_data = (
            i, cellid, 100, 5, 0.02, 99, 0.5, kh, 0.03, ncon, ustrf, ndiv
        )
        package_data.append(rch_data)

    connection_data = []
    for i in range(nreaches):
        if i == 0:
            cd = [i, -1 * (i + 1)]
        elif i == 9:
            cd = [i, i - 1]
        else:
            cd = [i, i - 1, -1 * (i + 1)]
        connection_data.append(cd)

    period_data = {
        i: [(0, "INFLOW", 50),] for i in range(12)
    }

    sfr = flopy.mf6.ModflowGwfsfr(gwf,
                                  nreaches=nreaches,
                                  packagedata=package_data,
                                  connectiondata=connection_data,
                                  perioddata=period_data,
                                  unit_conversion=86400,
                                  save_flows=True,
                                  mover=True)


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
                    4,
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
        simulate_gwseep=True,
        mover=True
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

    # create mvr package for wells
    period_data = {}
    for i in range(12):
        mvr_rec = []
        for col in range(23, 25):
            rec = ("sfr_0", 2, "uzf_0", col, "UPTO", 10.)
            mvr_rec.append(rec)
        for col in range(43, 45):
            rec = ("sfr_0", 4, "uzf_0", col, "UPTO", 10.)
            mvr_rec.append(rec)
        period_data[i] = mvr_rec

    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=4,
        maxpackages=2,
        packages=[("sfr_0",), ("uzf_0",)],
        perioddata=period_data
    )

    sim.write_simulation()
    model_ws = gwf.model_ws
    uzf_name = gwf.uzf.filename
    sfr_name = gwf.sfr.filename
    mf6_dev_no_final_check(model_ws, uzf_name)
    mf6_dev_no_final_check(model_ws, sfr_name)
    return sim, gwf


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


def compare_model_output(nwt, mf6, model):
    nwt_div1 = os.path.join(nwt, f"{model}.diversion11.txt")
    nwt_div2 = os.path.join(nwt, f"{model}.diversion12.txt")
    mf6_cbc = os.path.join(mf6, f"{model}.cbc")

    nwt_div1 = pd.read_csv(nwt_div1, delim_whitespace=True)
    nwt_div2 = pd.read_csv(nwt_div2, delim_whitespace=True)

    mf6_cbc = flopy.utils.CellBudgetFile(mf6_cbc)

    mf6_pump = mf6_cbc.get_data(text="WEL-TO-MVR")
    mf6_pump2 = mf6_cbc.get_data(text="WEL")

    fig, ax = plt.subplots(figsize=(8, 8))
    for wl in (45, 23):
        nwt_well = []
        mf6_well = []
        mf6_well2 = []
        for ix, recarray in enumerate(mf6_pump):
            idx = np.where(recarray["node"] == wl)[0]
            nwt_well.append(recarray[idx]['q'])
            idx = np.where(mf6_pump[ix]["node"] == wl)[0]
            mf6_well.append(mf6_pump[ix][idx]["q"])
            idx = np.where(mf6_pump2[ix]["node"] == wl)[0]
            mf6_well2.append(mf6_pump2[ix][idx]["q"])
        ax.plot(range(1, len(nwt_well) + 1), nwt_well, label=f"nwt well, node {wl}", lw=2)
        ax.plot(range(1, len(mf6_well) + 1), mf6_well, label=f"mf6 well to mvr, node {wl}", ls="--")
        ax.plot(range(1, len(mf6_well) + 1), mf6_well2, label=f"mf6 well, node {wl}", ls="-.")
        ax.plot()
        print("MF6: ", np.sum(mf6_well) * 0.000810714)
        print("NWT: ", np.sum(np.abs(nwt_well)) * 0.000810714)

    plt.legend(loc=0)
    plt.show()


if __name__ == "__main__":
    load_existing = False
    run_model = True
    dll = os.path.join("..", "modflow-bmi", "libmf6.dll")
    mf6_ws = os.path.join(sws, "..", "data", "mf6_etdemand_test_problems")
    nwt_ws = os.path.join(sws, "..", "data", "nwt_etdemand_test_problems")
    if run_model:
        if not load_existing:
            sim, gwf = build_mf6("etdemand_div")
        else:
            sim = flopy.mf6.MFSimulation.load(
                sim_ws=os.path.join(sws, "..", "data", "mf6_etdemand_test_problems")
            )
        mfag = ModflowAgmvr(sim, ag_type="etdemand", mvr_name="mvr")
        mfag.run_model(dll)

    compare_model_output(nwt_ws, mf6_ws, "etdemand_div")
