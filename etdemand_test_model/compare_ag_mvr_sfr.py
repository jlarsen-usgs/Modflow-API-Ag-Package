import flopy
from flopy.plot import styles
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress
sws = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(sws, "..", "develop_AG_mvr"))
from mf6_agmvr import ModflowAgmvr
from get_mvr_budget import MvrBudget

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
        outer_dvclose=5e-02,
        outer_maximum=213.16,
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

    obs_dict = {
        "etdemand_div.obs.out": [["sw_outflow", "outflow", 9]],
        "filename": "etdemand_div.obs"
    }
    sfr = flopy.mf6.ModflowGwfsfr(gwf,
                                  nreaches=nreaches,
                                  packagedata=package_data,
                                  connectiondata=connection_data,
                                  perioddata=period_data,
                                  unit_conversion=86400,
                                  save_flows=True,
                                  mover=True,
                                  observations=obs_dict)

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
            rec = ("sfr_0", 4, "uzf_0", col, "UPTO", 5.)
            mvr_rec.append(rec)
        period_data[i] = mvr_rec

    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=4,
        maxpackages=2,
        packages=[("sfr_0",), ("uzf_0",)],
        perioddata=period_data,
        budgetcsv_filerecord=f"{name}_mvr.csv",
        budget_filerecord=f"{name}_mvr.cbc"
    )

    sim.write_simulation()
    model_ws = gwf.model_ws
    uzf_name = gwf.uzf.filename
    sfr_name = gwf.sfr.filename
    # mf6_dev_no_final_check(model_ws, uzf_name)
    # mf6_dev_no_final_check(model_ws, sfr_name)
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
    mf6_lst = os.path.join(mf6, f"{model}.lst")
    mf6_cbc = os.path.join(mf6, f"etdemand_well.cbc")

    nwt_div1 = pd.read_csv(nwt_div1, delim_whitespace=True)
    nwt_div2 = pd.read_csv(nwt_div2, delim_whitespace=True)

    mf6_cbc = flopy.utils.CellBudgetFile(mf6_cbc)
    mf6_pump = mf6_cbc.get_data(text="WEL-TO-MVR")
    well_1 = []
    for recarray in mf6_pump:
        well_1.append(recarray[0]["q"])

    mf6_div = MvrBudget(mf6_lst).inc
    mf6_div = mf6_div.groupby(by=["provider", "pid", "totim"], as_index=False)[["qa", "qp"]].sum()
    mf6_div1 = mf6_div[mf6_div.pid == 2]
    mf6_div2 = mf6_div[mf6_div.pid == 4]

    nwt_total = nwt_div1["SW-DIVERSION"].values + nwt_div2["SW-DIVERSION"].values
    mf6_total = mf6_div1.qp.values + mf6_div2.qp.values

    r2 = linregress(mf6_total, np.abs(nwt_total))[2] ** 2
    err = np.sum(mf6_total - nwt_total)
    perr = err / np.sum(nwt_total)

    with styles.USGSPlot():
        mpl.rcParams["ytick.labelsize"] = 9
        mpl.rcParams["xtick.labelsize"] = 9
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.plot(range(1, len(nwt_div1) + 1), nwt_div1["SW-DIVERSION"], "k", label=f"nwt diversion, 1", lw=2)
        ax.plot(range(1, len(nwt_div2) + 1), nwt_div2["SW-DIVERSION"], "dimgray", label=f"nwt diversion, 2", lw=2)
        ax.plot(range(1, len(mf6_div1) + 1), mf6_div1.qp, color="skyblue",  label=f"mf6 diversion 1", ls="--", lw=2.5, zorder=3)
        ax.plot(range(1, len(mf6_div2) + 1), mf6_div2.qp, color="darkblue", label=f"mf6 diversion 2", ls="--", lw=2.5, zorder=3)
        styles.heading(ax=ax, heading="Comparison of MF6 API AG and MF-NWT AG diversions")
        styles.xlabel(ax=ax, label="Days", fontsize=10)
        styles.ylabel(ax=ax, label="Applied irrigation, in " + r"$m^{3}$", fontsize=10)
        styles.graph_legend(ax=ax, loc=2, fancybox=True, shadow=True, frameon=True, fontsize=10)
        styles.add_text(ax=ax, text=r"$r^{2}$" + f" = {r2 :.2f}", x=0.03,
                        y=0.78, fontsize=10)
        styles.add_text(ax=ax, text=f"dif. = {err :.2f} " + r"$m^{3}$", x=0.03,
                        y=0.75, fontsize=10)

        print("MF6: ", np.sum(mf6_div1.qp) * 0.000810714)
        print("NWT: ", np.sum(nwt_div1["SW-DIVERSION"]) * 0.000810714)
        print("MF6: ", np.sum(mf6_div2.qp) * 0.000810714)
        print("NWT: ", np.sum(nwt_div2["SW-DIVERSION"]) * 0.000810714)

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
    else:
        gwf = None
    compare_model_output(nwt_ws, mf6_ws, "etdemand_div")
