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
sys.path.append(os.path.join(sws, "..", "mf6api_agmvr"))
from mf6_agmvr import ModflowAgmvr
mpl.use("TKagg")

from math import log10, floor


def round_to_n(x, n):
    if x == 0:
        return 0
    t = round(x, -int(floor(log10(abs(x))) - (n - 1)))
    return t


def build_mf6(name, headtol=None, fluxtol=None):
    sim_ws = os.path.join(sws, "..", "data", "mf6_lak_etdemand_test_problem")
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
        outer_dvclose=None,
        outer_maximum=None,
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

    package_data = []
    for lk in range(3):
        rec = (lk, 100, 1,)
        package_data.append(rec)

    connection_data = []
    for lk, (node, gwf_cell) in enumerate(([2, 27], [4, 47], [6, 67])):
        cellid = (0, node, 7)
        lak_data = [(lk, 0, cellid, "vertical", 1e-06, 90, 100, 63.6, 63.6),]

        connection_data += lak_data

    outlets = []
    for i in range(3):
        rec = (i, i, -1, "specified", 0, 0, 0, 0)
        outlets.append(rec)

    period_data = {}
    for per in range(12):
        spd = []
        for lk in range(2):
            rec1 = (
                lk,
                "rate",
                -50 * (lk + 1)
            )
            rec2 = (lk,
                     "status",
                     "constant")
            rec3 = (
                lk,
                "stage",
                110
            )

            # spd.append(rec)
            # spd.append(rec2)
            spd.append(rec1)
            spd.append(rec2)
            spd.append(rec3)
        # if per == 0:
        #     spd.append((2, "status", "constant"))
        #     spd.append((2, "stage", 105))

        period_data[per] = spd

    lak = flopy.mf6.ModflowGwflak(
        gwf,
        save_flows=True,
        print_stage=True,
        mover=True,
        nlakes=3,
        noutlets=3,
        ntables=0,
        packagedata=package_data,
        connectiondata=connection_data,
        outlets=outlets,
        perioddata=period_data
    )

    # build a UZF package
    cimis_data = os.path.join("..", "data", "davis_monthly_ppt_eto.txt")
    df = pd.read_csv(cimis_data)

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
        for col in range(43, 45):
            rec = ("lak_0", 1, "uzf_0", col, "UPTO", 25.)
            mvr_rec.append(rec)
        for col in range(21, 23):
            rec = ("lak_0", 0, "uzf_0", col, "UPTO", 25.)
            mvr_rec.append(rec)
        period_data[i] = mvr_rec

    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=4,
        maxpackages=2,
        packages=[("lak_0",), ("uzf_0",)],
        perioddata=period_data
    )

    sim.write_simulation()
    model_ws = gwf.model_ws
    uzf_name = gwf.uzf.filename
    # lak_name = gwf.lak.filename
    mf6_dev_no_final_check(model_ws, uzf_name)
    # mf6_dev_no_final_check(model_ws, lak_name)
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


def run_mf6_exe(fpsim):
    fpsim.run_simulation()
    return


def compare_model_output(nwt, mf6, model):

    mf6_cbc = os.path.join(mf6, f"{model}.cbc")
    nwt_cbc = os.path.join(nwt, "etdemand_well.cbc")


    nwt_cbc = flopy.utils.CellBudgetFile(nwt_cbc)
    mf6_cbc = flopy.utils.CellBudgetFile(mf6_cbc)

    nwt_pump = nwt_cbc.get_data(text="AG WE")
    mf6_pump = mf6_cbc.get_data(text="MAW")

    mf6_total = []
    nwt_total = []
    with styles.USGSPlot():
        mpl.rcParams["ytick.labelsize"] = 9
        mpl.rcParams["xtick.labelsize"] = 9
        fig, ax = plt.subplots(figsize=(8, 8))
        n = 1
        nwt_c = ["k", "yellow"]
        mf6_c = ["skyblue", "darkblue"]
        for wl in (13, 55):
            nwt_well = []
            mf6_well = []
            for ix, recarray in enumerate(nwt_pump):
                idx = np.where(recarray["node"] == wl)[0]
                nwt_well.append(recarray[idx]['q'])
                nwt_total.append(recarray[idx]['q'][0])
                idx = np.where(mf6_pump[ix]["node"] == wl)[0]
                mf6_well.append(mf6_pump[ix][idx]["q"])
                mf6_total.append(mf6_pump[ix][idx]["q"][0])

            ax.plot(range(1, len(nwt_well) + 1), np.abs(nwt_well), color=nwt_c[n - 1], label=f"nwt well {n} irrigation", lw=2)
            ax.plot(range(1, len(mf6_well) + 1), np.abs(mf6_well), color=mf6_c[n - 1], label=f"mf6 MAW {n} irrigation", lw=2.5, ls="--")
            n += 1

            print("MF6: ", np.sum(np.abs(mf6_well)) * 0.000810714)
            print("NWT: ", np.sum(np.abs(nwt_well)) * 0.000810714)

        r2 = linregress(np.abs(nwt_total), np.abs(mf6_total))[2] ** 2
        err = np.sum(np.abs(mf6_total) - np.abs(nwt_total))

        styles.heading(ax=ax,
                       heading="Comparison of MF6 API AG MAW and MF-NWT AG well pumping")
        styles.xlabel(ax=ax, label="Days", fontsize=10)
        styles.ylabel(ax=ax, label="Applied irrigation, in " + r"$m^{3}$",
                      fontsize=10)
        styles.graph_legend(ax=ax, loc=2, fancybox=True, shadow=True, frameon=True, fontsize=10)
        styles.add_text(ax=ax, text=r"$r^{2}$" + f" = {r2 :.2f}", x=0.03,
                        y=0.78, fontsize=10)
        styles.add_text(ax=ax, text=f"dif. = {err :.2f} " + r"$m^{3}$", x=0.03,
                        y=0.75, fontsize=10)
        plt.show()


def inspect_model_output(model_ws, name):
    mf6_cbc = os.path.join(model_ws, f"{name}.cbc")
    mf6_cbc = flopy.utils.CellBudgetFile(mf6_cbc)
    mf6_pump = mf6_cbc.get_data(text="MAW")

    with styles.USGSPlot():
        mpl.rcParams["ytick.labelsize"] = 9
        mpl.rcParams["xtick.labelsize"] = 9
        fig, ax = plt.subplots(figsize=(8, 8))
        n = 1
        nwt_c = ["k", "yellow"]
        mf6_c = ["skyblue", "darkblue"]
        for wl in (13, 55):
            mf6_well = []
            for ix, recarray in enumerate(mf6_pump):
                idx = np.where(recarray["node"] == wl)[0]
                mf6_well.append(recarray[idx]['q'][0])

            ax.plot(range(1, len(mf6_well) + 1), np.abs(mf6_well),
                    color=mf6_c[n - 1], label=f"mf6 MAW {n} irrigation",
                    lw=2.5, ls="--")
            n += 1

        styles.heading(ax=ax,
                       heading="MF6 API AG multi-aquifer well pumping")
        styles.xlabel(ax=ax, label="Days", fontsize=10)
        styles.ylabel(ax=ax, label="Applied irrigation, in " + r"$m^{3}$",
                      fontsize=10)
        styles.graph_legend(ax=ax, loc=2, fancybox=True, shadow=True,
                            frameon=True, fontsize=10)
        plt.show()


if __name__ == "__main__":
    # set dll path
    load_existing = False
    run_model = True
    dll = os.path.join("..", "modflow-bmi", "libmf6.dll")
    mf6_ws = os.path.join(sws, "..", "data", "mf6_lak_etdemand_test_problem")
    nwt_ws = os.path.join(sws, "..", "data", "nwt_etdemand_test_problems")
    model_name = "etdemand_lak"
    if run_model:
        if not load_existing:
            sim, gwf = build_mf6(model_name)
        else:
            sim = flopy.mf6.MFSimulation.load(
                sim_ws=os.path.join(sws, "..", "data", "mf6_lak_etdemand_test_problem")
            )
        mfag = ModflowAgmvr(sim, ag_type="etdemand", mvr_name="mvr")
        mfag.run_model(dll)

    inspect_model_output(mf6_ws, model_name)
    compare_model_output(nwt_ws, mf6_ws, model_name)