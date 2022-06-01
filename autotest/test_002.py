# test a simple 2 SFR diversion, 2 farm (4 receiver cells) problem
import os
import sys
import pandas as pd
import numpy as np
import flopy
sys.path.append(os.path.join("..", "mf6api_agmvr"))
from mf6_agmvr import ModflowAgmvr
from get_mvr_budget import MvrBudget
import common


def test_etdemand_sfr():
    name = "etdemand_div"
    sim_ws = os.path.join(".", "temp", "etdemand_div")
    sim = flopy.mf6.MFSimulation(name, sim_ws=sim_ws)

    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    period_data = [(i, i, 1.0) for i in perlen]
    tdis = flopy.mf6.ModflowTdis(
        sim,
        nper=12,
        perioddata=tuple(period_data),
        time_units="days"
    )

    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        complexity="COMPLEX",
        no_ptcrecord=["ALL"],
        outer_dvclose=None, # 5e-02,
        outer_maximum=None, # 213.16,
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
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True,
                                  icelltype=1)
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
        i: [(0, "INFLOW", 50), ] for i in range(12)
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
                    common.round_to_n(df.ppt_avg_m.values[i] / perlen[i], 5),
                    common.round_to_n(df.eto_avg_m.values[i] / perlen[i], 5),
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

    mfag = ModflowAgmvr(sim, ag_type="etdemand", mvr_name="mvr")
    mfag.run_model(common.dll_loc())

    nwt_div1 = os.path.join(common.nwt_output_path(), f"{name}.diversion11.txt")
    nwt_div2 = os.path.join(common.nwt_output_path(), f"{name}.diversion12.txt")
    mf6_lst = os.path.join(sim_ws, f"{name}.lst")

    nwt_div1 = pd.read_csv(nwt_div1, delim_whitespace=True)
    nwt_div2 = pd.read_csv(nwt_div2, delim_whitespace=True)

    mf6_div = MvrBudget(mf6_lst).inc
    mf6_div = mf6_div.groupby(by=["provider", "pid", "totim"], as_index=False)[
        ["qa", "qp"]].sum()
    mf6_div1 = mf6_div[mf6_div.pid == 2]
    mf6_div2 = mf6_div[mf6_div.pid == 4]

    nwt_total = nwt_div1["SW-DIVERSION"].values + nwt_div2[
        "SW-DIVERSION"].values
    mf6_total = mf6_div1.qp.values + mf6_div2.qp.values

    err = np.sum(mf6_total - nwt_total)
    if np.abs(err) > 0.1:
        raise AssertionError("Diversion calculation outside of acceptable tolerance")


if __name__ == "__main__":
    test_etdemand_sfr()