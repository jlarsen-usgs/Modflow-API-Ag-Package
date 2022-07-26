# test 2 MAW pumping wells problem that irrigates 2 "farms"
# (2 cells each "farm").
import flopy
import numpy as np
import pandas as pd
import sys
import os
from mf6api_agmvr import ModflowAgmvr
import common


def test_etdemand_maw():
    name = "etdemand_maw"
    sim_ws = os.path.join(".", "temp", "etdemand_maw")
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
    npf = flopy.mf6.ModflowGwfnpf(gwf, save_specific_discharge=True,
                                  icelltype=1)
    sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1)

    nmawwells = 2
    package_data = [
        (0, 0.5, 50, 95, "theim", 1),
        (1, 0.5, 75, 95, "theim", 1)
    ]
    connection_data = [
        (0, 0, (0, 1, 2), 95, 50, 0.1, 0.25),
        (1, 0, (0, 5, 4), 95, 75, 0.1, 0.25)
    ]

    period_data = {i: [(0, "rate", -50), (1, "rate", -100)] for i in range(12)}

    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        save_flows=True,
        mover=True,
        nmawwells=nmawwells,
        packagedata=package_data,
        connectiondata=connection_data,
        perioddata=period_data
    )

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
        for col in range(43, 45):
            rec = ("maw_0", 1, "uzf_0", col, "UPTO", 50.)
            mvr_rec.append(rec)
        for col in range(21, 23):
            rec = ("maw_0", 0, "uzf_0", col, "UPTO", 25.)
            mvr_rec.append(rec)
        period_data[i] = mvr_rec

    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=4,
        maxpackages=2,
        packages=[("maw_0",), ("uzf_0",)],
        perioddata=period_data
    )

    sim.write_simulation()
    model_ws = gwf.model_ws
    uzf_name = gwf.uzf.filename
    common.mf6_dev_no_final_check(model_ws, uzf_name)

    mfag = ModflowAgmvr(sim, ag_type="etdemand", mvr_name="mvr")
    mfag.run_model(common.dll_loc())

    # compare output to the mfnwt well pumping problem
    mf6_ag_out = os.path.join(model_ws, f"{name}_ag.out")
    mf6_ag_out = ModflowAgmvr.load_output(mf6_ag_out)
    mf6_ag_out = mf6_ag_out.groupby(by=["pkg", "pid", "kstp"], as_index=False)[["q_from_provider", "q_to_receiver"]].sum()

    nwt_cbc = os.path.join(common.nwt_output_path(), f"etdemand_well.cbc")
    nwt_cbc = flopy.utils.CellBudgetFile(nwt_cbc)

    nwt_pump = nwt_cbc.get_data(text="AG WE")

    mf6_total = mf6_ag_out.q_from_provider.values
    nwt_total = []
    for wl in (13, 55):
        for ix, recarray in enumerate(nwt_pump):
            idx = np.where(recarray["node"] == wl)[0]
            nwt_total.append(recarray[idx]['q'][0])

    err = np.sum(np.abs(mf6_total) - np.abs(nwt_total))
    if np.abs(err) > 0.1:
        raise AssertionError("Error greater than acceptable tolerance")


if __name__ == "__main__":
    test_etdemand_maw()