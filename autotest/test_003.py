# test a simple 2 SFR diversion + 2 supplemental pumping well problem
# that irrigates 2 "farms" (2 cells each "farm").
import flopy
import numpy as np
import pandas as pd
import sys
import os
from mf6api_ag import ModflowApiAg
import common



def test_conjunctive_sfr_wel():
    name = "etdemand_sup"
    sim_ws = os.path.join(".", "temp", "etdemand_sup")
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

    # create well package
    stress_period_data = {}
    for i in range(12):
        if i == 2:
            stress_period_data[i] = [[(0, 1, 2), -50.],
                                     [(0, 5, 4), -100.], ]
        else:
            stress_period_data[i] = [[(0, 1, 2), -50.],
                                     [(0, 5, 4), -100.], ]

    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=stress_period_data,
        mover=True
    )

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
        "etdemand_sup.obs.out": [["sw_outflow", "outflow", 9]],
        "filename": "etdemand_sup.obs"
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

    # create mvr package
    period_data = {}
    for i in range(12):
        mvr_rec = []
        for col in range(23, 25):
            sfrrec = ("sfr_0", 2, "uzf_0", col, "UPTO", 10.)
            wellrec = ("wel_0", 0, "uzf_0", col, "UPTO", 20)
            mvr_rec.append(sfrrec)
            mvr_rec.append(wellrec)
        for col in range(43, 45):
            sfrrec = ("sfr_0", 4, "uzf_0", col, "UPTO", 5.)
            wellrec = ("wel_0", 1, "uzf_0", col, "UPTO", 20)
            mvr_rec.append(sfrrec)
            mvr_rec.append(wellrec)

        period_data[i] = mvr_rec

    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=8,
        maxpackages=3,
        packages=[("sfr_0",), ("uzf_0",), ("wel_0",)],
        perioddata=period_data,
        budgetcsv_filerecord=f"{name}_mvr.csv",
        budget_filerecord=f"{name}_mvr.cbc"
    )

    sim.write_simulation()
    mfag = ModflowApiAg(sim, ag_type="etdemand", mvr_name="mvr")
    mfag.run_model(common.dll_loc())

    nwt_div1 = os.path.join(common.nwt_output_path(), f"{name}.diversion11.txt")
    nwt_div2 = os.path.join(common.nwt_output_path(), f"{name}.diversion12.txt")
    nwt_cbc = os.path.join(common.nwt_output_path(), f"{name}.cbc")

    mf6_ag_out = os.path.join(sim_ws, f"{name}_ag.out")
    mf6_ag_out = ModflowApiAg.load_output(mf6_ag_out)
    mf6_ag_out = mf6_ag_out.groupby(by=["pkg", "pid", "kstp"], as_index=False)[["q_from_provider", "q_to_receiver"]].sum()
    mf6_div1 = mf6_ag_out[mf6_ag_out.pid == 2]
    mf6_div2 = mf6_ag_out[mf6_ag_out.pid == 4]
    mf6_well = mf6_ag_out[mf6_ag_out.pid == 0]
    mf6_well2 = mf6_ag_out[mf6_ag_out.pid == 1]

    nwt_div1 = pd.read_csv(nwt_div1, delim_whitespace=True)
    nwt_div2 = pd.read_csv(nwt_div2, delim_whitespace=True)

    nwt_cbc = flopy.utils.CellBudgetFile(nwt_cbc)
    nwt_pump = nwt_cbc.get_data(text="AG WE")

    nwt_well = []
    nwt_well2 = []
    for ix, recarray in enumerate(nwt_pump):
        idx = np.where(nwt_pump[ix]["node"] == 13)[0]
        nwt_well.append(nwt_pump[ix][idx]["q"][0])
        idx = np.where(nwt_pump[ix]["node"] == 55)[0]
        nwt_well2.append(nwt_pump[ix][idx]["q"][0])

    total_ag1 = mf6_div1.q_from_provider.values + mf6_well.q_from_provider.values
    total_ag2 = mf6_div2.q_from_provider.values + mf6_well2.q_from_provider.values
    nwt_total1 = nwt_div1["SW-DIVERSION"].values + np.abs(nwt_well)
    nwt_total2 = nwt_div2["SW-DIVERSION"].values + np.abs(nwt_well2)

    err1 = np.sum(total_ag1 - nwt_total1)
    err2 = np.sum(total_ag2 - nwt_total2)

    if np.abs(err1) > 0.1:
        raise AssertionError(
            "supplemental pumping calculation out of acceptable tolerance"
        )

    if np.abs(err2) > 0.1:
        raise AssertionError(
            "supplemental pumping calculation out of acceptable tolerance"
        )


if __name__ == "__main__":
    test_conjunctive_sfr_wel()