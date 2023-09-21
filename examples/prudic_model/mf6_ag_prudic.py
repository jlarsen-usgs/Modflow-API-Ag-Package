"""
Modflow 6 API version of the Green Valley example problem described in:

Prudic, D. E., Konikow, L. F., Banta, E. A. 2004. A New Streamflow-Routing
    (SFR1) Package to simulate stream-aquifer interaction with MODFLOW-2000:
    U.S. Geological Survey Tech. Methods 6-A13. DOI: 10.3133/ofr20041042

Niswonger, R. G., Prudic, D. E., Regan, R. S. 2006. Documentation of the
    Unsaturated-Zone Flow (UZF) Package for modeling unsaturated flow between
    the land surface and the water table with MODFLOW-2005. U.S. Geological
    Survey Techniques and Methods 6-A19. DOI: 10.3133/tm6A19

Niswonger, R. G. 2020. An agricultural water use package for MODFLOW and
    GSFLOW. Environmental Modelling and Software 125.
    doi: 10.1016/j.envsoft.2019.104617
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import flopy
import sys
from mf6api_ag import ModflowApiAg


def build_model(name, sim_ws):
    sim = flopy.mf6.MFSimulation(name, sim_ws=sim_ws)
    data_pth = os.path.join(".", "arrays")

    # build tdis package
    nper = 49
    perlen = [(1, 1, 1)]
    perlen += [(30.5, 30, 1) for _ in range(nper - 1)]
    tdis = flopy.mf6.ModflowTdis(
        sim,
        nper=nper,
        perioddata=tuple(perlen),
        time_units="days",
    )

    # build ims package
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        complexity="COMPLEX",
        # no_ptcrecord=["ALL"],
        outer_dvclose=None,
        outer_maximum=None,
        # rcloserecord=[1e-10, "L2NORM_RCLOSE"],
        # scaling_method="L2NORM",
        # linear_acceleration="BICGSTAB",
        # under_relaxation="DBD",
        # under_relaxation_gamma=0.0,
        # under_relaxation_theta=0.97,
        # under_relaxation_kappa=0.0001
    )

    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname=name,
        save_flows=True,
        print_input=True,
        print_flows=True,
        newtonoptions="NEWTON UNDER RELAXATION")

    # create dis package
    nlay, nrow, ncol = 1, 15, 10
    top = np.loadtxt(os.path.join(data_pth, "top.txt"), dtype=float)
    botm = np.loadtxt(os.path.join(data_pth, "bottom.txt"), dtype=float)
    idomain = np.loadtxt(os.path.join(data_pth, "idomain.txt"), dtype=int)
    dis = flopy.mf6.ModflowGwfdis(
        gwf,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=5000,
        delc=5000,
        top=top,
        botm=botm,
        idomain=idomain,
        length_units="feet"
    )

    # create ic package
    strt = np.loadtxt(os.path.join(data_pth, "strt.txt"), dtype=float)
    ic = flopy.mf6.ModflowGwfic(
        gwf,
        strt=strt
    )

    # create npf package
    hk = np.loadtxt(os.path.join(data_pth, "hk.txt"), dtype=float)
    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        save_flows=True,
        save_specific_discharge=True,
        k=hk,
        k33=1e-06
    )

    # create storage package
    sy = np.loadtxt(os.path.join(data_pth, "sy.txt"), dtype=float)
    steady_state = {0: True}
    transient = {0: False}
    for i in range(1, nper):
        steady_state[i] = False
        transient[i] = True

    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        save_flows=True,
        iconvert=1,
        ss=1e-06,
        sy=sy,
        steady_state=steady_state,
        transient=transient
    )

    # create ghb package
    period_data = {}
    for kper in range(49):
        ghb_spd = [
            ((0, 12, 0), 988.0, 0.038),
            [(0, 13, 8), 1045.0, 0.038],
        ]
        period_data[kper] = ghb_spd

    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        save_flows=True,
        stress_period_data=period_data
    )

    # create well package data
    period_data = {}
    for kper in range(nper):
        if kper in (0, 1, 2, 3):
            period_data[kper] = [
                ((0, 5, 3), 0),
                ((0, 5, 4), 0),
                ((0, 6, 3), 0),
                ((0, 6, 4), 0),
                ((0, 7, 3), 0),
                ((0, 7, 4), 0)
            ]
        else:
            spd = [
                ((0, 5, 3), -100),
                ((0, 5, 4), -100),
                ((0, 6, 3), -100),
                ((0, 6, 4), -100),
                ((0, 7, 3), -100),
                ((0, 7, 4), -100)
            ]
            period_data[kper] = spd

    wel = flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=period_data,
        mover=True
    )

    # create uzf package
    nuzfcells = 116
    ntrailwaves = 200
    nwavesets = 500
    package_data = []
    cnt = 0
    cid = 0
    id = idomain.ravel()
    finf = np.loadtxt(os.path.join(data_pth, "finf.txt"), dtype=float) * 1e-10
    finf = finf.ravel()
    pet = np.loadtxt(os.path.join(data_pth, "pet.txt"), dtype=float)
    for i in range(nrow):
        for j in range(ncol):
            vks = 4.6e-05
            if id[cid] == 0:
                cid += 1
                continue

            rec = (cnt, (0, i, j), 1, 0, 1.0, vks, 0.2, 0.38, 0.2, 7.5)
            package_data.append(rec)
            cid += 1
            cnt += 1

    period_data = {}
    for per in range(nper):
        spd = []
        cnt = 0
        for i in range(nrow * ncol):
            if id[i] == 0:
                continue
            rec = (cnt, finf[i], pet[per], 0.5, 0.2, -1.1, -75, 1.0)
            spd.append(rec)
            cnt += 1

        period_data[per] = spd

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

    # SFR Package
    strk = 1.5e-05
    sfr_pakdata = [
        [0, (0, 0, 0), 4500.0, 12, 8.6767896e-04, 1093.048, 3.0, strk, 0.030, 1, 1.0, 0,],
        [1, (0, 1, 1), 7000.0, 12, 8.6767896e-04, 1088.059, 3.0, strk, 0.030, 2, 1.0, 0,],
        [2, (0, 2, 2), 6000.0, 12, 8.6767896e-04, 1082.419, 3.0, strk, 0.030, 2, 1.0, 0,],
        [3, (0, 2, 3), 5550.0, 12, 8.6767896e-04, 1077.408, 3.0, strk, 0.030, 3, 1.0, 1,],
        [4, (0, 3, 4), 6500.0, 12, 9.4339624e-04, 1071.934, 3.0, strk, 0.030, 2, 1.0, 0,],
        [5, (0, 4, 5), 5000.0, 12, 9.4339624e-04, 1066.509, 3.0, strk, 0.030, 2, 1.0, 0,],
        [6, (0, 5, 5), 5000.0, 12, 9.4339624e-04, 1061.792, 3.0, strk, 0.030, 2, 1.0, 0,],
        [7, (0, 6, 5), 5000.0, 12, 9.4339624e-04, 1057.075, 3.0, strk, 0.030, 2, 1.0, 0,],
        [8, (0, 7, 5), 5000.0, 12, 9.4339624e-04, 1052.359, 3.0, strk, 0.030, 2, 1.0, 0,],
        [9, (0, 2, 4), 5000.0, 10, 5.4545456e-04, 1073.636, 2.0, strk, 0.030, 2, 0.0, 0,],
        [10, (0, 2, 5), 5000.0, 10, 5.4545456e-04, 1070.909, 2.0, strk, 0.030, 2, 1.0, 0,],
        [11, (0, 2, 6), 4500.0, 10, 5.4545456e-04, 1068.318, 2.0, strk, 0.030, 2, 1.0, 0,],
        [12, (0, 3, 7), 6000.0, 10, 5.4545456e-04, 1065.455, 2.0, strk, 0.030, 2, 1.0, 0,],
        [13, (0, 4, 7), 5000.0, 10, 5.4545456e-04, 1062.455, 2.0, strk, 0.030, 2, 1.0, 0,],
        [14, (0, 5, 7), 2000.0, 10, 5.4545456e-04, 1060.545, 2.0, strk, 0.030, 2, 1.0, 0,],
        [15, (0, 4, 9), 2500.0, 10, 1.8181818e-03, 1077.727, 3.0, strk, 0.030, 1, 1.0, 0,],
        [16, (0, 4, 8), 5000.0, 10, 1.8181818e-03, 1070.909, 3.0, strk, 0.030, 2, 1.0, 0,],
        [17, (0, 5, 7), 3500.0, 10, 1.8181818e-03, 1063.182, 3.0, strk, 0.030, 2, 1.0, 0,],
        [18, (0, 5, 7), 4000.0, 15, 1.0000000e-03, 1058.000, 3.0, strk, 0.030, 3, 1.0, 0,],
        [19, (0, 6, 6), 5000.0, 15, 1.0000000e-03, 1053.500, 3.0, strk, 0.030, 2, 1.0, 0,],
        [20, (0, 7, 6), 3500.0, 15, 1.0000000e-03, 1049.250, 3.0, strk, 0.030, 2, 1.0, 0,],
        [21, (0, 7, 5), 2500.0, 15, 1.0000000e-03, 1046.250, 3.0, strk, 0.030, 2, 1.0, 0,],
        [22, (0, 8, 5), 5000.0, 12, 9.0909092e-04, 1042.727, 3.0, strk, 0.030, 3, 1.0, 0,],
        [23, (0, 9, 6), 5000.0, 12, 9.0909092e-04, 1038.182, 3.0, strk, 0.030, 2, 1.0, 0,],
        [24, (0, 10, 6), 5000.0, 12, 9.0909092e-04, 1033.636, 3.0, strk, 0.030, 2, 1.0, 0,],
        [25, (0, 11, 6), 5000.0, 12, 9.0909092e-04, 1029.091, 3.0, strk, 0.030, 2, 1.0, 0,],
        [26, (0, 12, 6), 2000.0, 12, 9.0909092e-04, 1025.909, 3.0, strk, 0.030, 2, 1.0, 0,],
        [27, (0, 13, 8), 5000.0, 55, 9.6774194e-04, 1037.581, 3.0, 0.00006, 0.025, 1, 1.0, 0,],
        [28, (0, 12, 7), 5500.0, 55, 9.6774194e-04, 1032.500, 3.0, 0.00006, 0.025, 2, 1.0, 0,],
        [29, (0, 12, 6), 5000.0, 55, 9.6774194e-04, 1027.419, 3.0, 0.00006, 0.025, 2, 1.0, 0,],
        [30, (0, 12, 5), 5000.0, 40, 1.2500000e-03, 1021.875, 3.0, 0.00006, 0.025, 3, 1.0, 0,],
        [31, (0, 12, 4), 5000.0, 40, 1.2500000e-03, 1015.625, 3.0, 0.00006, 0.025, 2, 1.0, 0,],
        [32, (0, 12, 3), 5000.0, 40, 1.2500000e-03, 1009.375, 3.0, 0.00006, 0.025, 2, 1.0, 0,],
        [33, (0, 12, 2), 5000.0, 40, 1.2500000e-03, 1003.125, 3.0, 0.00006, 0.025, 2, 1.0, 0,],
        [34, (0, 12, 1), 5000.0, 40, 1.2500000e-03, 996.8750, 3.0, 0.00006, 0.025, 2, 1.0, 0,],
        [35, (0, 12, 0), 3000.0, 40, 1.2500000e-03, 991.8750, 3.0, 0.00006, 0.025, 1, 1.0, 0,],
    ]

    sfr_conn = [
        [0, -1],
        [1, 0, -2],
        [2, 1, -3],
        [3, 2, -4, -9],
        [4, 3, -5],
        [5, 4, -6],
        [6, 5, -7],
        [7, 6, -8],
        [8, 7, -22],
        [9, 3, -10],
        [10, 9, -11],
        [11, 10, -12],
        [12, 11, -13],
        [13, 12, -14],
        [14, 13, -18],
        [15, -16],
        [16, 15, -17],
        [17, 16, -18],
        [18, 14, 17, -19],
        [19, 18, -20],
        [20, 19, -21],
        [21, 20, -22],
        [22, 8, 21, -23],
        [23, 22, -24],
        [24, 23, -25],
        [25, 24, -26],
        [26, 25, -30],
        [27, -28],
        [28, 27, -29],
        [29, 28, -30],
        [30, 26, 29, -31],
        [31, 30, -32],
        [32, 31, -33],
        [33, 32, -34],
        [34, 33, -35],
        [35, 34],
    ]

    sfr_div = [[3, 0, 9, "UPTO"]]

    period_data = {}
    for per in range(nper):
        sfr_spd = [
            [0, "inflow", 50.0],
            [15, "inflow", 10.0],
            [27, "inflow", 150.0],
            [3, "diversion", 0, 10.0],
            # [9, "status", "simple"],
            # [10, "status", "simple"],
            # [11, "status", "simple"],
            # [12, "status", "simple"],
            # [13, "status", "simple"],
            # [14, "status", "simple"],
            # [9, "stage", 1075.545],
            # [10, "stage", 1072.636],
            # [11, "stage", 1069.873],
            # [12, "stage", 1066.819],
            # [13, "stage", 1063.619],
            # [14, "stage", 1061.581],
        ]
        period_data[per] = sfr_spd

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=len(sfr_pakdata),
        packagedata=sfr_pakdata,
        connectiondata=sfr_conn,
        diversions=sfr_div,
        perioddata=period_data,
        unit_conversion=None,
        save_flows=True,
        mover=True
    )

    # build output control package
    budget_file = f"{name}.cbc"
    head_file = f"{name}.hds"
    saverecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(nper)}
    printrecord = {i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(nper)}
    oc = flopy.mf6.ModflowGwfoc(gwf,
                                budget_filerecord=budget_file,
                                head_filerecord=head_file,
                                saverecord=saverecord,
                                printrecord=printrecord)

    # build AG MVR package
    maxpackages = 3
    packages=[("sfr_0",), ("uzf_0",), ("wel_0")]

    period_data = {}
    # total SFR diversion amt is 100
    div_p_cell = 100 / 6
    ag_pump = 100
    sfr_reach = 8
    uzfcells = (37, 38, 47, 48, 55, 56)
    for per in range(nper):
        spd = []
        if per <= 3:
            period_data[per] = []
        else:
            for cell in uzfcells:
                rec = ("sfr_0", sfr_reach, "uzf_0", cell, "UPTO", div_p_cell)
                spd.append(rec)

            for ix, cell in enumerate(uzfcells):
                rec = ("wel_0", ix, "uzf_0", cell, "UPTO", ag_pump)
                spd.append(rec)

            period_data[per] = spd

    mvr = flopy.mf6.ModflowGwfmvr(
        gwf,
        maxmvr=len(period_data[4]),
        maxpackages=maxpackages,
        packages=packages,
        perioddata=period_data,
        budgetcsv_filerecord=f"{name}_mvr.csv",
        budget_filerecord=f"{name}_mvr.cbc"
    )

    sim.write_simulation()
    return sim, gwf


if __name__ == "__main__":
    sim_ws = os.path.join("..", "..", "data", "mf6_prudic_ag")
    name = "prudic_ag"
    sim, gwf = build_model(name, sim_ws)

    mfag = ModflowApiAg(sim, ag_type="etdemand", mvr_name="mvr")
    mfag.run_model(develop=False)

