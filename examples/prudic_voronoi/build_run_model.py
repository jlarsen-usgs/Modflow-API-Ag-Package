import flopy
from flopy.utils.triangle import Triangle
from flopy.utils.voronoi import VoronoiGrid
from flopy.discretization import StructuredGrid, VertexGrid
from flopy.plot.plotutil import UnstructuredPlotUtilities
from mf6api_ag import ModflowApiAg


import shapefile
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from scipy.interpolate import griddata


def shoelace_area(modelgrid):
    """
    Use shoelace algorithm for non-self-intersecting polygons to
    calculate area.
    """
    xverts, yverts = modelgrid.cross_section_vertices
    xverts, yverts = UnstructuredPlotUtilities.irregular_shape_patch(
        xverts, yverts
    )
    area_x2 = np.zeros((1, len(xverts)))
    for i in range(xverts.shape[-1]):
        # calculate the determinant of each line in polygon
        area_x2 += xverts[:, i - 1] * yverts[:, i] - \
                   yverts[:, i - 1] * xverts[:, i]

    area = area_x2 / 2.
    return np.ravel(area)


def resample_structured_to_unstructured(sgrid, vgrid, array, method="nearest"):
    array = array.ravel()
    sxc = sgrid.xcellcenters.ravel()
    syc = sgrid.ycellcenters.ravel()
    sid = sgrid.idomain.ravel()
    sxc = sxc[sid > 0]
    syc = syc[sid > 0]

    vxc = vgrid.xcellcenters.ravel()
    vyc = vgrid.ycellcenters.ravel()

    array = array[sid > 0]

    new_array = griddata((sxc, syc), array, (vxc, vyc), method=method)

    if method != "nearest":
        buffer_array = griddata(
            (sxc, syc), array, (vxc, vyc), method="nearest"
        )
        new_array = np.where(np.isnan(new_array), buffer_array, new_array)

    return new_array

def build_simulation():
    ws = os.path.abspath(os.path.dirname(__file__))
    sim_ws = os.path.join(ws, "..", "..", "data", "mf6_prudic_voronoi")
    data_ws = os.path.join(ws, "data")

    boundary = os.path.join(data_ws, "green_valley_boundary.shp")
    lake = os.path.join(data_ws, "green_valley_lake.shp")
    ag = os.path.join(data_ws, "green_valley_ag.shp")

    with shapefile.Reader(boundary) as r:
        boundary = r.shape(0)

    with shapefile.Reader(lake) as r:
        lake = r.shape(0)

    with shapefile.Reader(ag) as r:
        ag = [shape for shape in r.shapes()]

    # Build a triangular boundary
    tri = Triangle(
        maximum_area=5000 * 2500, angle=30, model_ws=sim_ws
    )
    tri.add_polygon(boundary)
    tri.build()

    # create a voronoi grid from the tri object
    vor = VoronoiGrid(tri)
    gridprops = vor.get_gridprops_vertexgrid()

    # create temporary top and botm and then adjust later
    top = np.ones(vor.ncpl)
    botm = np.zeros((2, vor.ncpl))
    botm[-1] -= 1

    idomain = np.ones((2, vor.ncpl), dtype=int)

    vor_grid = VertexGrid(
        idomain=idomain,
        nlay=2,
        top=top,
        botm=botm,
        **gridprops
    )

    # Data for resampling from structured to vertex
    botm = os.path.join(data_ws, "bottom.txt")
    top = os.path.join(data_ws, "top.txt")
    strt = os.path.join(data_ws, "strt.txt")
    sy = os.path.join(data_ws, "sy.txt")
    hk = os.path.join(data_ws, "hk.txt")
    # uzf data
    surf = os.path.join(data_ws, "surf.txt")
    finf = os.path.join(data_ws, "finf.txt")
    pet = os.path.join(data_ws, "pet.txt")
    # sfr data

    # create a structured grid for resampling
    idomain = np.genfromtxt(os.path.join(data_ws, "idomain.txt"), dtype=int)
    nlay, nrow, ncol = 1, 15, 10
    stop = np.ones((nrow, ncol))
    sbotm = np.ones((nlay, nrow, ncol))
    delr = np.ones((ncol,)) * 5000
    delc = np.ones((nrow,)) * 5000
    struct_grid = StructuredGrid(
        delc=delc, delr=delr, top=stop, botm=sbotm, nlay=nlay, idomain=idomain
    )

    # Create voronoi model simulation
    sim = flopy.mf6.MFSimulation(sim_ws=sim_ws)

    nper = 13
    perlen = [(1, 1, 1)]
    perlen += [(30.5, 30, 1) for _ in range(nper - 1)]
    tdis = flopy.mf6.ModflowTdis(
        sim,
        nper=nper,
        perioddata=tuple(perlen),
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

    # create voronoi model
    gwf = flopy.mf6.ModflowGwf(
        sim,
        modelname="prudic_vor",
        save_flows=True,
        print_input=True,
        print_flows=True,
        newtonoptions="NEWTON UNDER RELAXATION"
    )

    # resample top and botm to voronoi grid and create DISV package
    top = np.genfromtxt(top)
    vor_top = resample_structured_to_unstructured(
        struct_grid, vor_grid, top, method="linear"
    )

    # consider doing mean([min(vor_top), max(vor_botm)])
    botm = np.genfromtxt(botm)
    vor_model_botm = resample_structured_to_unstructured(
        struct_grid, vor_grid, botm, method="linear"
    )

    vor_botm = np.zeros((2, vor_grid.ncpl))
    vor_botm[-1] = vor_model_botm
    vor_botm[0] = np.mean([vor_top, vor_model_botm], axis=0)

    disv = flopy.mf6.ModflowGwfdisv(
        gwf,
        nlay=vor_grid.nlay,
        top=vor_top,
        botm=vor_botm,
        idomain=vor_grid.idomain,
        **gridprops
    )

    # get the final modelgrid object
    vor_grid = gwf.modelgrid

    # create IC package
    strt = np.genfromtxt(strt)
    vor_strt = resample_structured_to_unstructured(
        struct_grid, vor_grid, strt, method="linear"
    )

    ic = flopy.mf6.ModflowGwfic(
        gwf,
        strt=[vor_strt, vor_strt]
    )

    # create NPF package
    hk = np.genfromtxt(hk)
    vor_hk = resample_structured_to_unstructured(
        struct_grid, vor_grid, hk, method="linear"
    )

    npf = flopy.mf6.ModflowGwfnpf(
        gwf,
        icelltype=1,
        k=[vor_hk, vor_hk],
        k33=1e-06,
    )

    # create STO package: this will take some math to calculate an equivalent SY
    sy = np.genfromtxt(sy)
    vor_sy = resample_structured_to_unstructured(
        struct_grid, vor_grid, sy, method="linear"
    )

    ss = vor_sy / (vor_top - vor_botm[0]) # should this be updated???

    steady_state = {0: True}
    transient = {0: False}
    for i in range(1, nper):
        steady_state[i] = False
        transient[i] = True


    sto = flopy.mf6.ModflowGwfsto(
        gwf,
        iconvert=1,
        ss=1e-06,
        sy=[vor_sy, vor_sy],
        steady_state=steady_state,
        transient=transient
    )

    # create GHB package
    struct_rec = [
        (0, 12, 0, 988.0, 0.038),
        (0, 13, 8, 1045.0, 0.038)
    ]
    vor_rec = []
    for _, r, c, elev, cond in struct_rec:
        xc = struct_grid.xcellcenters[r, c]
        yc = struct_grid.ycellcenters[r, c]
        node = vor_grid.intersect(xc, yc)
        for lay in range(2):
            vor_rec.append(((lay, node), elev, cond))

    ghb = flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data={i: vor_rec for i in range(nper)}
    )

    # create the UZF package
    finf = np.genfromtxt(finf) * 1e-10
    pet = np.genfromtxt(pet)

    vor_finf = resample_structured_to_unstructured(
        struct_grid, vor_grid, finf, method="nearest"
    )

    packagedata = []
    for node in range(vor_grid.ncpl):
        rec = (node, (0, node), 1, 0, 1.0, 4.6e-05, 0.2, 0.38, 0.2, 7.5)
        packagedata.append(rec)

    perioddata = {}
    for per in range(nper):
        spd = []
        for node in range(vor_grid.ncpl):
            rec = (node, vor_finf[node], pet[per], 0.5, 0.2, -1.1, -75, 1.0)
            spd.append(rec)

        perioddata[per] = spd

    uzf = flopy.mf6.ModflowGwfuzf(
        gwf,
        simulate_et=True,
        nuzfcells=vor_grid.ncpl,
        ntrailwaves=15,
        nwavesets=100,
        packagedata=packagedata,
        perioddata=perioddata,
        unsat_etwc=True,
        linear_gwet=True,
        simulate_gwseep=True,
        mover=True
    )

    # Build a LAK package
    gix = flopy.utils.GridIntersect(vor_grid, method="vertex")
    result = gix.intersect(lake, contains_centroid=True)
    lak_cellids = result.cellids.ravel()

    lak_strt = np.ceil(np.mean([vor_grid.top[cid] for cid in lak_cellids]))
    packagedata = [(0, lak_strt, len(lak_cellids)),]

    connectiondata = []
    for ix, cid in enumerate(lak_cellids):
        rec = (
            0,
            ix,
            (0, cid),
            "vertical",
            1e-06,
            vor_grid.top[cid] - 1,
            np.where(lak_strt > vor_grid.top[cid], lak_strt, vor_grid.top[cid]),
            1,
            1
        )
        connectiondata.append(rec)

    outlets = [(0, 0, -1, "specified", 0, 0, 0, 0),]

    spd = [(0, "rate", -50), (0, "status", "constant"), (0, "stage", lak_strt)]
    perioddata = {i: spd for i in range(nper)}

    lak = flopy.mf6.ModflowGwflak(
        gwf,
        mover=True,
        nlakes=len(packagedata),
        noutlets=len(outlets),
        ntables=0,
        packagedata=packagedata,
        connectiondata=connectiondata,
        outlets=outlets,
        perioddata=perioddata
    )


    # build MAW package
    maw_cells = []
    gix = flopy.utils.GridIntersect(vor_grid, method="vertex")
    for shp in ag:
        result = gix.intersect(shp, contains_centroid=True)
        cellids = [i for i in result.cellids.ravel()]
        xc = np.mean(vor_grid.xcellcenters[cellids])
        yc = np.mean(vor_grid.ycellcenters[cellids])
        cid = vor_grid.intersect(xc, yc)
        maw_cells.append(cid)

    # todo: perioddata
    packagedata = []
    for ix, cid in enumerate(maw_cells):
        rec = (ix, 0.5, vor_grid.botm[1, cid], vor_strt[cid], "SKIN", 2)
        packagedata.append(rec)

    connectiondata = []
    cnt = 0
    for cid in maw_cells:
        for lay in range(2):
            rec = (
                cnt, lay, (lay, cid), vor_strt[cid], vor_grid.botm[lay, cid], 0.001, 1
            )
            connectiondata.append(rec)
        cnt += 1

    perioddata = {}
    for i in range(nper):
        spd = []
        rate = -100
        if i < 4:
            rate = 0
        for ix in range(len(packagedata)):
            spd.append([ix, "rate", rate])
        perioddata[i] = spd

    maw = flopy.mf6.ModflowGwfmaw(
        gwf,
        mover=True,
        nmawwells=len(packagedata),
        packagedata=packagedata,
        connectiondata=connectiondata,
        perioddata=perioddata
    )

    # build the SFR package
    df = pd.read_csv(os.path.join('data', 'prudic_voronoi_sfr.txt'))
    # node 49, 5, 1,
    print('break')

    df.node -= 1

    arr = np.zeros(vor_grid.top.shape)
    arr[[int(v) for v in df.node.values]] = df.sfr_seg.values
    pmv = flopy.plot.PlotMapView(gwf)
    pmv.plot_array(arr, masked_values=[0])
    pmv.plot_bc("LAK")
    pmv.plot_grid()
    plt.show()

    nreaches = len(df)

    # build packagedata
    packagedata = []
    for iloc, row in df.iterrows():
        if row.mf6_rchno in (0, 26, 35, 49):
            ncon = 1
        elif row.mf6_rchno in (11, 18, 31):
            ncon = 3
        else:
            ncon = 2
        div = 0
        if row["div"] != 0:
            div = 1
            ncon += 1

        ustrf = 1
        if row.mf6_rchno == 42:
            ustrf = 0

        rec = [int(row.mf6_rchno), (0, int(row.node)), row.rlen, row.rwid, row.rgrd, row.strtop,
               row.rbth, row.strk, row.man, ncon, ustrf, div]
        packagedata.append(rec)

    # Build connection data
    # - is downstream, + is upstream
    connectiondata = []
    for iloc, row in df.iterrows():
        if row.mf6_rchno in (0, 26, 35, 49,):
            if row.mf6_rchno in (0, 35):
                rec = [int(row.mf6_rchno), -1 * (int(row.mf6_rchno) + 1)]
            elif row.mf6_rchno == 26:
                rec = [int(row.mf6_rchno), int(row.mf6_rchno) - 1]
            else:
                rec = [int(row.mf6_rchno), -27]

        elif row.mf6_rchno == 3:
            rec = [3, 2, -4, -42]
        elif row.mf6_rchno == 18:
            rec = [18, 17, 41, -19]
        elif row.mf6_rchno == 27:
            rec = [27, 49, -28]
        elif row.mf6_rchno == 31:
            rec = [31, 30, 48, -32]
        elif row.mf6_rchno == 11:
            rec = [11, 10, 34, -12]
        elif row.mf6_rchno == 42:
            rec = [42, 3, -43]
        else:
            dnrch = int(row.mf6_rchno + 1)
            if not np.isnan(row.rchto):
                dnrch = int(row.rchto)
            rec = [int(row.mf6_rchno), int(row.mf6_rchno) - 1, -1 * dnrch]

        connectiondata.append(rec)

    sfr_div = [[3, 0, 42, "UPTO"],]

    perioddata = {}
    for per in range(nper):
        sfr_spd = [
            [0, "inflow", 50],
            [49, "inflow", 10],
            [35, "inflow", 150],
            [3, "diversion", 0, 10.0]
        ]

    sfr = flopy.mf6.ModflowGwfsfr(
        gwf,
        nreaches=nreaches,
        packagedata=packagedata,
        connectiondata=connectiondata,
        diversions=sfr_div,
        perioddata=perioddata,
        unit_conversion=None,
        save_flows=True,
        mover=True
    )

    # build a simple MVR package to simulate AG
    # todo: add LAK SFR record to MVR @ 50 or 100
    ag_cells = []
    ag_cells_flat = []
    gix = flopy.utils.GridIntersect(vor_grid, method="vertex")
    ncells = 0
    for shp in ag:
        result = gix.intersect(shp, contains_centroid=True)
        ag_cells_flat += list(result.cellids.ravel())
        ag_cells.append(result.cellids.ravel())
        ncells += result.cellids.size

    arr = np.zeros(vor_grid.top.shape)
    arr[[int(v) for v in df.node.values]] = df.sfr_seg.values
    arr2 = np.zeros(vor_grid.top.shape)
    arr2[ag_cells_flat] = 1
    pmv = flopy.plot.PlotMapView(gwf)
    pmv.plot_array(arr, masked_values=[0])
    pmv.plot_array(arr2, masked_values=[0], cmap="summer")
    # pmv.plot_bc("SFR")
    pmv.plot_bc("LAK")
    pmv.plot_grid()
    plt.show()


    agpump = 100
    div_per_cell = 50 / ncells
    app_eff = 0.872
    packages = [
        (lak.package_name,),
        (sfr.package_name,),
        (uzf.package_name,),
        (maw.package_name,)
    ]

    spd = [(lak.package_name, 0, sfr.package_name, 0, 50, 1, 1)]
    for cid in ag_cells_flat:
        rec = (lak.package_name, 0, uzf.package_name, cid, div_per_cell, 1, app_eff) # / 2)
        spd.append(rec)

    # for cid in ag_cells_flat:
    #     rec = (sfr.package_name, 9, uzf.package_name, cid, "UPTO", div_per_cell / 2)
    #     spd.append(rec)

    for welno, ag_cids in enumerate(ag_cells):
        pump_per_cell = agpump / len(ag_cids)
        for cid in list(ag_cids):
            rec = (
                maw.package_name,
                welno,
                uzf.package_name,
                cid,
                pump_per_cell,
                1,
                app_eff
            )
            spd.append(rec)

    perioddata = {}
    for per in range(nper):
        if per < 4:
            perioddata[per] = []
        else:
            perioddata[per] = spd

    mvr = flopy.mf6.ModflowGwfapiag(
        gwf,
        maxmvr=len(perioddata[4]),
        maxpackages=len(packages),
        packages=packages,
        perioddata=perioddata,
        budgetcsv_filerecord="gv_mvr_bud.csv",
        pname='agmvr'
    )

    sim.write_simulation()
    # sim.run_simulation()
    return sim


sim = build_simulation()

mfag = ModflowApiAg(sim, ag_type="etdemand", mvr_name="apiag")
mfag.run_model(develop=False)
