import flopy as fp
import flopy.mf6
import matplotlib.colors
from flopy.plot import styles
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from matplotlib.collections import PathCollection
from matplotlib.path import Path
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import matplotlib as mpl
mpl.use("TKAgg")


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


def get_outputs():
    # now load up outputs from models
    data_ws = os.path.join("..", "..", "data")
    nwt_ws = os.path.join(data_ws, "nwt_prudic_ag")
    mf6_ws = os.path.join(data_ws, "mf6_prudic_voronoi")

    nwt_file = os.path.join(nwt_ws, "prudic_ag.diversion9.txt")
    mf6_file = os.path.join(mf6_ws, "prudic_vor_ag.out")

    nwt_ag = pd.read_csv(nwt_file, delim_whitespace=True)
    mf6_ag = pd.read_csv(mf6_file, delim_whitespace=True)

    nwt_cbc = os.path.join(nwt_ws, "prudic_ag.cbc")
    nwt_cbc = flopy.utils.CellBudgetFile(nwt_cbc)

    conv = 0.0283168  # ft3 to m3
    well_recs = nwt_cbc.get_data(text="AG WE")

    mf6_ag = mf6_ag.groupby(by=["pid", "pkg", "kstp"], as_index=False).sum()

    mf6_ag_lak = mf6_ag[mf6_ag.pkg == "lak_0"]
    mf6_ag_lak.reset_index(drop=True, inplace=True)

    total_mf6_div = mf6_ag_lak.q_from_provider.values

    mf6_ag_pump = mf6_ag[mf6_ag.pkg == "maw_0"]
    mf6_ag_pump = mf6_ag_pump.groupby(by=["pkg", "kstp"], as_index=False)[
        "q_from_provider"].sum()

    nwt_total = (nwt_ag["SW-DIVERSION"].values + nwt_ag[
        "SUP-PUMPING"].values) * conv
    mf6_total = (total_mf6_div + mf6_ag_pump.q_from_provider.values) * conv

    err = np.sum(nwt_total[:361] - np.array([0] * 91 + list(mf6_total)))

    print("SUM of NWT pumpage:", np.sum(nwt_total))

    print("volume error: ", err)
    print("percent error: ", (err / np.sum(nwt_total)) * 100)
    return mf6_ag_lak, mf6_total, nwt_ag, nwt_total


def plot_mvr(ax, gwf, provider, kper, arrows=False, head_width=10, **kwargs):
    modelgrid = gwf.modelgrid
    idomain = modelgrid.idomain[0].ravel()
    xc = modelgrid.xcellcenters.ravel()
    yc = modelgrid.ycellcenters.ravel()
    parray = np.zeros(xc.shape)
    rarray = np.zeros(xc.shape)
    uzfcell = np.ones(xc.shape) * np.nan
    cnt = 0
    for ix, ib in enumerate(idomain):
        if ib == 0:
            continue
        uzfcell[ix] = cnt
        cnt += 1

    recarray = gwf.mvr.perioddata.data[kper]
    idx = np.where((recarray["pname1"] == provider) & (recarray["pname2"] != "sfr_0") )[0]
    recarray = recarray[idx]
    nodes = []

    if "lak" in provider:
        providerdata = gwf.lak.connectiondata.array
        nodes = [330,] * len(recarray)
    elif "maw" in provider:
        tmp = gwf.maw
        providerdata = gwf.maw.connectiondata.array
        for rec in recarray:
            mawno = rec["id1"]
            ra = providerdata[providerdata.wellno == mawno]
            nodes.append(ra.cellid[0][-1])
    else:
        return

    for rid in recarray["id1"]:
        if "lak" in provider:
            pass
        else:
            rec = providerdata[[rid]]
        if not nodes:
            cellids = rec['cellid']
            if len(cellids[0]) == 3:
                node_list = modelgrid.get_node(list(cellids))
                for node in node_list:
                    while node > modelgrid.ncpl:
                        node -= modelgrid.ncpl
                    nodes.append(node)
            else:
                for node in cellids:
                    while node[0] < modelgrid.ncpl:
                        node -= modelgrid.ncpl
                    nodes.append(node)

    # todo: need 1:1 association of nodes to recievers

    parray[nodes] = 1
    pxy = [[xc[i], yc[i]] for i in nodes]

    # todo: need to link reciever to uzf node number
    adj_id = []
    for rid in recarray["id2"]:
        tid = np.where(uzfcell == rid)[0][0]
        adj_id.append(tid)

    rarray[adj_id] = 1
    rxy = [[xc[i], yc[i]] for i in adj_id]
    parray = np.ma.masked_equal(parray, 0)
    rarray = np.ma.masked_equal(rarray, 0)
    # parray[parray != 1] = np.nan
    # rarray[rarray != 1] = np.nan

    # use cached patch collection for plotting
    if not arrows:
        polygons = modelgrid.map_polygons
        if isinstance(polygons, dict):
            polygons = polygons[0]

        if len(polygons) == 0:
            return

        if not isinstance(polygons[0], Path):
            collection0 = ax.pcolormesh(
                modelgrid.xvertices, modelgrid.yvertices, parray.reshape(modelgrid.nrow, modelgrid.ncol)
            )
            collection1 = ax.pcolormesh(
                modelgrid.xvertices, modelgrid.yvertices, rarray.reshape(modelgrid.nrow, modelgrid.ncol)
            )

        else:
            parray = parray.ravel()
            collection0 = PathCollection(polygons)
            collection0.set_array(parray)

            rarray = rarray.ravel()
            collection1 = PathCollection(polygons)
            collection1.set_array(rarray)

        if "lak" in provider:
            cmap = matplotlib.colors.ListedColormap(["0", "c"])
            bounds = [0, 1, 2]
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        else:
            cmap = matplotlib.colors.ListedColormap(["0", "darkred"])
            bounds = [0, 1, 2]
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

        collection0.set(cmap=cmap, norm=norm)

        cmap = matplotlib.colors.ListedColormap(["0", "darkgreen"])
        bounds = [0, 1, 2]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

        collection1.set(cmap=cmap, norm=norm)
        ax.add_collection(collection0)
        ax.add_collection(collection1)
        return collection0, collection1
    else:
        if "lak" in provider:
            color = "darkblue"
        else:
            color = "indianred"
        arrows = []
        for ix, (x, y) in enumerate(pxy):
            rx, ry = rxy[ix]
            dx = rx - x
            dy = ry - y

            if dx < 1 and dy < 1:
                continue

            if dx < 0:
                dx += (head_width * 1.1)
            elif dx > 0:
                dx -= (head_width * 1.1)
            else:
                pass

            if dy < 0:
                dy += (head_width * 1.1)
            elif dy > 0:
                dy -= (head_width * 1.1)
            else:
                pass

            arrw = ax.arrow(x, y, dx, dy, color=color, head_width=head_width, zorder=3)
            arrows.append(arrw)

        return arrows


def nwt_plot_diversion(ax, ml, provider, kper, arrows=False, head_width=10, **kwargs):
    modelgrid = ml.modelgrid
    xc = modelgrid.xcellcenters
    yc = modelgrid.ycellcenters
    parray = np.zeros((modelgrid.nrow, modelgrid.ncol))
    rarray = np.zeros((modelgrid.nrow, modelgrid.ncol))

    if provider == "sfr":
        recarray = ml.ag.irrdiversion[kper]
        pdata = ml.sfr.reach_data
        rdata = None
    elif provider == "well":
        recarray = ml.ag.irrwell[kper]
        pdata = ml.ag.well_list
        rdata = None
    else:
        recarray = ml.ag.supwell[kper]
        pdata = ml.ag.well_list
        rdata = ml.ag.irrdiversion[kper]

    pnode = []
    if provider == "sfr":
        segids = recarray["segid"]
        for ix, segid in enumerate(segids):
            idx = np.where(pdata["iseg"] == segid)
            rec = pdata[idx]
            parray[rec["i"], rec["j"]] = 1
            for i in range(recarray[ix]["numcell"]):
                pnode.append([rec["i"][0], rec["j"][0]])
    else:
        wellids = recarray["wellid"]
        for ix, wellid in enumerate(wellids):
            rec = pdata[wellid]
            parray[rec["i"], rec["j"]] = 1
            if provider == "sup":
                segids = [recarray[ix][f"segid{i}"] for i in range(recarray[ix]["numcell"])]
                for segid in segids:
                    idx = np.where(rdata["segid"] == segid)[0]
                    n = rdata[idx]["numcell"]
                    for i in range(n[0]):
                        pnode.append([rec["i"], rec["j"]])
            else:
                for i in range(recarray[ix]["numcell"]):
                    pnode.append([rec["i"], rec["j"]])

    rnode = []
    if provider == "sfr" or provider == "well":
        for record in recarray:
            n = record["numcell"]
            for i in range(n):
                rarray[record[f"i{i}"], record[f"j{i}"]] = 1
                rnode.append([record[f"i{i}"], record[f"j{i}"]])
    else:
        for record in recarray:
            n = record["numcell"]
            segids = []
            for i in range(n):
                segids.append(record[f"segid{i}"])

            for segid in segids:
                idx = np.where(rdata["segid"] == segid)
                rec = rdata[idx]
                n = rec["numcell"][0]
                for i in range(n):
                    rarray[rec[f"i{i}"], rec[f"j{i}"]] = 1
                    rnode.append([rec[f"i{i}"][0], rec[f"j{i}"][0]])

    rarray[rarray != 1] = np.nan
    parray[parray != 1] = np.nan

    pxy = [[xc[i, j], yc[i, j]] for i, j in pnode]
    rxy = [[xc[i, j], yc[i, j]] for i, j in rnode]

    if not arrows:
        polygons = modelgrid.map_polygons
        if isinstance(polygons, dict):
            polygons = polygons[0]

        if len(polygons) == 0:
            return

        if not isinstance(polygons[0], Path):
            collection0 = ax.pcolormesh(
                modelgrid.xvertices, modelgrid.yvertices,
                parray.reshape(modelgrid.nrow, modelgrid.ncol)
            )
            collection1 = ax.pcolormesh(
                modelgrid.xvertices, modelgrid.yvertices,
                rarray.reshape(modelgrid.nrow, modelgrid.ncol)
            )

        else:
            parray = parray.ravel()
            collection0 = PathCollection(polygons)
            collection0.set_array(parray)

            rarray = rarray.ravel()
            collection1 = PathCollection(polygons)
            collection1.set_array(rarray)

        if "sfr" in provider:
            cmap = matplotlib.colors.ListedColormap(["0", "c"])
            bounds = [0, 1, 2]
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        else:
            cmap = matplotlib.colors.ListedColormap(["0", "darkred"])
            bounds = [0, 1, 2]
            norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

        collection0.set(cmap=cmap, norm=norm)

        cmap = matplotlib.colors.ListedColormap(["0", "darkgreen"])
        bounds = [0, 1, 2]
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

        collection1.set(cmap=cmap, norm=norm)
        ax.add_collection(collection0)
        ax.add_collection(collection1)
        return collection0, collection1
    else:
        if "lak" in provider:
            color = "darkblue"
        else:
            color = "indianred"
        arrows = []
        for ix, (x, y) in enumerate(pxy):
            rx, ry = rxy[ix]
            dx = rx - x
            dy = ry - y

            if dx < 0:
                dx += (head_width * 1.5)
            elif dx > 0:
                dx -= (head_width * 1.5)
            else:
                pass

            if dy < 0:
                dy += (head_width * 1.5)
            elif dy > 0:
                dy -= (head_width * 1.5)
            else:
                pass

            arrw = ax.arrow(x, y, dx, dy, color=color, head_width=head_width, zorder=3)
            arrows.append(arrw)

        return arrows


def plot_explanation(ax, wells=True, sfr=True):
    patches = []
    labels = []
    if wells:
        patches.append(ax.scatter([-1000], [-1000], c='darkred', marker='o', s=75))
        labels.append("AG provider well")
    if sfr:
        patches.append(mpatches.Patch(color="c"))
        labels.append("LAK AG provider")

    patches.append(mpatches.Patch(color="darkgreen"))
    labels.append("AG receiver cells")
    if sfr:
        patches.append(mpatches.Patch(color="teal"))
        labels.append("SFR cells")

    patches.append(mpatches.Patch(color="b", alpha=0.5))
    labels.append("LAK cells")

    if wells:
        patches.append(
            ax.scatter([-1000], [-1000], c="indianred", marker=r'$\leftarrow$', s=80)
        )
        labels.append("Irrigation from well")
    if sfr:
        patches.append(
            ax.scatter([-1000], [-1000], c="darkblue", marker=r'$\leftarrow$', s=100)
        )
        labels.append("Irrigation from lake")

    patches.append(
        Line2D([-1000], [-1000], color="saddlebrown", lw=1.5)
    )
    labels.append("Land surface contour line")
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2,
                     box.width, box.height * 0.8])

    legend = ax.legend(
        patches, labels, fancybox=True, shadow=True, loc="upper center", frameon=True, fontsize=6, ncol=2, bbox_to_anchor=(0.5, -0.10)
    )
    legend = styles.graph_legend_title(legend, fontsize=6)


def plot_well_points(ax, gwf, kper, size=10):
    recarray = gwf.maw.connectiondata.array
    modelgrid = gwf.modelgrid
    xc = modelgrid.xcellcenters
    yc = modelgrid.ycellcenters

    xypts = [[xc[cellid[-1]], yc[cellid[-1]]] for cellid in recarray["cellid"]]
    x, y = np.array(xypts).T

    ax.scatter(x, y, c="darkred", s=size, marker="o")


sim_ws = os.path.join("..", "..", "data", "mf6_prudic_voronoi")
sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws)
gwf = sim.get_model()
print(gwf.modelgrid.ncpl)

with styles.USGSMap():
    mpl.rcParams["ytick.labelsize"] = 6
    mpl.rcParams["xtick.labelsize"] = 6
    # mpl.rcParams["figure.dpi"] = 300
    fig, (ax, ax2) = plt.subplots(ncols=2, figsize=cm2inch(16.5, 13))
    pmv = fp.plot.PlotMapView(model=gwf, ax=ax)
    pmv.plot_inactive(alpha=0.5)
    pmv.plot_bc("lak_0")
    pmv.plot_bc("sfr_0")
    plot_mvr(ax, gwf, "lak_0", 5)  # debug this. Plotting all LAK and UZF cells instead of only ag
    plot_mvr(ax, gwf, "maw_0", 5)
    pmv.plot_grid(lw=0.5)
    plot_mvr(ax, gwf, "lak_0", 5, arrows=True, head_width=1000)
    plot_mvr(ax, gwf, "maw_0", 5, arrows=True, head_width=1000)
    plot_well_points(ax, gwf, 5, size=15)
    top = gwf.dis.top.array
    idomain = gwf.modelgrid.idomain[0]
    top[idomain == 0] = np.nan
    levels = [1010, 1030, 1050, 1070, 1090]
    contour_set = pmv.contour_array(top, colors="saddlebrown", levels=levels, lw=1.5, tri_mask=True)
    plt.clabel(contour_set, fmt="%.0f", colors="saddlebrown", fontsize=7)

    plot_explanation(ax)
    styles.heading(
        ax=ax, letter="A. ", heading="MF6 API Agricultural Water Use Example 3", fontsize=6
    )
    """
    styles.add_text(
        ax,
        text="Green Creek",
        x=0.13,
        y=0.80,
        bold=False,
        rotation=-45,
        fontsize=6
    )
    styles.add_text(
        ax,
        text="Canal",
        x=0.54,
        y=0.82,
        bold=False,
        fontsize=6
    )
    styles.add_text(
        ax,
        text="Little Creek",
        x=0.65,
        y=0.56,
        bold=False,
        rotation=45,
        fontsize=6
    )
    styles.add_text(
        ax,
        text="Blue River",
        x=0.40,
        y=0.16,
        bold=False,
        fontsize=6
    )
    """
    styles.xlabel(ax=ax, label="Model x-coordinate", fontsize=6)
    styles.ylabel(ax=ax, label="Model y-coordinate", fontsize=6)
    plt.xticks()

    # todo: do compare output here as figure
    conv = 0.0283168
    mf6_ag_lak, mf6_total, nwt_ag, nwt_total = get_outputs()
    ax2.plot(
        mf6_ag_lak.kstp.values + 1,
        mf6_total,
        ":",
        zorder=2,
        label="MF6 total irrigation",
        lw=2,
        color="cyan"
    )
    ax2.plot(
        range(1, len(nwt_ag) + 1),
        nwt_total,
        zorder=1,
        label="NWT total irrigation",
        lw=2.5,
        color="dimgray"
    )

    ax2.plot(
        mf6_ag_lak.kstp.values + 1,
        mf6_ag_lak.q_from_provider.values * conv,
        ":",
        zorder=4,
        label="MF6 LAK water irrigation",
        lw=2.5,
        color="yellowgreen"
    )

    ax2.plot(
        range(1, len(nwt_ag) + 1),
        nwt_ag["SW-DIVERSION"].values * conv,
        zorder=3,
        label="NWT surface water irrigation",
        lw=2,
        color="k"
    )

    styles.xlabel(ax=ax2, label="Timestep", fontsize=7)
    styles.ylabel(ax=ax2, label="Applied irrigation, in " + r"$m^{3}/day$",
                  fontsize=7)

    legend = ax2.legend(
        fancybox=True, shadow=True, loc="upper right",
        frameon=True, fontsize=6, ncol=1
    )
    legend = styles.graph_legend_title(legend, fontsize=6)

    ax2.set_xlim([0, 365])
    ax2.set_ylim([0, 2.25])
    styles.heading(
        ax=ax2, letter="B. ",
        heading="MF6 API Agricultural Water Use irrigation", fontsize=6
    )

    plt.savefig("prudic_voronoi_model.tiff")
    plt.show()