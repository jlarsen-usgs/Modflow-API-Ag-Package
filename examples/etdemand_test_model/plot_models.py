import flopy as fp
import matplotlib.colors
from flopy.plot import styles
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.collections import PathCollection
from matplotlib.path import Path
import matplotlib.patches as mpatches
import matplotlib as mpl
mpl.use("TKAgg")


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


def plot_mvr(ax, gwf, provider, kper, arrows=False, head_width=10, **kwargs):
    modelgrid = gwf.modelgrid
    xc = modelgrid.xcellcenters.ravel()
    yc = modelgrid.ycellcenters.ravel()
    parray = np.zeros(xc.shape)
    rarray = np.zeros(xc.shape)

    recarray = gwf.mvr.perioddata.data[kper]
    idx = np.where(recarray["pname1"] == provider)[0]
    recarray = recarray[idx]

    if "sfr" in provider:
        providerdata = gwf.sfr.packagedata.array
    elif "wel" in provider:
        providerdata = gwf.wel.stress_period_data.data[kper]
    else:
        return
    nodes = []
    for rid in recarray["id1"]:
        if "sfr" in provider:
            idx = np.where(providerdata["rno"] == rid)[0]
            rec = providerdata[idx]
        else:
            rec = providerdata[[rid]]
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

    parray[nodes] = 1
    pxy = [[xc[i], yc[i]] for i in nodes]

    rarray[recarray["id2"]] = 1
    rxy = [[xc[i], yc[i]] for i in recarray["id2"]]
    parray[parray != 1] = np.nan
    rarray[rarray != 1] = np.nan

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
        if "sfr" in provider:
            color = "darkblue"
        else:
            color = "k"
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
        if "sfr" in provider:
            color = "darkblue"
        else:
            color = "k"
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
        patches.append(mpatches.Patch(color="darkred"))
        labels.append("AG provider well")
    if sfr:
        patches.append(mpatches.Patch(color="c"))
        labels.append("AG provider reach")

    patches.append(mpatches.Patch(color="darkgreen"))
    labels.append("AG cells")
    if sfr:
        patches.append(mpatches.Patch(color="teal"))
        labels.append("SFR cells")

    if wells:
        patches.append(
            ax.scatter([-100], [-100], c='k', marker=r'$\leftarrow$', s=100)
        )
        labels.append("Irrigation from well")
    if sfr:
        patches.append(
            ax.scatter([-100], [-100], c="darkblue", marker=r'$\leftarrow$', s=100)
        )
        labels.append("Irrigation from SFR")


    legend = ax.legend(
        patches, labels, fancybox=True, shadow=True, loc=3, frameon=True, fontsize=7
    )
    legend = styles.graph_legend_title(legend, fontsize=7)


mf6_ws = os.path.join("..", "..", "data", "mf6_etdemand_test_problems")
nwt_ws = os.path.join("..", "..", "data", "nwt_etdemand_test_problems")

sim = fp.mf6.MFSimulation.load(sim_ws=os.path.join(mf6_ws, "etdemand_sup"))
gwf = sim.get_model("etdemand_sup")

ml = fp.modflow.Modflow.load(
    "etdemand_sup.nam",
    version="mfnwt",
    model_ws=nwt_ws
)
"""
with styles.USGSMap():
    mpl.rcParams["ytick.labelsize"] = 9
    mpl.rcParams["xtick.labelsize"] = 9
    fig, ax = plt.subplots(figsize=(8, 8))
    pmv = fp.plot.PlotMapView(model=ml)
    pmv.plot_bc("sfr")
    nwt_plot_diversion(ax, ml, "sfr", kper=0)
    nwt_plot_diversion(ax, ml, "sup", kper=0)
    pmv.plot_grid()
    nwt_plot_diversion(ax, ml, "sfr", kper=0, arrows=True)
    nwt_plot_diversion(ax, ml, "sup", kper=0, arrows=True)
    plot_explanation(ax)
    styles.heading(
        ax=ax, heading="MF-NWT Agricultural supplemental pumping example"
    )
    styles.xlabel(ax=ax, label="Model x-coordinate", fontsize=11)
    styles.ylabel(ax=ax, label="Model y-coordinate", fontsize=11)
    plt.xticks()
    plt.show()
    plt.show()


with styles.USGSMap():
    mpl.rcParams["ytick.labelsize"] = 9
    mpl.rcParams["xtick.labelsize"] = 9
    fig, ax = plt.subplots(figsize=(8, 8))
    pmv = fp.plot.PlotMapView(model=ml)
    # pmv.plot_bc("sfr")
    # nwt_plot_diversion(ax, ml, "sfr", kper=0)
    nwt_plot_diversion(ax, ml, "sup", kper=0)
    pmv.plot_grid()
    # nwt_plot_diversion(ax, ml, "sfr", kper=0, arrows=True)
    nwt_plot_diversion(ax, ml, "sup", kper=0, arrows=True)
    plot_explanation(ax, sfr=False)
    styles.heading(
        ax=ax, heading="MF-NWT Agricultural pumping example"
    )
    styles.xlabel(ax=ax, label="Model x-coordinate", fontsize=11)
    styles.ylabel(ax=ax, label="Model y-coordinate", fontsize=11)
    plt.xticks()
    plt.show()
    plt.show()


with styles.USGSMap():
    mpl.rcParams["ytick.labelsize"] = 9
    mpl.rcParams["xtick.labelsize"] = 9
    fig, ax = plt.subplots(figsize=(8, 8))
    pmv = fp.plot.PlotMapView(model=ml)
    pmv.plot_bc("sfr")
    nwt_plot_diversion(ax, ml, "sfr", kper=0)
    pmv.plot_grid()
    nwt_plot_diversion(ax, ml, "sfr", kper=0, arrows=True)
    plot_explanation(ax, wells=False)
    styles.heading(
        ax=ax, heading="MF-NWT Agricultural streamflow diversion example"
    )
    styles.xlabel(ax=ax, label="Model x-coordinate", fontsize=11)
    styles.ylabel(ax=ax, label="Model y-coordinate", fontsize=11)
    plt.xticks()
    plt.show()
    plt.show()


with styles.USGSMap():
    mpl.rcParams["ytick.labelsize"] = 9
    mpl.rcParams["xtick.labelsize"] = 9
    fig, ax = plt.subplots(figsize=(8, 8))
    pmv = fp.plot.PlotMapView(model=gwf)
    pmv.plot_bc("sfr_0")
    pmv.plot_bc("wel_0")
    plot_mvr(ax, gwf, "sfr_0", 0)
    plot_mvr(ax, gwf, "wel_0", 0)
    pmv.plot_grid()
    plot_mvr(ax, gwf, "sfr_0", 0, arrows=True)
    plot_mvr(ax, gwf, "wel_0", 0, arrows=True)
    plot_explanation(ax)
    styles.heading(
        ax=ax, heading="MF6 API Agricultural supplemental pumping example"
    )
    styles.xlabel(ax=ax, label="Model x-coordinate", fontsize=11)
    styles.ylabel(ax=ax, label="Model y-coordinate", fontsize=11)
    plt.xticks()
    plt.show()

with styles.USGSMap():
    mpl.rcParams["ytick.labelsize"] = 9
    mpl.rcParams["xtick.labelsize"] = 9
    fig, ax = plt.subplots(figsize=(8, 8))
    pmv = fp.plot.PlotMapView(model=gwf)
    pmv.plot_bc("wel_0")
    plot_mvr(ax, gwf, "wel_0", 0)
    pmv.plot_grid()
    plot_mvr(ax, gwf, "wel_0", 0, arrows=True)
    plot_explanation(ax, sfr=False)
    styles.heading(
        ax=ax, heading="MF6 API Agricultural pumping example"
    )
    styles.xlabel(ax=ax, label="Model x-coordinate", fontsize=11)
    styles.ylabel(ax=ax, label="Model y-coordinate", fontsize=11)
    plt.xticks()
    plt.show()

with styles.USGSMap():
    mpl.rcParams["ytick.labelsize"] = 9
    mpl.rcParams["xtick.labelsize"] = 9
    fig, ax = plt.subplots(figsize=(8, 8))
    pmv = fp.plot.PlotMapView(model=gwf)
    pmv.plot_bc("sfr_0")
    plot_mvr(ax, gwf, "sfr_0", 0)
    pmv.plot_grid()
    plot_mvr(ax, gwf, "sfr_0", 0, arrows=True)
    plot_explanation(ax, wells=False)
    styles.heading(ax=ax,
                   heading="MF6 API Agricultural streamflow diversion example")
    styles.xlabel(ax=ax, label="Model x-coordinate", fontsize=11)
    styles.ylabel(ax=ax, label="Model y-coordinate", fontsize=11)
    plt.xticks()
    plt.show()
"""

with styles.USGSMap():
    mpl.rcParams["ytick.labelsize"] = 6
    mpl.rcParams["xtick.labelsize"] = 6
    mpl.rcParams["figure.dpi"] = 170
    fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=cm2inch(17.15, 10.16))
    pmv = fp.plot.PlotMapView(model=gwf, ax=ax0)
    pmv.plot_bc("sfr_0")
    pmv.plot_bc("wel_0")
    plot_mvr(ax0, gwf, "sfr_0", 0)
    plot_mvr(ax0, gwf, "wel_0", 0)
    pmv.plot_grid(lw=0.5)
    plot_mvr(ax0, gwf, "sfr_0", 0, arrows=True)
    plot_mvr(ax0, gwf, "wel_0", 0, arrows=True)
    plot_explanation(ax0)
    styles.heading(
        letter="A",
        ax=ax0,
        heading="MF6 API Agricultural Water Use Package",
        fontsize=7
    )
    styles.xlabel(ax=ax0, label="Model x-coordinate", fontsize=7)
    styles.ylabel(ax=ax0, label="Model y-coordinate", fontsize=7)

    pmv = fp.plot.PlotMapView(model=ml, ax=ax1)
    pmv.plot_bc("sfr")
    nwt_plot_diversion(ax1, ml, "sfr", kper=0)
    nwt_plot_diversion(ax1, ml, "sup", kper=0)
    pmv.plot_grid(lw=0.75)
    nwt_plot_diversion(ax1, ml, "sfr", kper=0, arrows=True)
    nwt_plot_diversion(ax1, ml, "sup", kper=0, arrows=True)
    # plot_explanation(ax1)
    styles.heading(
        letter="B",
        ax=ax1,
        heading="NWT Agricultural Water Use Package",
        fontsize=7
    )
    styles.xlabel(ax=ax1, label="Model x-coordinate", fontsize=7)
    # styles.ylabel(ax=ax, label="Model y-coordinate", fontsize=11)

    plt.tight_layout()
    plt.show()