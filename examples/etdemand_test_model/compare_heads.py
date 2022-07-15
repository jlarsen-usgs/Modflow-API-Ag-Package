import flopy as fp
from flopy.plot import styles
from flopy.plot import PlotMapView
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import os
import numpy as np


def head_from_list_file(f):
    with open(f) as foo:
        read_heads = False
        header = False

        heads = []
        for line in foo:
            line = line.strip()
            if "HEAD IN LAYER" in line:
                header = True
            elif header:
                if "........" in line:
                    header = False
                    read_heads = True
                    ts_head = []
            elif read_heads and line == "":
                read_heads = False
                heads.append([ts_head])
            elif read_heads:
                t = line.strip().split()
                ts_head.append([float(i) for i in t[1:]])
            else:
                pass

    return np.array(heads)


def plot_heads(nwt, mf6, model, plot_nwt=True, vmin=None, vmax=None):
    mf6_hds = os.path.join(mf6, model, f"{model}.hds")
    mf6_cbc = os.path.join(mf6, model, f"{model}.cbc")
    nwt_hds = os.path.join(nwt, f"{model}.hds")
    nwt_cbc = os.path.join(nwt, f"{model}.cbc")
    nwt_nam = f"{model}.nam"
    grid_file = os.path.join(mf6, model, f"{model}.dis.grb")

    nwt_ml = fp.modflow.Modflow.load(nwt_nam, model_ws=nwt, version="mfnwt")
    mf6_hds = fp.utils.HeadFile(mf6_hds)
    mf6_cbc = fp.utils.CellBudgetFile(mf6_cbc)
    nwt_cbc = fp.utils.CellBudgetFile(nwt_cbc)
    nwt_hds = fp.utils.HeadFile(nwt_hds)
    grb = fp.mf6.utils.MfGrdFile(grid_file)
    modelgrid = grb.modelgrid

    mf6_head = mf6_hds.get_alldata()[-1]
    nwt_head = nwt_hds.get_alldata()[-1]
    spdis = mf6_cbc.get_data(text="DATA-SPDIS")[-1]
    qx = spdis["qx"].reshape(modelgrid.shape)
    qy = spdis["qy"].reshape(modelgrid.shape)

    frf = nwt_cbc.get_data(text="FLOW RIGHT FACE")[-1]
    fff = nwt_cbc.get_data(text="FLOW FRONT FACE")[-1]
    flf = np.zeros(fff.shape)
    qqx, qqy, qqz = fp.utils.postprocessing.get_specific_discharge(
        (frf, fff, flf), nwt_ml
    )

    if vmin is None and vmax is None:
        if plot_nwt:
            vmin = 10000000
            vmax = -10000000
            if np.min(mf6_head) < vmin:
                vmin = np.min(mf6_head)
            if np.min(nwt_head) < vmin:
                vmin = np.min(nwt_head)

            if np.max(mf6_head) > vmax:
                vmax = np.max(mf6_head)
            if np.max(nwt_head) > vmax:
                vmax = np.max(nwt_head)
        else:
            vmin, vmax = np.min(mf6_head), np.max(mf6_head)

    vmin, vmax = None, None

    ttext = "irrigation well"
    if model == "etdemand_div":
        ttext = "surface water irrigation"
    elif model =="etdemand_sup":
        ttext = "conjunctive use"

    if plot_nwt:
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 5))
    else:
        fig, ax0 = plt.subplots(figsize=(8, 8))

    with styles.USGSMap():
        pmv = PlotMapView(modelgrid=modelgrid)
        pc = pmv.plot_array(mf6_head, ax=ax0, cmap="plasma", vmin=vmin, vmax=vmax)
        cs = pmv.contour_array(mf6_head, ax=ax0, linewidths=0.5, linestyles="-", colors="black", zorder=3) #, levels=np.arange(vmin, vmax, 10))
        ax0.clabel(cs, fmt="%.3f")
        lc1 = pmv.plot_grid(ax=ax0, zorder=2)
        qv = pmv.plot_vector(qx, qy, ax=ax0, normalize=True, color="0.75")
        ax0.set_aspect("equal")
        styles.heading(ax=ax0, letter="A", heading=f"MF6 {ttext} example")

        if plot_nwt:
            pmv = PlotMapView(modelgrid=modelgrid)
            pc = pmv.plot_array(nwt_head, ax=ax1, cmap="plasma", vmin=vmin,
                                vmax=vmax)
            cs = pmv.contour_array(nwt_head, ax=ax1, linewidths=0.5,
                                   linestyles="-", colors="black",
                                   zorder=3)
            ax0.clabel(cs, fmt="%.3f")
            lc1 = pmv.plot_grid(ax=ax1, zorder=2)
            qv = pmv.plot_vector(qqx, qqy, ax=ax1, normalize=True, color="0.75")
            ax1.set_aspect("equal")
            styles.heading(ax=ax1, letter="B",heading=f"MF-NWT {ttext} example")

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        cbar = fig.colorbar(pc, cax=cbar_ax, label="Hydraulic Head")

        # create proxy artists
        lc = LineCollection([[[0, 0],[1, 1]],], colors="k", linestyle="solid")
        qv = ax0.scatter([-1], [-1], c="0.75", marker=r'$\leftarrow$', s=100)
        leg = fig.legend(
            [lc, qv],
            ["hydraulic head contour line", "flow direction vector"],
            fancybox=True,
            shadow=True,
            loc=8,
            frameon=True,
            ncol=2
        )
        styles.graph_legend_title(leg)
        plt.show()

    plt.imshow(mf6_head[0] - nwt_head[0], interpolation=None)
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    mf6_ws = os.path.join("..", "..", "data", "mf6_etdemand_test_problems")
    nwt_ws = os.path.join("..", "..", "data", "nwt_etdemand_test_problems")

    models = ["etdemand_well", "etdemand_div", "etdemand_sup"]
    plot_heads(nwt_ws, mf6_ws, models[0], plot_nwt=True, )
