import os
import flopy
from flopy.plot import styles
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("TKAgg")


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


data_ws = os.path.join("..", "..", "data")
nwt_ws = os.path.join(data_ws, "nwt_prudic_ag")
mf6_ws = os.path.join(data_ws, "mf6_prudic_voronoi")

nwt_file = os.path.join(nwt_ws, "prudic_ag.diversion9.txt")
mf6_file = os.path.join(mf6_ws, "prudic_vor_ag.out")

nwt_ag = pd.read_csv(nwt_file, delim_whitespace=True)
mf6_ag = pd.read_csv(mf6_file, delim_whitespace=True)

nwt_cbc = os.path.join(nwt_ws, "prudic_ag.cbc")
nwt_cbc = flopy.utils.CellBudgetFile(nwt_cbc)

conv = 0.0283168 # ft3 to m3
well_recs = nwt_cbc.get_data(text="AG WE")

mf6_ag = mf6_ag.groupby(by=["pid", "pkg", "kstp"], as_index=False).sum()
# mf6_ag_div = mf6_ag[mf6_ag.pkg == "sfr_0"]
# mf6_ag_div.reset_index(drop=True, inplace=True)

mf6_ag_lak = mf6_ag[mf6_ag.pkg == "lak_0"]
mf6_ag_lak.reset_index(drop=True, inplace=True)


total_mf6_div = mf6_ag_lak.q_from_provider.values # + mf6_ag_div.q_from_provider.values

# print(len(mf6_ag_div), len(mf6_ag_div.q_from_provider.values))
mf6_ag_pump = mf6_ag[mf6_ag.pkg == "maw_0"]
mf6_ag_pump = mf6_ag_pump.groupby(by=["pkg", "kstp"], as_index=False)["q_from_provider"].sum()

nwt_total = (nwt_ag["SW-DIVERSION"].values + nwt_ag["SUP-PUMPING"].values) * conv
mf6_total = (total_mf6_div + mf6_ag_pump.q_from_provider.values) * conv

tmp = np.array([0]*91 + list(mf6_total))
err = np.sum(nwt_total[:361] - np.array([0]*91 + list(mf6_total)))

print(np.sum(nwt_total[:361]))
print(np.sum(mf6_total))

print(err)
print((err / np.sum(nwt_total[:361])) * 100)
with styles.USGSPlot():
    mpl.rcParams["ytick.labelsize"] = 6
    mpl.rcParams["xtick.labelsize"] = 6
    mpl.rcParams["figure.dpi"] = 170
    fig, ax = plt.subplots(figsize=cm2inch(8.25, 10.25))
    ax.plot(
        mf6_ag_lak.kstp.values + 1,
        mf6_total,
        ":",
        zorder=2,
        label="MF6 total irrigation",
        lw=2,
        color="skyblue"
    )
    ax.plot(
        range(1, len(nwt_ag) + 1),
        nwt_total,
        zorder=1,
        label="NWT total irrigation",
        lw=2.5,
        color="dimgray"
    )

    ax.plot(
        mf6_ag_lak.kstp.values + 1,
        mf6_ag_lak.q_from_provider.values * conv,
        ":",
        zorder=4,
        label="MF6 lak water irrigation",
        lw=2.5,
        color="dodgerblue"
    )
    """
    ax.plot(
        mf6_ag_div.kstp.values + 1,
        total_mf6_div * conv,
        "-",
        zorder=4,
        label="MF6 total surface water irrigation",
        lw=2.5,
        color="limegreen"
    )
    """
    ax.plot(
        range(1, len(nwt_ag) + 1),
        nwt_ag["SW-DIVERSION"].values * conv,
        zorder=3,
        label="NWT surface water irrigation",
        lw=2,
        color="k"
    )

    styles.xlabel(ax=ax, label="Timestep", fontsize=7)
    styles.ylabel(ax=ax, label="Applied irrigation, in " + r"$m^{3}/day$", fontsize=7)
    styles.graph_legend(ax=ax, loc=1, fancybox=True, shadow=True, frameon=True,
                        fontsize=6, ncol=1)
    ax.set_xlim([0, 365])
    ax.set_ylim([0, 2.25])
    plt.savefig("prudic_results.tiff")
    plt.show()

