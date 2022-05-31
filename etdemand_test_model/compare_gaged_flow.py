import os

import matplotlib.pyplot as plt
import matplotlib as mpl
from flopy.plot import styles
import numpy as np
import pandas as pd
from scipy.stats import linregress


def plot_gaged_flows(nwt, mf6, model):
    nwt_gage = os.path.join(nwt, f"{model}_gage10.txt")
    mf6_gage = os.path.join(mf6, f"{model}.obs.out")

    nwt_gage = pd.read_csv(
        nwt_gage,
        delim_whitespace=True,
        names=["time", "stage", "flow"],
        skiprows=2
    )

    mf6_gage = pd.read_csv(
        mf6_gage,
        names=["time", "flow"],
        skiprows=1
    )

    r2 = linregress(nwt_gage.flow.values, np.abs(mf6_gage.flow.values)[1:])[2] ** 2

    if model == "etdemand_div":
        ttext = "surface water irrigation"
    else:
        ttext = "conjunctive use"

    with styles.USGSPlot():
        mpl.rcParams["ytick.labelsize"] = 9
        mpl.rcParams["xtick.labelsize"] = 9
        fig, ax = plt.subplots(figsize=(8, 8))

        lc0 = ax.plot(nwt_gage.time.values, nwt_gage.flow.values, label="MF-NWT simulated flows", color="k", ls="-", lw=2)
        lc1 = ax.plot(mf6_gage.time.values, np.abs(mf6_gage.flow.values), label="MF6 simulated flows", color="lightblue", ls="--", lw=2.5)
        styles.ylabel(ax=ax, label="Simulated streamflow, in " + r"$m^{3}/day$", fontsize=10)
        styles.xlabel(ax=ax, label="Day", fontsize=10)
        styles.heading(ax=ax, heading=f"Comparison of simulated outflows for {ttext} example")
        styles.add_text(ax=ax, text=r"$r^{2}$" + f" = {r2 :.2f}", x=0.5, y=0.88, fontsize=10)
        leg = ax.legend(fancybox=True, shadow=True, loc=0, frameon=True, fontsize=10)
        plt.ylim([15, 55])
        styles.graph_legend_title(leg)
        plt.show()


if __name__ == "__main__":
    mf6_ws = os.path.join("..", "data", "mf6_etdemand_test_problems")
    nwt_ws = os.path.join("..", "data", "nwt_etdemand_test_problems")

    models = ["etdemand_div", "etdemand_sup"]
    plot_gaged_flows(nwt_ws, mf6_ws, models[1])
