import os
import flopy
from mf6api_agmvr import ModflowAgmvr
import pandas as pd


# point to the flopy example model UZF_3lay
flopy_ws = os.path.split(flopy.__file__)[0]
sim_ws = os.path.join(
    flopy_ws, "..", "examples", "data", "mf6", "test001e_UZF_3lay"
)

# load the existing model
sim = flopy.mf6.modflow.MFSimulation.load(sim_ws=sim_ws)
gwf = sim.get_model("gwf_1")

# enable the mover option in UZF
uzf = gwf.uzf
gwf.uzf.mover = True

# add a WEL package to the model and enable MVR
perioddata = {
    i : [
        ((2, 0, 2), -50),
        ((2, 0, 3), -50),
    ]
    for i in range(gwf.nper)
}

wel = flopy.mf6.modflow.ModflowGwfwel(
    gwf,
    mover=True,
    save_flows=True,
    stress_period_data=perioddata
)

# build the AG package
packages = [("wel_0",), ("uzf-1",)]
perioddata = {
    i : [
        ("wel_0", 0, "uzf-1", 3, 50, 1, 1),
        ("wel_0", 1, "uzf-1", 3, 50, 1, 1)
    ]
    for i in range(gwf.nper)
}

ag = flopy.mf6.modflow.ModflowGwfagmvr(
    gwf,
    maxmvr=len(packages),
    maxpackages=2,
    packages=packages,
    perioddata=perioddata
)

# set new simulation path and write model to file
out_ws = os.path.join("..", "data", "docs_example")
sim.set_sim_path(out_ws)
sim.write_simulation()

# create a ModflowAgmvr object and run model
mf6ag = ModflowAgmvr(sim, mvr_name="agmvr")
mf6ag.run_model()

output_file = os.path.join(out_ws, "gwf_1_ag.out")
ag_df = pd.read_csv(output_file, delim_whitespace=True)
print(ag_df.head())