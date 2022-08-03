import flopy
import os
import pandas as pd
import matplotlib.pyplot as plt
from mf6api_agmvr import ModflowAgmvr


sim_ws = os.path.join("..", "data", "quickstart")
name = "quickstart"

nper = 1
perioddata = [(10, 10, 1.0) for _ in range(nper)]
sim = flopy.mf6.MFSimulation(name, sim_ws=sim_ws)
tdis = flopy.mf6.ModflowTdis(sim, nper=nper, perioddata=perioddata)
ims = flopy.mf6.ModflowIms(sim, complexity="COMPLEX")

nrow, ncol, delc, delr = 10, 10, 10, 10
top = 25
strt = 20
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
dis = flopy.mf6.ModflowGwfdis(gwf, nrow=nrow, ncol=ncol, delc=delc, delr=delr, top=top)
ic = flopy.mf6.ModflowGwfic(gwf, strt=strt)
npf = flopy.mf6.ModflowGwfnpf(gwf, icelltype=1, save_specific_discharge=True)
sto = flopy.mf6.ModflowGwfsto(gwf, iconvert=1)
oc = flopy.mf6.ModflowGwfoc(
    gwf,
    budget_filerecord=f"{name}.cbc",
    head_filerecord=f"{name}.hds",
    saverecord={i: [("HEAD", "ALL"), ("BUDGET", "ALL")] for i in range(nper)}
)

# build a UZF package
nuzfcells = 100
packagedata = []
cnt = 0
for i in range(nrow):
    for j in range(ncol):
        rec = (cnt, (0, i, j), 1, 0, 0.33, 8.64, 0.05, 0.35, 0.08, 5)
        packagedata.append(rec)
        cnt += 1

perioddata = {}
for i in range(nper):
    spd = []
    for node in range(nuzfcells):
        rec = (node, 1e-10, 8.0e-01, 6, 0.06, -1.1, -75, 1.0)
        spd.append(rec)
    perioddata[i] = spd

uzf = flopy.mf6.ModflowGwfuzf(
    gwf,
    simulate_et=True,
    unsat_etwc=True,
    linear_gwet=True,
    simulate_gwseep=True,
    nuzfcells=nuzfcells,
    packagedata=packagedata,
    perioddata=perioddata,
    mover=True
)

# create a well package to simulate AG pumping from wells
stress_period_data = {}
for per in range(nper):
    # cellid, maxq for ag pumping!
    stress_period_data[per] = [[(0, 4, 5), -100.],]

wel = flopy.mf6.ModflowGwfwel(
    gwf,
    stress_period_data=stress_period_data,
    mover=True
)

# finally build the AGMVR package
perioddata = {}
for per in range(nper):
    # rec: pak1, id1, pak2, id2, value, irr_eff, app_eff
    rec = ("wel_0", 0, "uzf_0", 55, 100, 1., 1.)
    perioddata[per] = [rec,]

mvr = flopy.mf6.ModflowGwfagmvr(
    gwf,
    maxmvr=1,
    maxpackages=2,
    packages=[("wel_0",), ("uzf_0",)],
    perioddata=perioddata
)

sim.write_simulation()

mfag = ModflowAgmvr(sim, mvr_name="agmvr")
mfag.run_model()

ag_out = os.path.join(sim_ws, f"{name}_ag.out")
ag_out = pd.read_csv(ag_out, delim_whitespace=True)

ag_out = plt.plot(ag_out.kstp.values, ag_out.q_to_receiver.values)
plt.ylabel("Agricultural irrigation provided, in " + r"$ft^{3}$")
plt.xlabel("Time step")
plt.show()
