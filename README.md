[![Modflow-API-Ag CI](https://github.com/jlarsen-usgs/mf6api_ag/actions/workflows/ci.yml/badge.svg)](https://github.com/jlarsen-usgs/mf6api_ag/actions/workflows/ci.yml)

# Modflow API Ag Package 
This repository contains the ModflowApiAg class that interfaces with the 
modflowapi to simulate irrigated agriculture with MODFLOW 6. The ModflowApiAg 
api package overrides the MODFLOW 6 water mover (MVR) package and calculates
irrigation requirement from potential and actual evapotranspiration in the 
soil zone. Provider packages supported by API-AG include WEL, SFR, MAW, and LAK.
The receiver package must be UZF. More information can be found in the 
[documentation](https://github.com/jlarsen-usgs/mf6api_ag/blob/main/docs/documentation.md) 
and example problems are located in the 
[examples directory](https://github.com/jlarsen-usgs/mf6api_ag/tree/main/examples) 
of this repository.

## Software requirements
Python >= 3.8  
flopy >= 3.3.6 (`pip install flopy`)  
modflowapi (`pip install modflowapi`)  
numpy  
pandas

## Supported Platforms
Windows 10  
Linux  
MacOS  

## Installing mf6api_ag
Open a command-line or anaconda prompt and run the following

```commandline
python -m pip install https://github.com/jlarsen-usgs/mf6api_ag/archive/refs/heads/main.zip
```

Alternatively, the user can download the repository cd into the main directory
and run the command
```commandline
python -m pip install .
```

## Importing ModflowApiAg

```python
from mf6api_ag import ModflowApiAg
```

## Quickstart guide
Quickstart usage for a single layer, single well AG MVR model

```python
import flopy
import os
import pandas as pd
import matplotlib.pyplot as plt
from mf6api_ag import ModflowApiAg

# build a new modflow-6 model 
sim_ws = os.path.join(".")
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
dis = flopy.mf6.ModflowGwfdis(gwf, nrow=nrow, ncol=ncol, delc=delc, delr=delr,
                              top=top)
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
    stress_period_data[per] = [[(0, 4, 5), -100.], ]

wel = flopy.mf6.ModflowGwfwel(
    gwf,
    stress_period_data=stress_period_data,
    mover=True
)

# finally build the AGMVR package
perioddata = {}
for per in range(nper):
    # rec: pak1, id1, pak2, id2, value (max volume of water applied), irr_eff, app_eff
    rec = ("wel_0", 0, "uzf_0", 55, 100, 1., 1.)
    perioddata[per] = [rec, ]

mvr = flopy.mf6.ModflowGwfapiag(
    gwf,
    maxmvr=1,
    maxpackages=2,
    packages=[("wel_0",), ("uzf_0",)],
    perioddata=perioddata
)

sim.write_simulation()

mfag = ModflowApiAg(sim, mvr_name="apiag")
mfag.run_model()

ag_out = os.path.join(sim_ws, f"{name}_ag.out")
ag_out = pd.read_csv(ag_out, delim_whitespace=True)

ag_out = plt.plot(ag_out.kstp.values, ag_out.q_to_receiver.values)
plt.ylabel("Agricultural irrigation provided, in " + r"$ft^{3}$")
plt.xlabel("Time step")
plt.show()
```

## Examples
Example problems can be found [here](https://github.com/jlarsen-usgs/Modflow-API-Ag-Package/tree/main/examples)

   - `docs_example.py`: example problem included in the documentation file
   - `quickstart.py`: example problem listed in the Quickstart section of readme.md
   - `etdemand_gwet_test_model`: folder that includes problems that test combined gwet and uzfet
calculations for the Modflow-API-Ag-Package.
     - `build_nwt_models.py`: builds modflow-nwt equivalent models for output comparison
     - `compare_ag_mvr_conjunctive.py`: builds and runs a conjunctive use (surface-water and groundwater) irrigation problem
and then compares output with the modflow-nwt AG package. This problem corresponds to example 1 in the 
accompanying Groundwater paper (Larsen et. al., xxxx; in review)
     - `compare_ag_mvr_sfr.py`: builds and runs a surface-water irrigation problem
and compares output with the modflow-nwt AG package
     - `compare_ag_mvr_well.py`: builds and runs a ground-water irrigation problem using the WEL package
and compares output with the modflow-nwt AG package
   - `etdemand_test_model`: folder that includes problems that test uzfet
calculations for the Modflow-API-Ag-Package
     - `compare_ag_mvr_conjunctive.py`: builds and runs a conjunctive use (surface-water and groundwater) irrigation problem
and then compares output with the modflow-nwt AG package
     - `compare_ag_mvr_sfr.py`: builds and runs a surface-water irrigation problem
and compares output with the modflow-nwt AG package
     - `compare_ag_mvr_well.py`: builds and runs a ground-water irrigation problem using the WEL package
and compares output with the modflow-nwt AG package
   - `flopy_agmvr_test_model`: folder that includes a problem that tests the
flopy compatible `ModflowGwfApiAg` class.
     - `compare_agmvr_conjunctive.py`: script that tests the flopy compatible
`ModflowGwfApiAg` class.
   - `lak_test_model`: folder that includes a test problem that tests surface-water
irrigation supplied by the LAK package
     - `ag_mvr_lak.py`: script that tests surface-water irrigation from the LAK package
   - `maw_test_model`: folder that includes a test problem that tests ground water
irrigation supplied by the MAW package
     - `ag_mvr_maw.py`: script that tests ground water irrigation from the MAW package
   - `prudic_model`: folder that contains scripts that build, run, and compare output
for the Green Valley model (Prudic et. al., 2004; Niswonger et. al., 2006; Niswonger et. al., 2020).
These example problems are presented as Example 2 in Larsen et. al. (20XX; in review).
     - `mf6_ag_prudic.py`: script that builds a MODFLOW6 API AG version of 
the Green Valley model described in Niswonger (2020) and compares output to 
the MODFLOW-NWT equivalent of the model.
     - `mf6_ag_prudic_scenarios.py`: script that builds three versions of the 
Green Valley model (Niswonger, 2020) and applies application efficiency factors
to the simulated irrigation 
     - `nwt_ag_prudic.py`: script that builds the MODFLOW-NWT version of the
Green Valley model and includes irrigation through the AG package (Niswonger, 2020)

## Documentation
Documentation can be found [here](https://github.com/jlarsen-usgs/mf6api_ag/blob/main/docs/documentation.md)

## Authors
Joshua D. Larsen

## Version
1.0.0

## IPDS number
IP-149615

## Suggested Citation
Larsen, J.D., 2023, Agricultural water use package for the MODFLOW API. 
U.S. Geological Survey Software Release. DOI:[10.5066/P9K6UW9F](https://doi.org/10.5066/P9K6UW9F)

## Disclaimer
This software is preliminary or provisional and is subject to revision. It is 
being provided to meet the need for timely best science. The software has not 
received final approval by the U.S. Geological Survey (USGS). No warranty, 
expressed or implied, is made by the USGS or the U.S. Government as to the 
functionality of the software and related material nor shall the fact of 
release constitute any such warranty. The software is provided on the condition 
that neither the USGS nor the U.S. Government shall be held liable for any 
damages resulting from the authorized or unauthorized use of the software.