# Guide to the mf6api_agmvr package

The mf6api_agmvr package simulates irrigated agricultural processes by 
extending the functionality of MODFLOW-6 through MODFLOW’s Extended Model 
Interface (XMI). The package is written in pure python, is tightly coupled with 
MODFLOW-6, and calculates irrigation demand and applied irrigation water at the 
outer iteration of a model’s time-step. The package also interfaces with FloPy 
(Bakker and others, 2016; Bakker and others, 2022) and provides FloPy support 
for constructing Agricultural water mover (AGMVR) input files.

This guide covers the installation process, input file construction and 
formatting, and provides links to example problems.

## Installation
The mf6api_agmvr package can be installed from the command line or terminal 
window by running the following command:
```commandline
python -m pip install https://github.com/jlarsen-usgs/mf6api_agmvr/archive/refs/heads/main.zip
```

## The Agricultural Mover Package (AGMVR)
The AGMVR package can be used to calculate irrigation demand and transfer water 
from an agricultural provider to agricultural receiver nodes. Providers are 
extraction wells, streamflow routing reaches, and lakes that have water 
available for extraction. The list of packages that can provide water to the 
AGMVR package are
   - Well Package (WEL)
   - Multi-Aquifer Well Package (MAW)
   - Streamflow Routing Package (SFR)
   - Lake Package (LAK)

Receiver nodes are unsaturated zone packages that simulate evapotranspiration 
and solve a continuity equation of inflows, outflows, and change in storage. 
These nodes receive water from a provider or multiple provider(s). The only 
receiver package supported is:
   - Unsaturated Zone Flow Package (UZF)

The program will terminate with an error if the AGMVR is used with an 
unsupported package type.

The AGMVR Package is based on the calculation of irrigation demand presented by 
Niswonger (2020) and is limited by the available water that can be moved from a 
provider to receiver nodes. Irrigation demand is calculated by minimizing the 
difference between potential crop evapotranspiration $\left( ET_o K_c \right)$ 
and actual crop evapotranspiration $\left( ET_a \right)$. The total volume of 
water demanded and consumed by a given crop under efficient conditions 
$\left( Q_{ET} \right)$ is calculated using:

$$Q_{ET} = \displaystyle\sum_{n=1}^{ncell} ET_o K_c A_n$$ 

where $ET_o$ is the potential evapotranspiration flux, $K_c$ is the crop 
coefficient, and $A_n$ is the area of the model cell. The amount of water that 
must be supplied by provider nodes to meet this condition is calculated by 
minimizing the ET deficit $\left( ET_d \right)$.

$$min \left( ET_d \right) = ET_o K_c - ET_a$$

Niswonger (2020) presented a volumetric solution to this minimization problem. 
This equation has been adapted and modified to allow for the simulation of 
deficit and flush irrigation by including an application efficiency factor 
$\left( e_a \right)$ as follows:

$$Q_{c,i+1} = \left( Q_{c,i} + \frac{Q_{ET,i+1} - Q_{ET,i}}{\delta Q_{ET,i} / \delta Q_{c,i}} \right) * e_a $$

where i is the outer iteration counter for a given MODFLOW timestep and $Q_c$ 
is the calculated volumetric water demand with adjustment for application 
efficiency. The amount of water available to be moved from a provider to 
receiver nodes $\left( Q_A \right)$ is then calculated as:

$$Q_A =
    \begin{cases}
      Q_{c,i+1} & \quad Q_{p} > Q_{c,i+1}\\
      Q{p} & \quad Q{p} \leq  Q_{c,i+1}
    \end{cases}
$$

where Q_p is the amount of water available from a provider.

Once available water for irrigation has been calculated, the actual amount of 
irrigation that is applied to receiver nodes $\left( q_a \right)$ is calculated 
from the equation:

$$q_a = \frac{Q_A * e_i}{A_n}$$

The irrigation efficiency factor $\left( e_i \right)$ is used to adjust for 
inefficient irrigation methods. Irrigation losses $\left( Q_L \right)$ due to 
canopy interception or other inefficiencies is calculated by

$$Q_L = Q_A - (Q_A * e_a)$$

and is removed from the model.

Input to the Agricultural Water Mover (AGMVR) Package is read from the file 
type “AGMVR” in the Name File by the “ModflowGwfagmvr” python class and 
processed by the “ModflowAgmvr” class.

### Structure of Blocks
_FOR EACH SIMULATION_

BEGIN **OPTIONS**  
[ETDEMAND]  
END **OPTIONS**


BEGIN **DIMENSIONS**  
MAXMVR <maxmvr<n>>  
MAXPACKAGES <maxpackages<n>>  
END **DIMENSIONS**


BEGIN **PACKAGES**  
<pname<n>>  
<pname<n>>  
....  
END **PACKAGES**


_FOR ANY STRESS PERIOD_

BEGIN **PERIOD**  
<pname1<n>> <id1<n>> <pname2<n>> <id2<n>> <irr_eff<n>> <app_eff<n>>  
<pname1<n>> <id1<n>> <pname2<n>> <id2<n>> <irr_eff<n>> <app_eff<n>>  
....  
END **PERIOD**

Note: All of the information supplied within a **PERIOD** block will apply 
only for the stress period that it is specified for. If no **PERIOD** block is 
specified for a specific stress period, no agricultural processes will be 
simulated in that stress period. This behavior differs from the behavior of 
the standard MVR package.

### Explanation of Variables
#### Block: OPTIONS
ETDEMAND — optional keyword to simulate agricultural processes using the 
evapotranspiration deficit calculation. This keyword is active by default even 
if not supplied by the user.

#### Block: DIMENSIONS
maxmvr — integer value specifying the maximum number of water mover entries that 
will be specified for any stress period

maxpackages — integer value specifying the maximum number of unique packages that 
are included in this agricultural water mover package

#### Block: PACKAGES
pname — the name of the package that may be included in a subsequent stress 
period block. The package name is assigned in the name file for the groundwater 
flow model. Package names are optionally provided in the name file. If they are 
not provided by the user, then packages are assigned a default value, which is 
the package acronym followed by a hyphen and the package number. For example, 
the first WELL Package is named WEL-1.

#### Block: PERIOD
iper — integer value specifying the stress period number for which the data 
specified in the PERIOD block apply. IPER must be less than or equal to NPER in 
the TDIS package and greater than zero. The IPER value assigned to a stress 
period block must be greater than the IPER value assigned for the previous 
PERIOD block.

pname1 — package name for the provider. The package PNAME1 must be designated 
to provide water through the MVR Package by specifying the keyword “MOVER” in 
the package’s OPTIONS block. 

id1 — is the identifier for the provider. For the standard boundary packages, 
the provider identifier is the number of the boundary as it is listed in the 
package input file. (Note that the order of these boundaries may change by 
stress period, which must be accounted for in the Mover Package.) The first 
well has an identifier of one, the second is two, and so forth. For the 
advanced packages, the identifier is the reach number (SFR Package) or well 
number (MAW Package). For the Lake Package, ID1 is the lake outlet number.

pname2 — package name for the receiver. Note: The UZF package is the only 
package supported at this time.

id2 — is the identifier for the receiver. For the UZF package this is the UZF 
variable IUZNO.

value — is the maximum value of water to move from a provider to a receiver 
node. This is specified as a volumetric flow rate.

eff_fact — is the irrigation efficiency factor. This factor can range from 0 
to 1, where 1 represents 100% efficiency in irrigation, which translates to no 
irrigation losses.

app_fact — is the application efficiency factor. Values greater than 1 can be 
used to simulate the application of irrigation water beyond the crop’s 
requirement based on the ET deficit (e.g., flush irrigation). Values less than 
1 can be used to simulate deficit irrigation.   

### Example Input File
<p align="center">
<img src="https://raw.githubusercontent.com/jlarsen-usgs/mf6api_agmvr/main/docs/agmvr.png" alt="input file"/>
</p>

### Building an AGMVR package with Python
Installation of the mf6api_package dynamically links with the FloPy package to 
create input reading and writing routines for the AGMVR package. Example code 
is presented here that shows how to build an AGMVR package and add it to an 
existing MODFLOW model. 

   * import the required packages and set the path to the existing model
```python
import os
import flopy


# point to the flopy example model UZF_3lay
flopy_ws = os.path.split(flopy.__file__)[0]
sim_ws = os.path.join(
    flopy_ws, "..", "examples", "data", "mf6", "test001e_UZF_3lay"
)
```
   * Load the existing model with FloPy and enable the mover option in the UZF 
package
```python
# load the existing model
sim = flopy.mf6.modflow.MFSimulation.load(sim_ws=sim_ws)
gwf = sim.get_model("gwf_1")

# enable the mover option in UZF
gwf.uzf.mover = True
```
   * Add a WEL package for agricultural providers and enable the mover option
```python
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
```
   * Build the AGMVR package using FloPy functions
```python
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
```
   * Write the simulation to a new directory
```python
# set new simulation path and write model to file
out_ws = os.path.join("..", "data", "docs_example")
sim.set_sim_path(out_ws)
sim.write_simulation()
```

### Running MODFLOW with the AGMVR package
The AGMVR package cannot be run with the traditional MODFLOW-6 executable 
because the calculation code has been implemented in python for use with the 
MODFLOW-6 XMI. As a result, a model that uses the AGMVR package must be run in 
python using the mf6api_agmvr package. Once a simulation has been loaded into 
FloPy the user can pass the simulation object to the ModflowAgmvr class and 
run the model.
```python
from mf6api_agmvr import ModflowAgmvr

mf6ag = ModflowAgmvr(sim, mvr_name="agmvr")
mf6ag.run_model()
```

### Reading output from the AGMVR package
The AGMVR package writes a single ascii output file upon completion of the 
model. This file is named {model_name}_ag.out, where model_name is the name
of the model specified in the MODFLOW6 sim file that has been run. For example,
if the name of the model is "gwf_1" the output file name will be 
"gwf_1_ag.out". The output file contains information about the potential et,
actual et, demand, the volumetric flux (units, $\left( L^{3}/T \right)$ ) 
removed from a provider, and the volumetric flux 
(units, $\left( L^{3}/T \right)$ ) supplied to receiver nodes. This output is 
written for every time step.

Output can be loaded and processed in python using pandas
```python
import pandas as pd

output_file = os.path.join(out_ws, "gwf_1_ag.out")
ag_df = pd.read_csv(output_file, delim_whitespace=True)
```

### Example problems
Example problems are distributed with the mf6api_agmvr repository and can be 
found in the [examples directory](https://github.com/jlarsen-usgs/mf6api_agmvr/tree/main/examples).

### References
Bakker, Mark, Post, Vincent, Langevin, C. D., Hughes, J. D., White, J. T., 
Starn, J. J. and Fienen, M. N., 2016, Scripting MODFLOW Model Development Using 
Python and FloPy: Groundwater, v. 54, p. 733–739, doi:10.1111/gwat.12413.

Bakker, Mark, Post, Vincent, Hughes, J. D., Langevin, C. D., White, J. T., 
Leaf, A. T., Paulinski, S. R., Bellino, J. C., Morway, E. D., Toews, M. W., 
Larsen, J. D., Fienen, M. N., Starn, J. J., and Brakenhoff, Davíd, 2022, 
FloPy v3.3.5: U.S. Geological Survey Software Release, 07 March 2022, 
http://dx.doi.org/10.5066/F7BK19FH

Niswonger, R. G., 2020, An Agricultural Water Use Package for MODFLOW and 
GSFLOW. Environmental Modelling and Software 125. 
doi: 10.1016/j.envsoft.2019.104617