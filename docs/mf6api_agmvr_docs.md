# Guide to the mf6api_agmvr package

The mf6api_agmvr package simulates irrigated agricultural processes by 
extending the functionality of MODFLOW-6 through MODFLOW’s Extended Model 
Interface (XMI). The package is written in pure python, is tightly coupled with 
MODFLOW-6, and calculates irrigation demand and applied irrigation water at the 
outer iteration of a model’s time-step. The package also interfaces with FloPy 
(Bakker and others, 20xx; Bakker and others, 2022) and provides FloPy support 
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

The AGMVR Package is based on the calculation of irrigation demand and is 
limited by the available water that can be moved from a provider to receiver 
nodes. Irrigation demand is calculated by minimizing the difference between 
potential crop evapotranspiration $\left( ET_o K_c \right)$ and actual crop evapotranspiration 
( $ET_a$ ). The total volume of water demanded and consumed by a given crop under 
efficient conditions ( $Q_ET$ ) is calculated using:

