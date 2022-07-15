[![ModflowAgmvr CI](https://github.com/jlarsen-usgs/modflow6api_agMVR/actions/workflows/ci.yml/badge.svg)](https://github.com/jlarsen-usgs/modflow6api_agMVR/actions/workflows/ci.yml)

# modflowapi_Agmvr
This repository contains the ModflowAgmvr class that interfaces with the modflowapi to simulate
irrigated agriculture with Modflow6. 

## Importing ModflowAgmvr
From the base of this repository: 

```python
import sys
sys.path.append("mf6api_agmvr")

from mf6_agmvr import ModflowAgmvr
```

## Software requirements
Python >= 3.7  
flopy >= 3.3.5 (`pip install flopy`)  
modflowapi (`pip install modflowapi`)  
numpy  
pandas  

## Authors
Joshua D. Larsen

## Version
1.0b
