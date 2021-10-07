# quick and dirty script to calculate monthly mean ppt and eto for Davis CIMIS
# station

import numpy as np
import pandas as pd

deto = {i: [] for i in range(1, 13)}
dppt = {i: [] for i in range(1, 13)}
lut = {'jan': 1,
       'feb': 2,
       'mar': 3,
       'apr': 4,
       "may": 5,
       "jun": 6,
       "jul": 7,
       "aug": 8,
       "sep": 9,
       "oct": 10,
       "nov": 11,
       "dec": 12}

with open('cimis_davis_monthly.csv') as foo:
    for ix, line in enumerate(foo):
        t = line.strip().split(",")
        if ix == 0:
            continue
        else:
            try:
                mo = t[3].split()[0].lower()
                mo = lut[mo]
            except:
                continue

            try:
                eto = float(t[4]) * 0.001
            except:
                eto = np.nan

            try:
                ppt = float(t[6]) * 0.001
            except:
                ppt = np.nan

            deto[mo].append(eto)
            dppt[mo].append(ppt)

    d = {'month': [], "eto_avg_m": [], "ppt_avg_m": []}
    for i in range(1, 13):
        d["month"].append(i)
        d["eto_avg_m"].append(np.nanmean(deto[i]))
        d["ppt_avg_m"].append(np.nanmean(dppt[i]))

    df = pd.DataFrame.from_dict(d)
    df.to_csv("davis_monthly_ppt_eto.txt", index=False)