import os
import flopy
import numpy as np
import pandas as pd
from math import log10, floor

sws = os.path.abspath(os.path.dirname(__file__))


def round_to_n(x, n):
    if x == 0:
        return 0
    t = round(x, -int(floor(log10(abs(x))) - (n - 1)))
    return t


def build_nwt_test_models(name):
    model_ws = os.path.join(sws, "..", "..", "data", "nwt_etdemand_gwet_test_problems")

    ml = flopy.modflow.Modflow(
        name,
        version="mfnwt",
        exe_name="mfnwt",
        model_ws=model_ws
    )

    nrow = 10
    ncol = 10

    perlen = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

    dis = flopy.modflow.ModflowDis(
        ml,
        nlay=1,
        nrow=10,
        ncol=10,
        nper=12,
        delr=63.6,
        delc=63.6,
        top=100,
        perlen=list(perlen),
        nstp=list(perlen),
        steady=False,
        itmuni=4,
        lenuni=2
    )

    bas = flopy.modflow.ModflowBas(
        ml,
        ibound=np.ones((nrow, ncol), dtype=int),
        strt=np.ones((nrow, ncol)) * 95,
    )

    upw = flopy.modflow.ModflowUpw(
        ml,
        ipakcb=52,
        laytyp=1,
        hk=1,
        ss=1e-05,
        sy=0.15,
    )

    if name in ("etdemand_div", "etdemand_sup"):
        # build a SFR package
        reach_data = flopy.modflow.ModflowSfr2.get_empty_reach_data(12)
        for i in range(0, 10):
            outreach = i + 2
            if outreach == 11:
                outreach = 0
            reach_data[i] = (
            10 * i + 7, 0, i, 7, i + 1, 1, 100, 99, 0.02, 0.5, 0.000015, 0.35, 0.08, 5, 8.64, i + 1, outreach
            )
        reach_data[10] = (10 * 2 + 6, 0, 2, 6, 11, 1, 100, 99, 0.02, 0.5, 0.0000015, 0.35, 0.08, 5, 8.64, 11, 0)
        reach_data[11] = (10 * 4 + 6, 0, 4, 6, 12, 1, 100, 99, 0.02, 0.5, 0.0000015, 0.35, 0.08, 5, 8.64, 12, 0)

        segdata = {}
        for per in range(12):
            segment_data = flopy.modflow.ModflowSfr2.get_empty_segment_data(12)
            for i in range(1, 11):
                outseg = i + 1
                if outseg == 11:
                    outseg = 0
                flow = 0
                if i == 1:
                    flow = 50
                segment_data[i - 1] = (
                i, 1, outseg, 0, 0, 0, flow, 0, 0, 0, 0.03, 0, 0, 0, 0, 0, 0.000015,
                0.5, 99 - ((i -1) * 0.02), 5, 0, 0, 0, 0, 0, 0.000015, 0.5,  99 - (i * 0.02), 5, 0, 0, 0, 0, 0)

            divflow = 10
            segment_data[10] = (
                11, 1, 0, 3, 0, 0, divflow * 2, 0, 0, 0, 0.03, 0, 0, 0, 0, 0, 0.000015,
                0.5, 98.94, 5, 0, 0, 0, 0, 0, 0.000015, 0.5, 98.92, 5, 0, 0, 0, 0, 0)
            segment_data[11] = (
                12, 1, 0, 5, 0, 0, divflow, 0, 0, 0, 0.03, 0, 0, 0, 0, 0, 0.000015,
                0.5, 98.90, 5, 0, 0, 0, 0, 0, 0.000015, 0.5, 98.88, 5, 0, 0, 0, 0,
                0)

            segdata[per] = segment_data

        sfr = flopy.modflow.ModflowSfr2(
            ml,
            nstrm=11,
            nss=11,
            const=86400,
            dleak=0.0001,
            ipakcb=52,
            reach_data=reach_data,
            segment_data=segdata,

        )

        gage = flopy.modflow.ModflowGage(
            ml,
            numgage=3,
            gage_data=[(11, 1, 55, 0), (12, 1, 56, 5), (10, 1, 57, 0)],
            files=["gage11.txt", "gage12.txt", f"{name}_gage10.txt"]
        )

    # build a UZF package
    cimis_data = os.path.join("..", "..", "data", "davis_monthly_ppt_eto.txt")
    df = pd.read_csv(cimis_data)

    finf = {}
    pet = {}
    extdp = {}
    extwc = {}
    for i in range(12):
        finf[i] = np.ones((nrow, ncol)) * round_to_n(
            df.ppt_avg_m.values[i] / perlen[i], 5)
        pet[i] = np.ones((nrow, ncol)) * round_to_n(
            df.eto_avg_m.values[i] / perlen[i], 5)
        extdp[i] = np.ones((nrow, ncol)) * 6
        extwc[i] = np.ones((nrow, ncol)) * 0.06

    uzf = flopy.modflow.ModflowUzf1(
        ml,
        nuztop=1,
        iuzfopt=1,
        irunflg=0,
        ipakcb=52,
        ietflg=1,
        ntrail2=7,
        nsets=40,
        surfdep=0.33,
        iuzfbnd=np.ones((nrow, ncol)),
        vks=8.64,
        eps=5,
        thts=0.35,
        thtr=0.05,
        thti=0.08,
        specifythtr=True,
        specifythti=True,
        finf=finf,
        pet=pet,
        extdp=extdp,
        extwc=extwc

    )

    stress_period_data = {}
    for kper, nts in enumerate(perlen):
        for ts in range(nts):
            stress_period_data[(kper, ts)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(
        ml,
        stress_period_data=stress_period_data
    )

    nwt = flopy.modflow.ModflowNwt(
        ml,
        headtol=2e-02,
        fluxtol=50,
        options="COMPLEX"
    )

    if name == "etdemand_sup":
        ag = build_agsupplemental_package(ml)
    elif name == "etdemand_div":
        ag = build_agdiversion_package(ml)
    else:
        ag = build_agwell_package(ml)

    ml.write_input()
    nam = os.path.join(model_ws, f"{name}.nam")

    if name in ("etdemand_div", "etdemand_sup"):
        with open(nam, "a") as foo:
            foo.write(f"DATA              53  {name}.diversion11.txt\n")
            foo.write(f"DATA              54  {name}.diversion12.txt")
    else:
        with open(nam, "a") as foo:
            foo.write(f"DATA              53  {name}.well1.txt\n")
            foo.write(f"DATA              54  {name}.well2.txt")
    return ml


def build_agdiversion_package(ml):

    dis = ml.dis

    options = flopy.utils.OptionBlock(
        f"ETDEMAND 1.0 IRRIGATION_DIVERSION 2 4 "
        "TIMESERIES_DIVERSION".lower(),
        flopy.modflow.ModflowAg
    )

    timeseries_diversion = flopy.modflow.ModflowAg.get_empty(2, block='time series')
    timeseries_diversion[0] = ("diversion", 11, 53)
    timeseries_diversion[1] = ("diversion", 12, 54)

    irrdiversion = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(2, 2, "irrdiversion")
        # todo: update these with segids
        spd[0] = (11, 2, 31, 0.2, 4, 3, 0.0, 0.5, 4, 4, 0.0, 0.5)
        spd[1] = (12, 2, 31, 0.2, 2, 3, 0.0, 0.5, 2, 4, 0.0, 0.5)
        irrdiversion[i] = spd

    ag = flopy.modflow.ModflowAg(
        ml,
        options=options,
        time_series=timeseries_diversion,
        irrdiversion=irrdiversion,
        nper=dis.nper
    )

    return ag


def build_agsupplemental_package(ml):

    dis = ml.dis

    options = flopy.utils.OptionBlock(
        f"ETDEMAND 1.0 IRRIGATION_DIVERSION 2 2 "
        "IRRIGATION_WELL 2 2 "
        "SUPPLEMENTAL_WELL 2 2 "
        "MAXWELLS 2 "
        "TIMESERIES_DIVERSION "
        "WELLCBC 52".lower(),
        flopy.modflow.ModflowAg
    )

    timeseries_diversion = flopy.modflow.ModflowAg.get_empty(2, block='time series')
    timeseries_diversion[0] = ("diversion", 11, 53)
    timeseries_diversion[1] = ("diversion", 12, 54)

    well_list = flopy.modflow.ModflowAg.get_empty(2, block="well")
    x = [[0, 1, 2, -50.], [0, 5, 4, -100.] ]

    for ix, rec in enumerate(well_list):
        well_list[ix] = tuple(x[ix])

    irrwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(2, 2, "irrwell")
        spd[0] = (0, 2, 31, 0.2, 2, 3, 0.0, 0.5, 2, 4, 0.0, 0.5)
        spd[1] = (1, 2, 31, 0.2, 4, 3, 0.0, 0.5, 4, 4, 0.0, 0.5)
        irrwell[i] = spd

    irrdiversion = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(2, 2, "irrdiversion")
        spd[0] = (11, 2, 31, 0.2, 2, 3, 0.0, 0.5, 2, 4, 0.0, 0.5)
        spd[1] = (12, 2, 31, 0.2, 4, 3, 0.0, 0.5, 4, 4, 0.0, 0.5)
        irrdiversion[i] = spd

    supwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(2, 1, "supwell")
        spd[0] = (0, 1, 11, 1, 1)
        spd[1] = (1, 1, 12, 1, 1)
        supwell[i] = spd

    ag = flopy.modflow.ModflowAg(
        ml,
        options=options,
        time_series=timeseries_diversion,
        well_list=well_list,
        irrwell=irrwell,
        irrdiversion=irrdiversion,
        supwell=supwell,
        nper=dis.nper
    )

    return ag


def build_agwell_package(ml):

    dis = ml.dis

    options = flopy.utils.OptionBlock(
        f"ETDEMAND 1.0 IRRIGATION_WELL 2 4 "
        "MAXWELLS 2 WELLCBC 52 TIMESERIES_WELL".lower(),
        flopy.modflow.ModflowAg
    )

    timeseries_well = flopy.modflow.ModflowAg.get_empty(2, block='time series')
    timeseries_well[0] = ("well", 1, 53)
    timeseries_well[1] = ("well", 2, 54)

    well_list = flopy.modflow.ModflowAg.get_empty(2, block="well")
    x = [[0, 5, 4, -100.], [0, 1, 2, -50.],]

    for ix, rec in enumerate(well_list):
        well_list[ix] = tuple(x[ix])

    irrwell = {}
    for i in range(12):
        spd = flopy.modflow.ModflowAg.get_empty(2, 2, "irrwell")
        spd[0] = (0, 2, 31, 0.2, 4, 3, 0.0, 0.5, 4, 4, 0.0, 0.5)
        spd[1] = (1, 2, 31, 0.2, 2, 3, 0.0, 0.5, 2, 4, 0.0, 0.5)
        irrwell[i] = spd

    ag = flopy.modflow.ModflowAg(
        ml,
        options=options,
        time_series=timeseries_well,
        well_list=well_list,
        irrwell=irrwell,
        nper=dis.nper
    )

    return ag


if __name__== "__main__":
    load_existing = False
    if not load_existing:
        ml = build_nwt_test_models("etdemand_well")
        ml2 = build_nwt_test_models("etdemand_div")
        ml3 = build_nwt_test_models("etdemand_sup")
    else:
        model_ws = os.path.join(sws, "..", "..", "data", "nwt_etdemand_gwet_test_problems")
        ml = flopy.modflow.Modflow.load(
            "etdemand_well.nam", model_ws=model_ws, version="mfnwt"
        )
        ml2 = flopy.modflow.Modflow.load(
            "etdemand_div.nam", model_ws=model_ws, version="mfnwt"
        )
        ml3 = flopy.modflow.Modflow.load(
            "etdemand_sup.nam", model_ws=model_ws, version="mfnwt"
        )

    ml.run_model()
    ml2.run_model()
    ml3.run_model()