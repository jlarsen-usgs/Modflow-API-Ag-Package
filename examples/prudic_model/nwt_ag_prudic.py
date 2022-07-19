import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import flopy
import pandas as pd


def build_model(name, model_ws):
    ml = flopy.modflow.Modflow(
        name, version="mfnwt", exe_name="mfnwt", model_ws=model_ws
    )
    data_pth = os.path.join(".", "arrays")

    # build NWT solver
    nwt = flopy.modflow.ModflowNwt(
        ml,
        headtol=1e-04,
        fluxtol=1.0,
        maxiterout=500,
        thickfact=1.0e-07,
        linmeth=2,
        iprnwt=1,
        ibotav=1,
        options="SPECIFIED",
        Continue=True,
        dbdtheta=0.97,
        dbdkappa=1e-05,
        dbdgamma=0.0,
        momfact=0.1,
        backflag=0,
        maxbackiter=30,
        backtol=1.05,
        backreduce=0.9,
        iacl=1,
        norder=0,
        level=7,
        north=7,
        iredsys=1,
        rrctols=0,
        idroptol=1,
        epsrn=5.0e-03,
        hclosexmd=5e-03,
        mxiterxmd=50
    )

    # build DIS package
    nlay, nrow, ncol = 1, 15, 10
    nper = 49
    itmuni, lenuni = 4, 1
    top = np.loadtxt(os.path.join(data_pth, "top.txt"), dtype=float)
    botm = np.loadtxt(os.path.join(data_pth, "bottom.txt"), dtype=float)
    perlen, nstp, tsmult, steady = list(), list(), list(), list()
    for per in range(nper):
        if per == 0:
            steady.append(True)
            perlen.append(1)
            nstp.append(1)
            tsmult.append(1)
        else:
            steady.append(False)
            perlen.append(30.5)
            nstp.append(30)
            tsmult.append(1)

    dis = flopy.modflow.ModflowDis(
        ml,
        nlay,
        nrow,
        ncol,
        nper,
        delr=5000,
        delc=5000,
        top=top,
        botm=botm,
        perlen=perlen,
        nstp=nstp,
        tsmult=tsmult,
        steady=steady,
        itmuni=itmuni,
        lenuni=lenuni
    )

    # build BAS package
    ibound = np.loadtxt(os.path.join(data_pth, "idomain.txt"), dtype=int)
    strt = np.loadtxt(os.path.join(data_pth, "strt.txt"), dtype=float)

    bas = flopy.modflow.ModflowBas(
        ml,
        ibound=ibound,
        strt=strt,
    )

    # build UPW package
    hk = np.loadtxt(os.path.join(data_pth, "hk.txt"), dtype=float)
    sy = np.loadtxt(os.path.join(data_pth, "sy.txt"), dtype=float)

    upw = flopy.modflow.ModflowUpw(
        ml,
        laytyp=1,
        ipakcb=102,
        hk=hk,
        vka=1e-06,
        ss=1e-06,
        sy=sy
    )

    # build GHB package
    stress_period_data = {}
    for per in range(nper):
        rec = [
            (0, 12, 0, 988., 0.038),
            (0, 13, 8, 1045., 0.038)
        ]
        stress_period_data[per] = rec

    ghb = flopy.modflow.ModflowGhb(
        ml,
        ipakcb=102,
        stress_period_data=stress_period_data,
    )

    # build UZF package
    irunbnd = np.loadtxt(os.path.join(data_pth, "irunbnd.txt"), dtype=int)
    finf = np.loadtxt(os.path.join(data_pth, "finf.txt"), dtype=float)
    pet = np.loadtxt(os.path.join(data_pth, "pet.txt"), dtype=float)
    pet = {i: x for i, x in enumerate(pet)}
    ntrail = 15
    nsets = 100
    surfdep = 1.0
    uzfbnd = np.where(ibound != 0, 1, 0)

    options = flopy.utils.OptionBlock("", flopy.modflow.ModflowUzf1)
    options.specifythti = True
    options.specifythtr = True
    options.rejectsurfk = True
    options.specifysurfk = True

    uzf = flopy.modflow.ModflowUzf1(
        ml,
        nuztop=1,
        iuzfopt=1,
        irunflg=1,
        ietflg=1,
        ipakcb=102,
        ntrail2=ntrail,
        nsets=nsets,
        surfdep=surfdep,
        iuzfbnd=uzfbnd,
        irunbnd=irunbnd,
        vks=4.6e-05,
        surfk=1e-05,
        eps=7.5,
        thts=0.38,
        thtr=0.2,
        thti=0.2,
        finf={i: finf for i in range(nper)},
        pet=pet,
        extdp=0.5,
        extwc=0.2,
    )

    # build SFR package
    reach_data = pd.read_csv(
        os.path.join(data_pth, "reach_data.txt"), delim_whitespace=True
    )
    nstrm = len(reach_data)
    nss = len(reach_data.seg.unique())
    const = 86400
    dleak = 1e-06

    reach_data = reach_data.to_list()
    # todo: segment_data SFR
    print('break')
    # todo: AG package

if __name__ == "__main__":
    model_ws = os.path.join("..", "..", "data", "nwt_prudic_ag")
    name = "prudic_ag"
    build_model(name, model_ws)


