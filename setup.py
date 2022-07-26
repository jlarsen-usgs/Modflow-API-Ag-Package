import setuptools


if __name__ == "__main__":
    setuptools.setup()

    import os
    import flopy
    import mf6api_agmvr
    import shutil
    import subprocess
    flopy_dir = os.path.split(flopy.__file__)[0]
    ag_dir = os.path.split(mf6api_agmvr.__file__)[0]

    flopy_ag_dfn = os.path.join(flopy_dir, "mf6", "data", "dfn", "gwf-agmvr.dfn")
    ag_dfn = os.path.join(ag_dir, "dfn", "gwf-agmvr.dfn")

    shutil.copy(ag_dfn, flopy_ag_dfn)

    flopy.mf6.utils.createpackages.create_packages()
    print('break')
