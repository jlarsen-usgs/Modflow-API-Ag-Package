import setuptools
import sys
import os
# from flopy.mf6.api.create_api_package import create_api_package


if __name__ == "__main__":
    # code for future plugin adaptation
    # sys.path.append(os.path.join("mf6api_agmvr", "dfn"))
    # from create_flopy_package import make_modflowagmvr_package

    # make_modflowagmvr_package()

    # ws = os.path.abspath(os.path.dirname(__file__))
    # ag_dfn = os.path.join(ws, "mf6api_agmvr", "dfn", "gwf-agmvr.dfn")
    # script_loc = os.path.join(ws, "mf6api_agmvr")
    # create_api_package(ag_dfn, script_loc)
    
    setuptools.setup()


    
    # current version of the setup, however this will be changed for plugin
    # compatibility
    import flopy
    import mf6api_agmvr
    import shutil
    flopy_dir = os.path.split(flopy.__file__)[0]
    ag_dir = os.path.split(mf6api_agmvr.__file__)[0]

    flopy_ag_dfn = os.path.join(flopy_dir, "mf6", "data", "dfn", "gwf-agmvr.dfn")
    ag_dfn = os.path.join(ag_dir, "dfn", "gwf-agmvr.dfn")

    shutil.copy(ag_dfn, flopy_ag_dfn)

    flopy.mf6.utils.createpackages.create_packages()
    
