import os
from flopy.mf6.api.create_api_package import create_api_package


def make_modflowagmvr_package():
    """
    Method to make a FloPy compatible ModflowGwfagmvr package
    for handling custom IO within FloPy

    Returns
    -------
        None
    """
    ws = os.path.abspath(os.path.dirname(__file__))
    create_api_package(
        dfn_file=os.path.join(ws, "gwf-agmvr.dfn"),
        script_path=os.path.join(ws, "..")
    )


if __name__ == "__main__":
    make_modflowagmvr_package()