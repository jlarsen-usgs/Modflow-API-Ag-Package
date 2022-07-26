# test importing ModflowAgmvr
import sys
import os
import shutil


def test_import_agmvr():
    from mf6api_agmvr import ModflowAgmvr


def test_setup():
    test_dir = os.path.join(".", "temp")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)

    os.makedirs(test_dir)


if __name__ == "__main__":
    test_setup()
    test_import_agmvr()
