import os
from math import log10, floor
import platform


def round_to_n(x, n):
    if x == 0:
        return 0
    t = round(x, -int(floor(log10(abs(x))) - (n - 1)))
    return t


def mf6_dev_no_final_check(model_ws, fname):
    contents = []
    with open(os.path.join(model_ws, fname)) as foo:
        for line in foo:
            if "options" in line.lower():
                contents.append(line)
                contents.append("  DEV_NO_FINAL_CHECK\n")
            else:
                contents.append(line)

    with open(os.path.join(model_ws, fname), "w") as foo:
        for line in contents:
            foo.write(line)


def dll_loc():
    cplatform = platform.platform().lower()
    if 'window' in cplatform:
        return os.path.join("..", "bin", "libmf6.dll")
    elif ('mac' in cplatform) and ('arm' in cplatform):
        return os.path.join("..", "bin", "mac_arm", "libmf6.dylib")
    elif ('mac' in cplatform) and ('arm' not in cplatform):
        return os.path.join("..", "bin", "mac_x86_64", "libmf6.dylib")
    else:
        return os.path.join("..", "bin", "libmf6.so")


def nwt_output_path():
    return os.path.join("..", "data", "nwt_etdemand_test_problems")


def nwt_prudic_output_path():
    return os.path.join("..", "data", "nwt_prudic_ag")
