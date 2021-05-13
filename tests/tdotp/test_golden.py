#!/usr/bin/env python3

import argparse
import numpy as np
import xarray as xr
import subprocess


def run_scene(scene_executable: str, input_file: str):
    """Run SCENE on a given input file"""
    if input_file.endswith(".dat"):
        input_file = input_file.replace(".dat", "")
    command = f"echo {input_file} | {scene_executable}"
    subprocess.run(command, shell=True, check=True, encoding="utf-8")


def check_golden(new_filename: str, golden_filename: str):
    """Check two files match using numpy.allclose"""
    new_df = xr.open_dataset(new_filename)
    golden_df = xr.open_dataset(golden_filename)

    failures = {}
    for name in golden_df:
        if name not in new_df:
            failures[name] = "not in new file"
            continue
        if not np.allclose(golden_df[name], new_df[name]):
            error_l_inf = np.max(np.abs(golden_df[name] - new_df[name]))
            failures[name] = f"difference ({error_l_inf}) larger than tolerance"

    return failures


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Check netCDF files match within some tolerance"
    )
    parser.add_argument(
        "-s", "--scene-executable", type=str, help="Path to SCENE executable"
    )
    parser.add_argument(
        "-i", "--input", type=str, help="SCENE input file", dest="input_file"
    )
    parser.add_argument("-n", "--new", type=str, help="New file")
    parser.add_argument("-o", "--old", type=str, help="Old file (golden answer)")

    args = parser.parse_args()

    run_scene(args.scene_executable, args.input_file)
    failures = check_golden(args.new, args.old)

    if failures:
        print(f"Found {len(failures)} errors!")
        for name, reason in failures.items():
            print(f"{name}: {reason}")
        exit(len(failures))

    print("Files match")
