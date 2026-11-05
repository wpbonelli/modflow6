"""Compile PRMS as CI does, run nhm.control on drb_2yr, make the inflow nc files.

Args: domain_dir
"""

import pathlib as pl
import subprocess
import sys

import numpy as np
import xarray as xr

import pywatershed as pws

domain = pl.Path(sys.argv[1]).resolve()
exe = pws.utils.compile_prms("5.2.1")
print("prms exe", exe, flush=True)
(domain / "output").mkdir(exist_ok=True)
r = subprocess.run(
    [str(exe), "nhm.control", "-MAXDATALNLEN", "60000"],
    cwd=domain,
    capture_output=True,
    text=True,
)
print(r.stdout[-300:], r.stderr[-300:], flush=True)
r.check_returncode()
for vv in ["sroff", "ssres_flow", "gwres_flow"]:
    out = domain / "output"
    pws.utils.CsvFile(out / f"{vv}.csv").to_netcdf(out / f"{vv}.nc")
    arr = xr.load_dataset(out / f"{vv}.nc")[vv].values
    print(f"INFLOW {vv} mean={np.nanmean(arr):.17e} shape={arr.shape}")
