import os
from typing import Union

import flopy
import numpy as np
import pandas as pd


def get_track_dtype(path: os.PathLike):
    """Get the dtype of the track data recarray from the ascii header file."""

    hdr_lns = open(path).readlines()
    hdr_lns_spl = [[ll.strip() for ll in l.split(",")] for l in hdr_lns]
    return np.dtype(list(zip(hdr_lns_spl[0], hdr_lns_spl[1])))


def check_track_data(
    track_bin: os.PathLike,
    track_hdr: os.PathLike,
    track_csv: os.PathLike,
):
    """Check that track data written to binary, CSV, and budget files are equal."""

    # get dtype from ascii header file
    dt = get_track_dtype(track_hdr)

    # read output files
    data_bin = np.fromfile(track_bin, dtype=dt)
    data_csv = np.genfromtxt(track_csv, dtype=dt, delimiter=",", names=True)
    if len(data_csv.shape) == 0:
        # https://stackoverflow.com/a/24943993/6514033
        data_csv = np.array([data_csv])

    assert (
        data_bin.shape == data_csv.shape
    ), f"Binary and CSV track data shapes do not match: {data_bin.shape} != {data_csv.shape}"

    # check particle tracks written to all output files are equal
    # check each column separately to avoid:
    # TypeError: The DType <class 'numpy._FloatAbstractDType'> could not be promoted by <class 'numpy.dtype[void]'>
    for k in data_bin.dtype.names:
        assert np.allclose(data_bin[k], data_csv[k], equal_nan=True)

    # make sure columns all have values in the expected range
    assert all(data_bin["iprp"] >= 1)
    assert all(data_bin["irpt"] >= 1)
    assert all(data_bin["kper"] >= 1)
    assert all(data_bin["kstp"] >= 1)
    assert all(data_bin["ilay"] >= 1)
    assert all(data_bin["icell"] >= 1)
    assert all(data_bin["istatus"] >= 0)
    assert all(data_bin["ireason"] >= 0)


def to_mp7_format(data: Union[pd.DataFrame, np.recarray]) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        data = data.to_records(index=False)

    mp7_dtypes = np.dtype(
        [
            ("particleid", np.int32),
            ("particlegroup", np.int32),
            ("sequencenumber", np.int32),
            ("particleidloc", np.int32),
            ("time", np.float32),
            ("x", np.float32),
            ("y", np.float32),
            ("z", np.float32),
            ("k", np.int32),
            ("node", np.int32),
            ("xloc", np.float32),
            ("yloc", np.float32),
            ("zloc", np.float32),
            ("stressperiod", np.int32),
            ("timestep", np.int32),
        ]
    )

    return pd.DataFrame(
        np.core.records.fromarrays(
            [
                data["irpt"],
                data["iprp"],
                np.zeros(
                    data.shape[0]
                ),  # todo use sequencenumber passed explicitly as particle name
                np.zeros(data.shape[0]),
                data["t"],
                data["x"],
                data["y"],
                data["z"],
                data["ilay"],  # todo add to PRT output?
                data["icell"],
                np.zeros(data.shape[0]),
                np.zeros(data.shape[0]),
                np.zeros(data.shape[0]),
                data["kper"],
                data["kstp"],
            ],
            dtype=mp7_dtypes,
        )
    )
