#!/usr/bin/env python3
"""P3.2 acceptance: a frehg2 HDF5 output file opens in h5py with all datasets/metadata
accessible. Usage: check_h5py.py <file.h5> <nx> <ny>."""
import sys

import h5py


def main() -> int:
    path, nx, ny = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    with h5py.File(path, "r") as f:
        meta = f["/simulation/metadata"]
        for attr in ("git_sha", "config_sha256", "io_mode", "mpi_ranks"):
            assert attr in meta.attrs, f"missing metadata attr {attr}"
        dom = f["/domain"]
        assert int(dom.attrs["nx"]) == nx, "domain nx mismatch"
        assert int(dom.attrs["ny"]) == ny, "domain ny mismatch"
        ds = f["/surface/water_depth/0"]
        assert ds.shape[0] == nx * ny, f"surface size {ds.shape[0]} != {nx*ny}"
        assert ds.attrs["units"].decode() if isinstance(ds.attrs["units"], bytes) \
            else ds.attrs["units"] == "m"
        _ = ds[:]  # readable
    print(f"h5py OK: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
