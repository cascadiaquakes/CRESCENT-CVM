#!/usr/bin/env python

import os
import sys
import getopt
import time
import h5py
import numpy as np
import netCDF4 as nc

"""
Input a parameter file and a CVM-compatible netCDF file, then output an HDF file.

Author: Manochehr Bahavar / EarthScope
"""
# Get the directory paths.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop
import writer_prop as writer_prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib
import writer_lib as meta_lib

# Set up the logger.
logger = lib.get_logger()
line_break = "\n"


def usage():
    logger.info(
        f"""
    Input  a netCDF file to generate an output as HDF5 file.

    Call arguments:
        -h, --help: this message.
        -d, --data: [required] a netCDF data filename 
        -o, --output: [required] the output filename without extension.
        -v, --verbose: [default False] turns the verbose mode on no parameter required.
"""
    )


def copy_attrs(nc_obj, h5_obj):
    """
    Copy attributes from a netCDF4 object to an HDF5 object.
    """
    for attr_name in nc_obj.ncattrs():
        h5_obj.attrs[attr_name] = nc_obj.getncattr(attr_name)


def netcdf_to_hdf5(netcdf_file, hdf5_file):
    # Open the netCDF file
    with nc.Dataset(netcdf_file, "r") as src:
        # Create a new HDF5 file
        with h5py.File(hdf5_file, "w") as dst:
            # Copy global attributes
            copy_attrs(src, dst)

            # Create the MODEL group
            model_group = dst.create_group("MODEL")

            # Copy variables and data from netCDF to HDF5 under MODEL group
            for var_name, var_data in src.variables.items():
                if (
                    var_name in src.dimensions
                ):  # skip dimensions, we'll handle them separately
                    continue
                data = var_data[:]
                # Create the dataset in HDF5 with the same shape and data type as the NetCDF variable
                model_var = model_group.create_dataset(
                    var_name, data=data, chunks=True, compression="gzip"
                )
                copy_attrs(var_data, model_var)
                # Ensure the dimension order is maintained
                model_var.dims.create_scale(model_var, var_name)
                for i, dim in enumerate(var_data.dimensions):
                    model_var.dims[i].attach_scale(model_var)
                if verbose:
                    logger.info(
                        f"[INFO] Copied variable '{var_name}' with shape {data.shape}"
                    )

                # Verification: Read back the data from HDF5 and compare with original
                read_back_data = model_var[:]
                if np.array_equal(data, read_back_data):
                    logger.info(
                        f"[INFO] Data for variable '{var_name}' matches between netCDF and HDF5."
                    )
                else:
                    logger.warning(
                        f"[WARN] Data mismatch for variable '{var_name}' between netCDF and HDF5."
                    )
                    max_difference = np.max(np.abs(data - read_back_data))
                    logger.warning(f"[WARN] Max difference: {max_difference}")

            # Copy dimensions
            for dim_name in src.dimensions:
                dim_data = src.variables[dim_name][:]
                dim_dataset = model_group.create_dataset(dim_name, data=dim_data)
                copy_attrs(src.variables[dim_name], dim_dataset)
                if verbose:
                    logger.info(
                        f"[INFO] Copied dimension '{dim_name}' with data: {dim_data[:5]}..."
                    )

            # Create the surface elevation group
            dst.create_group("SURFACES")

            logger.info(f"[INFO] Converted {netcdf_file} to {hdf5_file}")


if __name__ == "__main__":
    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv,
            "hd:o:v",
            ["help", "data=", "output=", "verbose"],
        )
    except getopt.GetoptError as err:
        # Print the error, and help information and exit:
        logger.error(err)
        usage()
        sys.exit(2)

    # Initialize the variables.
    output = None
    data_file = None
    verbose = False
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-d", "--data"):
            data_file = a
        elif o in ("-o", "--output"):
            output = a
        elif o in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # The output file is required.
    if output is None:
        usage()
        logger.error(f"[ERR] output file is required.")
        sys.exit(1)
    hdf5_file = output
    if not hdf5_file.endswith(".h5"):
        hdf5_file = f"{output}.h5"

    # Data file is required.
    if data_file is None:
        usage()
        logger.error(f"[ERR] data file is required.")
        sys.exit(1)
    else:
        data_file_list = data_file.strip().split(",")
    netcdf_file = data_file_list[0]

    # Initialize the timer.
    t0 = time.time()

    # Figure out the input's file type.
    file_type = lib.check_file_type(data_file)

    if file_type["engine"] == "netcdf" and file_type["valid"] == True:
        logger.info(f"[INFO] {data_file} is a netCDF file")

    netcdf_to_hdf5(netcdf_file, hdf5_file)
