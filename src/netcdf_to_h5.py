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
    Convert a CVM netCDF file to CVM HDF5 format.

    Call arguments:
        -h, --help: this message.
        -i, --input: [required] a CVM netCDF model filename.
        -o, --output: [required] the HDF5 output filename without extension.
        -g, --group: [default: MODEL] The main HDF5 group to store the netCDF data under.
        -s, --subgroup [default: ] The subgroup under the above group to store the netCDF data under.
        -v, --verbose: [default False] turns the verbose mode on no parameter required.
"""
    )


def copy_attrs(nc_obj, h5_obj):
    """
    Copy attributes from a netCDF4 object to an HDF5 object.
    """
    for attr_name in nc_obj.ncattrs():
        h5_obj.attrs[attr_name] = nc_obj.getncattr(attr_name)


def convert_netcdf_to_hdf5(netcdf_file, hdf5_file, group_name, subgroup_name):
    # Open the NetCDF file
    with nc.Dataset(netcdf_file, "r") as nc_data:
        # Open the HDF5 file
        with h5py.File(hdf5_file, "w") as hdf:
            # Create the main group
            main_group = hdf.require_group(group_name)

            # Create or access the specified subgroup within the group
            if len(subgroup_name) > 0:
                subgroup = main_group.require_group(subgroup_name)
            else:
                subgroup = main_group

            group_subgroup = subgroup.name

            # Copy dimensions from NetCDF to HDF5, including actual values
            dim_mapping = {}
            for dim_name, dim in nc_data.dimensions.items():
                if dim_name not in subgroup:
                    # Create dataset for the dimension using actual values
                    if (
                        dim_name in nc_data.variables
                    ):  # Check if dimension has a corresponding variable
                        dim_data = nc_data.variables[dim_name][:]
                    else:
                        dim_data = np.arange(
                            dim.size
                        )  # Fallback to indices if no variable is found

                    dim_dataset = subgroup.create_dataset(dim_name, data=dim_data)
                    dim_dataset.attrs["units"] = (
                        nc_data.variables[dim_name].units
                        if dim_name in nc_data.variables
                        and hasattr(nc_data.variables[dim_name], "units")
                        else ""
                    )
                    dim_dataset.attrs["long_name"] = (
                        nc_data.variables[dim_name].long_name
                        if dim_name in nc_data.variables
                        and hasattr(nc_data.variables[dim_name], "long_name")
                        else ""
                    )
                    dim_mapping[dim_name] = dim_dataset
                else:
                    logger.warning(
                        f"[WARN] Dimension '{dim_name}' already exists in the subgroup '{subgroup_name}'. Skipping."
                    )
                    dim_mapping[dim_name] = subgroup[dim_name]

            # Copy variables and their attributes, respecting dimension order
            for var_name, var in nc_data.variables.items():
                if var_name not in subgroup:
                    data = var[:]
                    dim_names = var.dimensions
                    hdf_var = subgroup.create_dataset(
                        var_name, data=data, compression="gzip", shape=data.shape
                    )

                    # Copy variable attributes
                    for attr_name in var.ncattrs():
                        hdf_var.attrs[attr_name] = var.getncattr(attr_name)

                    # Link the dimensions to the dataset
                    for i, dim_name in enumerate(dim_names):
                        if i == 0:
                            # Create the scale for the first dimension
                            subgroup[dim_name].make_scale(dim_name)
                        hdf_var.dims[i].attach_scale(subgroup[dim_name])
                else:
                    logger.warning(
                        f"[WARN] Variable '{var_name}' already exists in the subgroup '{subgroup_name}'. Skipping."
                    )

            # Copy global attributes from NetCDF to HDF5
            for attr_name in nc_data.ncattrs():
                if attr_name not in subgroup.attrs:
                    subgroup.attrs[attr_name] = nc_data.getncattr(attr_name)
                else:
                    logger.warning(
                        f"[WARN] Global attribute '{attr_name}' already exists in the subgroup '{subgroup_name}'. Skipping."
                    )

    logger.info(
        f"[INFO] Conversion from {netcdf_file} to {hdf5_file} under '{group_subgroup}' completed successfully."
    )


if __name__ == "__main__":
    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv,
            "hi:g:,s:o:v",
            ["help", "input=", "group=", "subgroup=", "output=", "verbose"],
        )
    except getopt.GetoptError as err:
        # Print the error, and help information and exit:
        logger.error(err)
        usage()
        sys.exit(2)

    # Initialize the variables.
    output = None
    input = None
    verbose = False
    group = "MODEL"
    subgroup = ""
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--input"):
            input = a
        elif o in ("-g", "--group"):
            group = a
        elif o in ("-s", "--subgroup"):
            subgroup = a
        elif o in ("-o", "--output"):
            output = a
        elif o in ("-v", "--verbose"):
            verbose = True
        else:
            assert False, "unhandled option"

    # The input file is required.
    if input is None:
        usage()
        logger.error(f"[ERR] data file is required.")
        sys.exit(1)
    else:
        netcdf_file = input

    # Figure out the input's file type.
    file_type = lib.check_file_type(netcdf_file)

    if file_type["engine"] == "netcdf" and file_type["valid"] == True:
        logger.info(f"[INFO] {netcdf_file} is a netCDF file")
    else:
        usage()
        logger.error(
            f"[ERR] bad input file type of {file_type['engine']} for {netcdf_file}."
        )
        sys.exit(1)

    # The output file is required.
    if output is None:
        usage()
        logger.error(f"[ERR] output file is required.")
        sys.exit(1)
    hdf5_file = output
    # Add the h5 extension if missing.
    if not hdf5_file.endswith(".h5"):
        hdf5_file = f"{output}.h5"

    # Initialize the timer.
    t0 = time.time()

    convert_netcdf_to_hdf5(netcdf_file, hdf5_file, group, subgroup)
