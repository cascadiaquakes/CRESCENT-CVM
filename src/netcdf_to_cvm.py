#!/usr/bin/env python

import os
import sys
import getopt
import time
import numpy as np
import xarray as xr


"""
Input a parameter file and a specified netCDF file, then output an upgraded netCDF file that complies with CVM standards.

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
    Input metadata parameter files and a specified netCDF file, then output an upgraded netCDF file that complies with CVM standards.

    Call arguments:
        -h, --help: this message.
        -m, --meta: [required] a text parameter file based on the parameter template files under the template directory.
        -g, --global: [required] a text file based on the parameter template files under the template directory that contains the global parameters for all datasets.
        -d, --data: [required] a netCDF data filename 
        -o, --output: [required] the output filename without extension.
        -v, --verbose: [default False] turns the verbose mode on no parameter required.
"""
    )


def dataset_to_dataframe_with_dims(xr_dset):
    """Function to convert dataset to dataframe."""
    df = xr_dset.to_dataframe().reset_index()
    dim_order = {var: list(xr_dset[var].dims) for var in xr_dset.variables}
    return df, dim_order


def dataframe_to_dataset_with_dims(df, dim_order):
    """Function to convert dataframe back to dataset."""
    ds = df.set_index(list(dim_order.keys())).to_xarray()
    return ds


def preserve_attributes(original_ds, new_ds):
    """Preserving global and variable attributes."""
    new_ds.attrs = original_ds.attrs
    for var in new_ds.variables:
        if var in original_ds.variables:
            new_ds[var].attrs = original_ds[var].attrs


def detect_initial_dimensionality(xr_dset):
    """Function to detect the initial dimensionality of the dataset."""
    if "z" in xr_dset.dims:
        return 3  # 3D dataset
    else:
        return 2  # 2D dataset


def main():
    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv,
            "hm:g:d:o:v",
            ["help", "meta=", "global=", "data=", "output=", "verbose"],
        )
    except getopt.GetoptError as err:
        # Print the error, and help information and exit:
        logger.error(err)
        usage()
        sys.exit(2)

    # Initialize the variables.
    meta_file = None
    global_meta_file = None
    output = None
    data_file = None
    verbose = False
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-m", "--meta"):
            meta_file = a
        elif o in ("-g", "--global"):
            global_meta_file = a
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

    # Metadata files are required.
    if meta_file is not None:
        meta_file_list = meta_file.strip().split(",")
        for mfile in meta_file_list:
            if not os.path.isfile(mfile):
                logger.error(f"[ERR] metadata file '{mfile}' not found!")
                sys.exit(1)
    else:
        usage()
        logger.error(f"[ERR] metadata file is required.")
        sys.exit(1)

    if global_meta_file is not None:
        global_meta_file_list = global_meta_file.strip().split(",")
        for mfile in global_meta_file_list:
            if not os.path.isfile(mfile):
                logger.error(f"[ERR] global metadata file '{mfile}' not found!")
                sys.exit(1)
    else:
        usage()
        logger.error(f"[ERR] metadata file is required.")
        sys.exit(1)

    # Data file is required.
    if data_file is None:
        usage()
        logger.error(f"[ERR] data file is required.")
        sys.exit(1)
    else:
        data_file_list = data_file.strip().split(",")
        if len(data_file_list) != len(meta_file_list):
            usage()
            logger.error(
                f"[ERR] metadata file count of {len(meta_file_list)} is not the same as the data file count of {len(data_file_list)}."
            )
            sys.exit(1)
        else:
            for dfile in data_file_list:
                if not os.path.isfile(dfile):
                    logger.error(f"[ERR] data file '{dfile}' not found!")
                    sys.exit(1)

    data_file = data_file_list[0]
    meta_file = meta_file_list[0]
    global_meta_file = global_meta_file_list[0]
    # Initialize the timer.
    t0 = time.time()
    # Figure out the input's file type.
    file_type = lib.check_file_type(data_file)

    if file_type["engine"] == "netcdf" and file_type["valid"] == True:
        logger.info(f"[INFO] {data_file} is a netCDF file")

    # Parse the metadata file.
    params = lib.read_model_metadata(global_meta_file)
    params = lib.read_model_metadata(meta_file, params=params)

    if not os.path.isfile(meta_file):
        logger.error(f"[ERR] metadata file '{meta_file}' not found!")
        sys.exit(1)
    else:
        metadata_dict, metadata, var_dict, data_variables = meta_lib.get_metadata(
            params
        )
    nc_format = params["netcdf_format"]
    if nc_format not in writer_prop.netcdf_format:
        usage()
        logger.error(
            f"[ERR] Invalid netcdf option: {nc_format}. Must be one of {list(writer_prop.netcdf_format.keys)}."
        )
        sys.exit(1)
    nc_format = writer_prop.netcdf_format[nc_format]
    logger.info(f"[INFO] nc_format: '{nc_format}'")
    with xr.open_dataset(data_file, engine="netcdf4") as xr_dset:
        if verbose:
            logger.info(f"\n\n[INFO] DataSet from netCDF: '{xr_dset}'\n\n")

        # Detect initial dimensionality of the dataset
        initial_dimensionality = detect_initial_dimensionality(xr_dset)
        if verbose:
            logger.info(f"[INFO] Checking the global attributes using {metadata_dict}")
        # Update/add global attributes based on the metadata file.
        for key in metadata_dict:
            attr = metadata_dict[key]
            if key not in xr_dset.attrs:
                logger.info(
                    f"[INFO] Added the missing global attribute {key}: '{attr}'"
                )
                xr_dset.attrs[key] = attr
            elif xr_dset.attrs[key] != attr:
                logger.info(
                    f"[INFO] Updated global attribute '{key}' from '{xr_dset.attrs[key]}' to '{attr}'"
                )
                xr_dset.attrs[key] = attr

        # Update/add variable attributes based on the metadata file.
        for key in var_dict:
            var = var_dict[key]["variable"]
            if var not in xr_dset:
                logger.warning(f"[WARN] variable {var} not in {xr_dset}")
            else:
                for attr_name, attr_value in var_dict[key].items():
                    if attr_name not in ["column", "dimensions", "variable"]:
                        if (
                            attr_name not in xr_dset[var].attrs
                            and attr_name not in xr_dset[var].encoding
                        ):
                            logger.info(
                                f"[INFO] Added the missing {var} attribute {attr_name}: '{attr_value}'"
                            )
                            xr_dset[var].attrs[attr_name] = attr_value

        # Convert the dataset to DataFrame and preserve dimension order
        df, dim_order = dataset_to_dataframe_with_dims(xr_dset)
        coords = meta_lib.get_coords(df, metadata, flat="2d")

        if verbose:
            logger.info(f"[INFO] Checking the coordinates")

        # Check the shape of primary coordinates
        y_shape = xr_dset[coords["y"]["var"]].shape
        x_shape = xr_dset[coords["x"]["var"]].shape
        z_shape = xr_dset[coords["z"]["var"]].shape if "z" in coords else (1,)

        # Calculate the expected shape for the 2D or 3D grid
        if initial_dimensionality == 3:
            expected_shape = (z_shape[0], y_shape[0], x_shape[0])
        else:
            expected_shape = (y_shape[0], x_shape[0])

        # Add the auxiliary coordinates, if necessary.
        if "x2" in coords and "y2" in coords:
            x2_flat = coords["x2"]["data"]
            y2_flat = coords["y2"]["data"]

            # Calculate the expected shape for the 2D grid
            expected_2d_shape = (y_shape[0], x_shape[0])

            # Reshape the auxiliary coordinates to match the 2D grid
            x2_data = np.array(x2_flat).reshape(expected_2d_shape)
            y2_data = np.array(y2_flat).reshape(expected_2d_shape)

            # Ensure the reshaped auxiliary coordinates match the 2D grid
            if (
                x2_data.shape == expected_2d_shape
                and y2_data.shape == expected_2d_shape
            ):
                # Extract the existing attributes.
                x2_attrs = None
                if coords["x2"]["var"] in xr_dset.coords:
                    x2_attrs = xr_dset.coords[coords["x2"]["var"]].attrs
                    y2_attrs = xr_dset.coords[coords["y2"]["var"]].attrs

                # Adding auxiliary coordinates directly to the dataset
                xr_dset.coords[coords["x2"]["var"]] = (
                    coords["y"]["var"],
                    coords["x"]["var"],
                ), x2_data
                xr_dset.coords[coords["y2"]["var"]] = (
                    coords["y"]["var"],
                    coords["x"]["var"],
                ), y2_data

                # Reassign the existing attributes.
                if x2_attrs is not None:
                    xr_dset.coords[coords["x2"]["var"]].attrs = x2_attrs
                    xr_dset.coords[coords["y2"]["var"]].attrs = y2_attrs

            else:
                raise ValueError(
                    "The shape of auxiliary coordinates does not match the expected dimensions."
                )

        # Reshape data values based on the initial dimensionality
        for var in xr_dset.data_vars:
            if xr_dset[var].ndim == 1:
                if initial_dimensionality == 3:
                    xr_dset[var].values = xr_dset[var].values.reshape(expected_shape)
                else:
                    xr_dset[var].values = xr_dset[var].values.reshape(expected_shape)

        # Preserve attributes from the original dataset
        reshaped_data = xr_dset.copy()
        preserve_attributes(xr_dset, reshaped_data)

        # Add the auxiliary coordinates attributes to the reshaped_data, if the variables exist.
        if "x2" in coords and "y2" in coords:
            y2_var = coords["y2"]["var"]
            x2_var = coords["x2"]["var"]
            for key in var_dict:
                var = var_dict[key]["variable"]
                if var in (x2_var, y2_var):
                    if var not in reshaped_data:
                        logger.warning(f"[WARN] variable {var} not in {reshaped_data}")
                    else:
                        for attr_name, attr_value in var_dict[key].items():
                            if attr_name not in ["column", "dimensions", "variable"]:
                                if (
                                    attr_name not in reshaped_data[var].attrs
                                    and attr_name not in xr_dset[var].encoding
                                ):
                                    logger.info(
                                        f"[INFO] Added the missing {var} attribute {attr_name}: '{attr_value}'"
                                    )
                                    reshaped_data[var].attrs[attr_name] = attr_value

        # Output the netCDF file.
        for _var in reshaped_data.coords:
            reshaped_data.coords[_var].encoding["_FillValue"] = None

        if verbose:
            logger.info(f"\n\n[INFO] Updated dataset: {reshaped_data}\n\n")

        message = meta_lib.write_netcdf_file(output, reshaped_data, nc_format)
        t0, time_txt = lib.time_it(t0)
        logger.info(f"[{time_txt}] {message}")


if __name__ == "__main__":
    main()
