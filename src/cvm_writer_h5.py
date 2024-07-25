#!/usr/bin/env python

import os
import sys
import getopt
import time
import h5py
import numpy as np
from pyproj import Proj, transform

"""
Create a CVM model in HDF5 format based on input parameter files and given CSV data file(s).

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
    Create a CVM model in HDF5 format based on input parameter files and given CSV data file(s).

    Call Arguments:
      -h, --help: Display this message.
      -m, --meta: [required] A comma-separated list of parameter files based on the parameter template files in the template directory.
      -g, --global: [required] A text file, based on the parameter template files in the template directory, that contains the global parameters for all datasets.
      -d, --data: [required] A comma-separated list of CSV data filenames. The number of files listed must match the number of files listed under the -m option.
      -o, --output: [required] The output filename without extension.
"""
    )


def convert_coordinates(x, y, metadata, xy_to_latlon):
    """Convert a given coordinate from geographic to UTM
    or from UTM to geographic"""
    # Define the UTM to Lat/Lon transformer
    utm_proj = Proj(
        proj="utm",
        zone=metadata["global_attrs"]["utm_zone"],
        ellps=metadata["global_attrs"]["ellipsoid"],
    )
    latlon_proj = Proj(proj="latlong", datum=metadata["global_attrs"]["ellipsoid"])

    # Convert UTM to Latitude and Longitude
    if xy_to_latlon:
        y2, x2 = transform(
            utm_proj,
            latlon_proj,
            x,
            y,
        )
    else:
        x2, y2 = utm_proj(x, y)

    return x2, y2


def copy_attrs(nc_obj, h5_obj):
    """
    Copy attributes from a netCDF4 object to an HDF5 object.
    """
    for attr_name in nc_obj.ncattrs():
        h5_obj.attrs[attr_name] = nc_obj.getncattr(attr_name)


def reshape_data_3d(df, y_unique, x_unique, z_unique, column, x, y, z):
    reshaped_data = np.full((len(z_unique), len(y_unique), len(x_unique)), np.nan)
    for idx, row in df.iterrows():
        z_idx = np.where(z_unique == row[z])[0][0]
        y_idx = np.where(y_unique == row[y])[0][0]
        x_idx = np.where(x_unique == row[x])[0][0]
        reshaped_data[z_idx, y_idx, x_idx] = row[column]
    return reshaped_data


def reshape_data_2d(
    df, y_unique, x_unique, column, x, y, metadata, is_x2=False, is_y2=False
):
    """Reshapes a 2D grid to make sure all rows and columns exist. In case of x2 and y2, NaN for
    the missing data is not desirable. If we have x2 or y2, we need to calculate aux coordinate values.
    """
    reshaped_data = np.full((len(y_unique), len(x_unique)), np.nan)
    for idx, row in df.iterrows():
        y_idx = np.where(y_unique == row[y])[0][0]
        x_idx = np.where(x_unique == row[x])[0][0]
        reshaped_data[y_idx, x_idx] = row[column]

    # Assign values to x2 and y2, they cannot be NaN.
    if is_x2:
        for i in range(reshaped_data.shape[0]):
            for j in range(reshaped_data.shape[1]):
                if np.isnan(reshaped_data[i, j]):
                    x2, _ = convert_coordinates(
                        x_unique[j], y_unique[i], metadata, xy_to_latlon=False
                    )
                    reshaped_data[i, j] = x2
    elif is_y2:
        for i in range(reshaped_data.shape[0]):
            for j in range(reshaped_data.shape[1]):
                if np.isnan(reshaped_data[i, j]):
                    _, y2 = convert_coordinates(
                        x_unique[j], y_unique[i], metadata, xy_to_latlon=False
                    )
                    reshaped_data[i, j] = y2
    return reshaped_data


def create_hdf5(
    output,
    data_dict_2d,
    meta_dict_2d,
    data_dict_3d,
    meta_dict_3d,
    global_metadata,
    global_params,
    metadata,
    verbose,
):
    with h5py.File(f"{output}.h5", "w") as f:
        # Global attributes
        for attr, key_value_pairs in global_metadata.items():
            for key, value in key_value_pairs.items():
                f.attrs[key] = value

        # Create volumes group for 3D data if there is any 3D data
        if data_dict_3d:
            root_group = "MODEL"
            if "groups" in global_params:
                if "MODEL" in global_params["groups"]:
                    root_group = global_params["groups"]["MODEL"]
            volumes_group = f.create_group(root_group)

            # Add 3D data variables to the volumes group
            for group, df in data_dict_3d.items():
                metadata, params, coords, var_dict = meta_dict_3d[group]
                # Get unique values to form the unified grid for this dataset
                y_var = coords["y"]["var"]
                x_var = coords["x"]["var"]
                z_var = coords["z"]["var"]
                if "y2" in coords:
                    y2_var = coords["y2"]["var"]
                    x2_var = coords["x2"]["var"]
                else:
                    y2_var = None
                    x2_var = None

                y_unique = np.unique(df[y_var])
                x_unique = np.unique(df[x_var])
                z_unique = np.unique(df[z_var])
                grid_shape = (len(z_unique), len(y_unique), len(x_unique))

                # Create a subgroup for this dataset
                if not group.strip():
                    dataset_group = volumes_group
                    logger.warning(
                        f"[WARN] Empty subgroup '{group}' under the '{root_group}' group. Will put the dataset under '{root_group}' group. Watch for possible conflicts."
                    )
                else:
                    if group not in volumes_group:
                        dataset_group = volumes_group.create_group(group)
                    else:
                        logger.warning(
                            f"[WARN] duplicate subgroup '{group}' under the '{root_group}' group. Watch for possible conflicts."
                        )
                        dataset_group = volumes_group[group]

                # Create coordinate datasets with actual names
                depth = dataset_group.create_dataset(z_var, data=z_unique)
                latitude = dataset_group.create_dataset(y_var, data=y_unique)
                longitude = dataset_group.create_dataset(x_var, data=x_unique)

                # Add attributes for the coordinate datasets
                for dataset, column in zip(
                    [depth, latitude, longitude],
                    [z_var, y_var, x_var],
                ):
                    if column in var_dict:
                        for attr_name, attr_value in var_dict[column].items():
                            if attr_name not in [
                                "column",
                                "dimensions",
                                "variable",
                            ]:
                                dataset.attrs[attr_name] = attr_value

                # Add the data variables
                for column in df.columns:
                    if column not in [z_var, y_var, x_var, x2_var, y2_var]:
                        reshaped_data = reshape_data_3d(
                            df,
                            y_unique,
                            x_unique,
                            z_unique,
                            column,
                            x_var,
                            y_var,
                            z_var,
                        )
                        data_dataset = dataset_group.create_dataset(
                            column, data=reshaped_data
                        )
                        depth.make_scale("phony_dim_0")
                        latitude.make_scale("phony_dim_1")
                        longitude.make_scale("phony_dim_2")
                        data_dataset.dims[0].attach_scale(depth)
                        data_dataset.dims[1].attach_scale(latitude)
                        data_dataset.dims[2].attach_scale(longitude)

                        for attr_name, attr_value in var_dict[column].items():
                            if attr_name not in [
                                "column",
                                "dimensions",
                                "variable",
                            ]:
                                data_dataset.attrs[attr_name] = attr_value

                # Add auxiliary coordinates
                if x2_var in df.columns and y2_var in df.columns:
                    aux_x_data = reshape_data_2d(
                        df,
                        y_unique,
                        x_unique,
                        x2_var,
                        x_var,
                        y_var,
                        metadata,
                        is_x2=True,
                    )
                    aux_y_data = reshape_data_2d(
                        df,
                        y_unique,
                        x_unique,
                        y2_var,
                        x_var,
                        y_var,
                        metadata,
                        is_y2=True,
                    )
                    aux_x_dataset = dataset_group.create_dataset(
                        x2_var, data=np.ma.masked_invalid(aux_x_data), fillvalue=np.nan
                    )
                    aux_y_dataset = dataset_group.create_dataset(
                        y2_var, data=np.ma.masked_invalid(aux_y_data), fillvalue=np.nan
                    )
                    for aux_dataset, aux_column in zip(
                        [aux_x_dataset, aux_y_dataset], [x2_var, y2_var]
                    ):
                        if aux_column in var_dict:
                            for attr_name, attr_value in var_dict[aux_column].items():
                                if attr_name not in [
                                    "column",
                                    "dimensions",
                                    "variable",
                                ]:
                                    aux_dataset.attrs[attr_name] = attr_value

        # Create surfaces group for 2D data if there is any 2D data
        if data_dict_2d:

            root_group = "SURFACES"
            if "groups" in global_params:
                if "SURFACES" in global_params["groups"]:
                    root_group = global_params["groups"]["SURFACES"]

            surfaces_group = f.create_group(root_group)
            # Add 2D data variables to the root_group group
            for group, df in data_dict_2d.items():
                metadata, params, coords, var_dict = meta_dict_2d[group]
                # Get unique values to form the unified grid for this dataset
                y_var = coords["y"]["var"]
                x_var = coords["x"]["var"]
                if "y2" in coords:
                    y2_var = coords["y2"]["var"]
                    x2_var = coords["x2"]["var"]
                else:
                    y2_var = None
                    x2_var = None

                y_unique = np.unique(df[y_var])
                x_unique = np.unique(df[x_var])
                grid_shape = (len(y_unique), len(x_unique))
                # Create a subgroup for this dataset
                if not group.strip():
                    dataset_group = surfaces_group
                    logger.warning(
                        f"[WARN] Empty subgroup '{group}' under the '{root_group}' group. Will put the dataset under '{root_group}' group. Watch for possible conflicts."
                    )
                else:
                    if group not in surfaces_group:
                        dataset_group = surfaces_group.create_group(group)
                    else:
                        logger.warning(
                            f"[WARN] duplicate subgroup '{group}' under the '{group}' group. Watch for possible conflicts."
                        )
                        dataset_group = surfaces_group[group]

                # Create coordinate datasets with actual names
                latitude = dataset_group.create_dataset(y_var, data=y_unique)
                longitude = dataset_group.create_dataset(x_var, data=x_unique)

                # Add attributes for the coordinate datasets
                for dataset, column in zip(
                    [latitude, longitude],
                    [y_var, x_var],
                ):
                    if column in var_dict:
                        for attr_name, attr_value in var_dict[column].items():
                            if attr_name not in [
                                "column",
                                "dimensions",
                                "variable",
                            ]:
                                dataset.attrs[attr_name] = attr_value

                # Add the data variables
                for column in df.columns:
                    if column not in [y_var, x_var, x2_var, y2_var]:
                        reshaped_data = reshape_data_2d(
                            df, y_unique, x_unique, column, x_var, y_var, metadata
                        )
                        data_dataset = dataset_group.create_dataset(
                            column, data=reshaped_data
                        )
                        latitude.make_scale("phony_dim_0")
                        longitude.make_scale("phony_dim_1")
                        data_dataset.dims[0].attach_scale(latitude)
                        data_dataset.dims[1].attach_scale(longitude)

                        for attr_name, attr_value in var_dict[column].items():
                            if attr_name not in [
                                "column",
                                "dimensions",
                                "variable",
                            ]:
                                data_dataset.attrs[attr_name] = attr_value

                # Add auxiliary coordinates
                if x2_var in df.columns and y2_var in df.columns:
                    aux_x_data = reshape_data_2d(
                        df,
                        y_unique,
                        x_unique,
                        x2_var,
                        x_var,
                        y_var,
                        metadata,
                        is_x2=True,
                    )
                    aux_y_data = reshape_data_2d(
                        df,
                        y_unique,
                        x_unique,
                        y2_var,
                        x_var,
                        y_var,
                        metadata,
                        is_y2=True,
                    )
                    aux_x_dataset = dataset_group.create_dataset(
                        x2_var, data=np.ma.masked_invalid(aux_x_data), fillvalue=np.nan
                    )
                    aux_y_dataset = dataset_group.create_dataset(
                        y2_var, data=np.ma.masked_invalid(aux_y_data), fillvalue=np.nan
                    )
                    for aux_dataset, aux_column in zip(
                        [aux_x_dataset, aux_y_dataset], [x2_var, y2_var]
                    ):
                        if aux_column in var_dict:
                            for attr_name, attr_value in var_dict[aux_column].items():
                                if attr_name not in [
                                    "column",
                                    "dimensions",
                                    "variable",
                                ]:
                                    aux_dataset.attrs[attr_name] = attr_value


def main():
    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv,
            "hm:d:o:g:v",
            ["help", "meta=", "data=", "output=", "global=", "verbose"],
        )
    except getopt.GetoptError as err:
        # Print the error, and help information and exit:
        logger.error(err)
        usage()
        sys.exit(2)

    # Initialize the variables.
    meta_file = None
    output = None
    data_file = None
    global_file = None
    verbose = False
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-m", "--meta"):
            meta_file = a
        elif o in ("-d", "--data"):
            data_file = a
        elif o in ("-g", "--global"):
            global_file = a
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

    # Global metadata file is required.
    if global_file is not None:
        if not os.path.isfile(global_file):
            logger.error(f"[ERR] metadata file '{global_file}' not found!")
            sys.exit(1)
    else:
        usage()
        logger.error(f"[ERR] global metadata file is required.")
        sys.exit(1)

    # Metadata file is required.
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

    # Initialize the timer.
    t0 = time.time()

    logger.info(f"[INFO] Creating output file '{output}.h5'.")

    # Initialize dictionaries to store data
    data_dict_2d = dict()
    meta_dict_2d = dict()
    data_dict_3d = dict()
    meta_dict_3d = dict()

    # Go through the global parameter file.
    global_params = lib.read_model_metadata(global_file)
    # All dataset must be under the same grid_mapping_name.
    if "grid_mapping_name" not in global_params:
        logger.error(
            f"[ERR] global metadata file {global_file} is missing the required 'grid_mapping_name' parameter."
        )
        sys.exit(1)
    grid_mapping_name = global_params["grid_mapping_name"]
    # Parse the  global metadata file.
    metadata_dict, global_metadata = meta_lib.get_h5_metadata(
        global_params, global_meta=True
    )

    # Go through each metadata file and check things first
    groups_list = list()
    params_list = list()

    for file_index, mfile in enumerate(meta_file_list):
        _params = lib.read_model_metadata(meta_file_list[file_index])
        params_list.append(_params)

        # dataset_group is required.
        if "dataset_group" not in params_list[-1]:
            logger.error(
                f"[ERR] metadata file {meta_file_list[file_index]} is missing the required 'dataset_group' parameter."
            )
            sys.exit(1)

        groups_list.append(_params["dataset_group"])

        # Parse the metadata file.
        metadata_v, var_dict, data_variables = meta_lib.get_h5_metadata(_params)
        metadata = lib.merge_dictionaries(metadata_v, [global_metadata])

        # MB Need to add a check to make sure the same x, y, z variable names are used.
        # Read the CSV data to Pandas DataFrame.
        df, _params = meta_lib.read_csv(
            data_file_list[file_index],
            _params,
        )

        coords = meta_lib.get_coords(df, metadata, flat="flat")
        # Add the auxiliary coordinates, if necessary.
        if "x2" in coords and "y2" in coords:
            x2_flat = coords["x2"]["data"]
            y2_flat = coords["y2"]["data"]
            if len(x2_flat) == len(df) and len(y2_flat) == len(df):
                df[coords["x2"]["var"]] = x2_flat
                df[coords["y2"]["var"]] = y2_flat
            else:
                raise ValueError(
                    "The length of auxiliary coordinates does not match the DataFrame index length."
                )

        # 2D or 3D?
        if "z" in coords:
            data_dict_3d[groups_list[-1]] = df
            meta_dict_3d[groups_list[-1]] = [metadata, _params, coords, var_dict]
        else:
            data_dict_2d[groups_list[-1]] = df
            meta_dict_2d[groups_list[-1]] = [metadata, _params, coords, var_dict]
    logger.info(f"[INFO] 3D datasets: {data_dict_3d}\n] 2D datasets: {data_dict_2d}")

    create_hdf5(
        output,
        data_dict_2d,
        meta_dict_2d,
        data_dict_3d,
        meta_dict_3d,
        global_metadata,
        global_params,
        metadata,
        verbose,
    )


if __name__ == "__main__":
    main()
