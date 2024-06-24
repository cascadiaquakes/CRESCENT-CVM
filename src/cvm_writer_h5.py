#!/usr/bin/env python

import os
import sys
import getopt
import time
import h5py
import numpy as np


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


def main():
    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv, "hm:d:o:g:", ["help", "meta=", "data=", "output=", "global="]
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

    # Start a new HDF5 file
    with h5py.File(f"{output}.h5", "w") as f:

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
            df = meta_lib.read_csv(
                data_file_list[file_index],
                _params,
            )

            coords = meta_lib.get_coords(df, metadata, flat="flat")
            # Add the auxiliary coordinates, if necessary. (skip for now)
            if coords["x2"]["var"] not in df and "x2" in coords:
                # Ensure correct shape and update DataFrame with auxiliary coordinates
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
        logger.info(
            f"[INFO] 3D datasets: {data_dict_3d}\n] 2D datasets: {data_dict_2d}"
        )

        # Function to reshape 3D data to 3D grid
        def reshape_data_3d(df, y_unique, x_unique, z_unique, column, x, y, z):
            reshaped_data = np.full(
                (len(y_unique), len(x_unique), len(z_unique)), np.nan
            )
            for idx, row in df.iterrows():
                y_idx = np.where(y_unique == row[y])[0][0]
                x_idx = np.where(x_unique == row[x])[0][0]
                z_idx = np.where(z_unique == row[z])[0][0]
                reshaped_data[y_idx, x_idx, z_idx] = row[column]
            return reshaped_data

        # Function to reshape 2D data to 2D grid
        def reshape_data_2d(df, y_unique, x_unique, column, x, y):
            reshaped_data = np.full((len(y_unique), len(x_unique)), np.nan)
            for idx, row in df.iterrows():
                y_idx = np.where(y_unique == row[y])[0][0]
                x_idx = np.where(x_unique == row[x])[0][0]
                reshaped_data[y_idx, x_idx] = row[column]
            return reshaped_data

        # Global attributes
        for attr, key_value_pairs in global_metadata.items():
            for key, value in key_value_pairs.items():
                f.attrs[key] = value

        # Create volumes group for 3D data if there is any 3D data
        if data_dict_3d:
            root_group = "volumes"
            if "groups" in global_params:
                if "volumes" in global_params["groups"]:
                    root_group = global_params["groups"]["volumes"]
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
                grid_shape = (len(y_unique), len(x_unique), len(z_unique))

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

                # Create dimension datasets
                y_dataset = dataset_group.create_dataset(
                    coords["y"]["var"], data=y_unique
                )
                x_dataset = dataset_group.create_dataset(
                    coords["x"]["var"], data=x_unique
                )
                z_dataset = dataset_group.create_dataset(
                    coords["z"]["var"], data=z_unique
                )

                # Add attributes for the dimension datasets
                for dataset, column in zip(
                    [y_dataset, x_dataset, z_dataset],
                    [y_var, x_var, z_var, x2_var, y2_var],
                ):

                    for attr_name, attr_value in var_dict[column].items():
                        if attr_name not in [
                            "column",
                            "dimensions",
                            "variable",
                        ]:
                            dataset.attrs[attr_name] = attr_value

                # Add the data variables
                for column in df.columns:
                    if column not in [y_var, x_var, z_var, x2_var, y2_var]:
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
                        data_dataset.dims[0].attach_scale(y_dataset)
                        data_dataset.dims[1].attach_scale(x_dataset)
                        data_dataset.dims[2].attach_scale(z_dataset)

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
                        df, y_unique, x_unique, x2_var, x_var, y_var
                    )
                    aux_y_data = reshape_data_2d(
                        df, y_unique, x_unique, y2_var, x_var, y_var
                    )
                    aux_x_dataset = dataset_group.create_dataset(
                        x2_var, data=aux_x_data
                    )
                    aux_y_dataset = dataset_group.create_dataset(
                        y2_var, data=aux_y_data
                    )
                    for aux_dataset, aux_column in zip(
                        [aux_x_dataset, aux_y_dataset], [x2_var, y2_var]
                    ):
                        for attr_name, attr_value in var_dict[aux_column].items():
                            if attr_name not in [
                                "column",
                                "dimensions",
                                "variable",
                            ]:
                                aux_dataset.attrs[attr_name] = attr_value

        # Create surfaces group for 2D data if there is any 2D data
        if data_dict_2d:

            root_group = "surfaces"
            if "groups" in global_params:
                if "surfaces" in global_params["groups"]:
                    root_group = global_params["groups"]["surfaces"]

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

                # Create dimension datasets
                y_dataset = dataset_group.create_dataset(
                    coords["y"]["var"], data=y_unique
                )
                x_dataset = dataset_group.create_dataset(
                    coords["x"]["var"], data=x_unique
                )

                # Add attributes for the dimension datasets
                for dataset, column in zip(
                    [y_dataset, x_dataset],
                    [y_var, x_var, x2_var, y2_var],
                ):

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
                            df,
                            y_unique,
                            x_unique,
                            column,
                            x_var,
                            y_var,
                        )
                        data_dataset = dataset_group.create_dataset(
                            column, data=reshaped_data
                        )
                        data_dataset.dims[0].attach_scale(y_dataset)
                        data_dataset.dims[1].attach_scale(x_dataset)

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
                        df, y_unique, x_unique, x2_var, x_var, y_var
                    )
                    aux_y_data = reshape_data_2d(
                        df, y_unique, x_unique, y2_var, x_var, y_var
                    )
                    aux_x_dataset = dataset_group.create_dataset(
                        x2_var, data=aux_x_data
                    )
                    aux_y_dataset = dataset_group.create_dataset(
                        y2_var, data=aux_y_data
                    )
                    for aux_dataset, aux_column in zip(
                        [aux_x_dataset, aux_y_dataset], [x2_var, y2_var]
                    ):
                        for attr_name, attr_value in var_dict[aux_column].items():
                            if attr_name not in [
                                "column",
                                "dimensions",
                                "variable",
                            ]:
                                aux_dataset.attrs[attr_name] = attr_value


if __name__ == "__main__":
    main()
