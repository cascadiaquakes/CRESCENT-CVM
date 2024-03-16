#!/usr/bin/env python

import os
import sys
import getopt
import importlib
import time

"""
Output metadata and model files based on an input parameter file and a given data file. The output can be one more of the following:
    -- Metadata in GeoCSV
    -- Metadata in JSON
    -- Model in GeoCSV
    -- Model in netCDF

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
    Output metadata and model based on an input parameter file and a given data file. The output can be one or more of the following:
        -- Metadata in GeoCSV
        -- Metadata in JSON
        -- Model in GeoCSV
        -- Model in netCDF

    Call arguments:
        -h, --help: this message.
        -m, --meta: [required] a Python parameter file based on the parameter template files under the template directory.
                        The templates files are fully documented. The filename should have the ".py" extension and should not include "." 
                        in the filename.
        -d, --data: [required] a CSV data filename with a header
        -o, --output: [required] the output filename without extension. The metadata output files will add "_metadata" to this output filename.
        -t, --output_types: [default netcdf] a comma separated string of output types. The valid output types are: "metadata", "geocsv", "netcdf"
             metadata: output metadata as GeoCSV and JSON files.
             geocsv: output the model in GeoCSV format.
             netcdf: output the model in netCDF format.
"""
    )


def main():
    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv, "hm:d:o:t:", ["help", "meta=", "data=", "output=", "output_types="]
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
    output_types = "netcdf"
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-m", "--meta"):
            meta_file = a
            if not os.path.isfile(meta_file):
                logger.error(
                    f"[ERR] Invalid metadata file: {meta_file}. File not found!"
                )
                sys.exit(2)
            params = lib.read_model_metadata(meta_file)
        elif o in ("-d", "--data"):
            data_file = a
        elif o in ("-o", "--output"):
            output = a
        elif o in ("-t", "--output_types"):
            output_types = a.strip()
        else:
            assert False, "unhandled option"

    # The output file is required.
    if output is None:
        usage()
        logger.error(f"[ERR] output file is required.")
        sys.exit(1)

    # Metadata file is required.
    if meta_file is None:
        usage()
        logger.error(f"[ERR] metadata file is required.")
        sys.exit(1)

    # Data file is required.
    if data_file is None:
        usage()
        logger.error(f"[ERR] data file is required.")
        sys.exit(1)

    # Check the output types.
    output_types = output_types.split(",")
    for _type in output_types:
        if _type not in writer_prop.valid_output_types:
            usage()
            logger.error(
                f"[ERR] Invalid output type [{_type}] requested. Must be one of {list(writer_prop.valid_output_types)}."
            )
            sys.exit(1)
    if "netcdf" in output_types:
        nc_format = params["netcdf_format"]
        if nc_format not in writer_prop.netcdf_format:
            usage()
            logger.error(
                f"[ERR] Invalid netcdf option: {nc_format}. Must be one of {list(writer_prop.netcdf_format.keys)}."
            )
            sys.exit(1)
        nc_format = writer_prop.netcdf_format[nc_format]

    # Parse the metadata file.
    if not os.path.isfile(meta_file):
        logger.error(f"[ERR] metadata file '{meta_file}' not found!")
        sys.exit(1)
    else:
        metadata_dict, metadata, var_dict, data_variables = meta_lib.get_metadata(
            params
        )

    # Initialize the timer.
    t0 = time.time()

    # Metadata outputs, if requested.
    if "metadata" in output_types:
        # Write metadata as JSON.
        message = meta_lib.write_json_metadata(
            f"{output}_metadata", metadata_dict, var_dict, data_variables
        )
        t0, time_txt = lib.time_it(t0)
        logger.info(f"[{time_txt}] {message}")

        # Write metadata as GeoCSV.
        message = meta_lib.write_geocsv_metadata(
            params, f"{output}_metadata", metadata_dict, var_dict
        )
        t0, time_txt = lib.time_it(t0)
        logger.info(f"[{time_txt}] {message}")

    # NetCDF or geocsv model files requested?
    if "netcdf" in output_types or "geocsv" in output_types:
        # Proceed only if a data file is given
        if data_file is None:
            logger.warning(f"[WARN] No data file was provided.")
            sys.exit(0)
        elif not os.path.isfile(data_file):
            logger.error(f"[ERR] data file '{data_file}' not found!")
            sys.exit(1)
        else:
            df = meta_lib.read_csv(
                data_file,
                params,
            )

            t0, time_txt = lib.time_it(t0)
            logger.info(f"[{time_txt}] data in DataFrame")

            coords = meta_lib.get_coords(df, metadata)
            t0, time_txt = lib.time_it(t0)
            logger.info(f"[{time_txt}] Got the coordinates")

            # Model output as netCDF.
            if "netcdf" in output_types:
                # Initialize a 2D or 3D grid.
                grid = meta_lib.init_grid(coords)
                t0, time_txt = lib.time_it(t0)
                logger.info(f"[{time_txt}] Initialized the grid")

                # Create a grid for each variable.
                var_grid = meta_lib.create_var_grid(df, data_variables, coords, grid)

                t0, time_txt = lib.time_it(t0)
                logger.info(f"[{time_txt}] Created a grid for each variable")

                # Create the Xarray DataFrame.
                xr_df = meta_lib.build_xarray_dataframe(
                    data_variables, coords, var_grid, metadata
                )
                t0, time_txt = lib.time_it(t0)
                logger.info(f"[{time_txt}] Created the DataFrame")

                # Convert the xr DataFrame to an xr Dataset
                xr_dset = meta_lib.build_xarray_dataset(xr_df, coords, metadata)
                t0, time_txt = lib.time_it(t0)
                logger.info(f"[{time_txt}] Converted the DataFrame to Dataset")

                # Output the netCDF file.
                # For the coordinate variables, CF convention requires excluding the _FillValue.
                # So, we have to  explicitly disable it via encoding.
                for _var in xr_dset.coords:
                    xr_dset.coords[_var].encoding["_FillValue"] = None
                message = meta_lib.write_netcdf_file(output, xr_dset, nc_format)
                t0, time_txt = lib.time_it(t0)
                logger.info(f"[{time_txt}] {message}")

            # Model as GeoCSV.
            if "geocsv" in output_types:
                df = meta_lib.add_aux_coord_columns(df, coords)

                # Write metadata as GeoCSV.
                message = meta_lib.write_geocsv_file(
                    df, params, output, metadata_dict, var_dict
                )
                t0, time_txt = lib.time_it(t0)
                logger.info(f"[{time_txt}] {message}")


if __name__ == "__main__":
    main()
