#!/usr/bin/env python

import sys
import os
import getopt

"""Convert a CVM netCDF file to GeoCSV:


Author: Manochehr Bahavar / EarthScope
"""

script = os.path.basename(sys.argv[0])
# Get the directory paths.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop
import netcdf_to_geocsv_prop as convert_prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib
import netcdf_to_geocsv_lib as convert_lib

# Set up the logger.
logger = lib.get_logger()
line_break = "\n"


def usage():
    logger.info(
        f"""
    Convert a CVM netCDF file to GeoCSV

    Call arguments:
        -h, --help: this message.
        -i, --input: [required] the input nefiletCDF filename.
        -o, --output: [optional] the output  GeoCSV filename. If not provided, it will have the same filename as the input file.
"""
    )


def main():
    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(argv, "hi:o:", ["help", "input=", "output="])
    except getopt.GetoptError as err:
        # Print the error, and help information and exit:
        logger.error(err)
        usage()
        sys.exit(2)
    # Initialize the variables.
    input_file = None
    output_file = None
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--input"):
            input_file = a
            if not os.path.isfile(input_file):
                logger.error(
                    f"[ERR] Invalid input netCDF file: [{input_file}]. File not found!"
                )
                sys.exit(2)
        elif o in ("-o", "--output"):
            output_file = a.strip()
        else:
            assert False, "unhandled option"

    # The input filename is required.
    if input_file is None:
        usage()
        logger.error("[ERR] missing -i or --input.")
        sys.exit(1)

    logger.info(f"[INFO] Working on input {input_file}")
    # Set the output filename the same as the input file.
    if output_file is None:
        output_file = f"{os.path.splitext(input_file)[0]}{prop.extension['geocsv']}"
    logger.info(f"[INFO] Will write the output to {output_file}")

    # Figure out the input's file type.
    file_type = lib.check_file_type(input_file)

    # Convert netCDF to GeoCSV.
    if file_type["engine"] == "netcdf" and file_type["valid"] == True:
        logger.info(f"[INFO] {input_file} is a netCDF file")
        metadata, data = convert_lib.netcdf_to_geocsv(input_file)
        with open(output_file, "w") as outfile:
            outfile.write(f"{metadata}\n{data}")
    else:
        logger.error(f"[ERR] {input_file} is not a valid netCDF file")
        sys.exit(1)


if __name__ == "__main__":
    main()
