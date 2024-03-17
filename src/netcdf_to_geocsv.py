#!/usr/bin/env python

import sys
import os
import getopt
import json

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
    Convert a CVM netCDF file to GeoCSV an/or output the metadata in JSON.

    Call arguments:
        -h, --help: this message.
        -i, --input: [required] the input nefiletCDF filename.
        -g, --geocsv: [true/false, default false] output the GeoCSV. The output will have the same name as
              the input file.
        -m, --metadata: [true/false, default true] output metadata in JSON,
              it will have the same filename as the input file.
"""
    )


def main():
    # Capture the input parameters.
    try:
        argv = sys.argv[1:]
        opts, args = getopt.getopt(
            argv, "hg:m:i:", ["help", "geocsv=", "metadata=", "input="]
        )
    except getopt.GetoptError as err:
        # Print the error, and help information and exit:
        logger.error(err)
        usage()
        sys.exit(2)
    # Initialize the variables.
    input_file = None
    output_file = None
    do_geocsv = False
    do_metadata = False
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
        elif o in ("-g", "--geocsv"):
            if a.lower() == "true":
                do_geocsv = True
            elif a.lower() == "false":
                do_geocsv = False
            else:
                usage()
                logger.error(f"[ERR] invalid -g option {a}.")
                sys.exit(1)
        elif o in ("-m", "--metadata"):
            if a.lower() == "true":
                do_metadata = True
            elif a.lower() == "false":
                do_metadata = False
            else:
                usage()
                logger.error(f"[ERR] invalid -m option {a}.")
                sys.exit(1)
        else:
            assert False, "unhandled option"

    # The input filename is required.
    if input_file is None:
        usage()
        logger.error("[ERR] missing -i or --input.")
        sys.exit(1)
    logger.info(f"[INFO] Working on input {input_file}")

    if not do_metadata and not do_geocsv:
        usage()
        logger.error("[ERR] No Heo.")
        sys.exit(1)
    # Set the output filename the same as the input file.
    output_file = os.path.splitext(input_file)[0]
    logger.info(f"[INFO] Will write the output to {output_file}")

    # Figure out the input's file type.
    file_type = lib.check_file_type(input_file)

    if file_type["engine"] == "netcdf" and file_type["valid"] == True:
        logger.info(f"[INFO] {input_file} is a netCDF file")
        # Convert netCDF to GeoCSV.
        if do_geocsv:
            metadata, data = convert_lib.netcdf_to_geocsv(input_file)
            with open(f"{output_file}{prop.extension['geocsv']}", "w") as outfile:
                outfile.write(f"{metadata}\n{data}")

        # Save metadata as JSON
        if do_metadata:
            output_json = convert_lib.json_metadata(input_file)
            with open(f"{output_file}{prop.extension['json']}", "w") as outfile:
                outfile.write(f"{json.dumps(output_json, indent=4)}")
    else:
        logger.error(f"[ERR] {input_file} is not a valid netCDF file")
        sys.exit(1)


if __name__ == "__main__":
    main()
