#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import json
import xarray as xr
import numpy as np
import pandas as pd
import sys
import getopt

# Get the root _directory.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop
import writer_prop as writer_prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib

logger = lib.get_logger()


def geojson_to_csv(geojson_data, output_csv):
    """
    Converts GeoJSON data to a CSV file with x, y, z coordinates,
    avoiding duplicate rows and ensuring proper handling of elevation (z).

    Parameters:
        geojson_data (dict): Parsed GeoJSON data.
        output_csv (str): Path to save the converted CSV file.
    """
    features = geojson_data.get("features", [])
    if not features:
        print("No features found in GeoJSON.")
        return

    rows = []
    for feature in features:
        properties = feature.get("properties", {})
        geometry = feature.get("geometry", {})
        coordinates = geometry.get("coordinates", [])

        # Handle different geometry types
        if geometry.get("type") == "Point":
            x, y, z = coordinates + [np.nan] * (
                3 - len(coordinates)
            )  # Pad missing z with NaN
            row = properties.copy()
            row.update({"x": x, "y": y, "z": z})
            rows.append(row)
        elif geometry.get("type") == "LineString":
            # Process each coordinate in LineString
            for coord in coordinates:
                x, y, z = coord + [np.nan] * (3 - len(coord))
                row = properties.copy()
                row.update({"x": x, "y": y, "z": z})
                rows.append(row)
        elif geometry.get("type") == "Polygon":
            # Process each point in all rings of the polygon
            for ring in coordinates:
                for point in ring:
                    x, y, z = point + [np.nan] * (3 - len(point))
                    row = properties.copy()
                    row.update({"x": x, "y": y, "z": z})
                    rows.append(row)
        else:
            print(f"Unsupported geometry type: {geometry.get('type')}")

    # Create DataFrame and handle unhashable types
    df = pd.DataFrame(rows)
    df = df.replace({None: np.nan, "": np.nan})

    # Convert unhashable list-like values to strings for hashable operations
    for col in df.columns:
        df[col] = df[col].apply(lambda v: str(v) if isinstance(v, list) else v)

    # Remove duplicates
    df = df.drop_duplicates()

    # Save to CSV
    df.to_csv(output_csv, index=False, na_rep="NaN")
    print(f"GeoJSON converted to CSV: {output_csv}")


def netcdf_to_csv(dataset, output_csv):
    """
    Converts a NetCDF dataset to a CSV file.

    Parameters:
        dataset (xarray.Dataset): Opened NetCDF dataset.
        output_csv (str): Path to save the converted CSV file.
    """
    df = dataset.to_dataframe().reset_index()

    # Replace missing values with NaN
    df = df.replace({None: np.nan, "": np.nan})  # Handle both None and empty strings
    df = df.fillna(np.nan)  # Ensure all blanks are converted to NaN

    # Save to CSV with NaN representation
    df.to_csv(output_csv, index=False, na_rep="NaN")
    logger.info(f"[INFO] NetCDF converted to CSV: {output_csv}")


def determine_file_type(input_file):
    """
    Determines the type of file based on its content.

    Parameters:
        input_file (str): Path to the input file.

    Returns:
        str: File type ('geojson' or 'netcdf').
    """
    try:
        # Try reading as GeoJSON
        with open(input_file, "r") as f:
            json_data = json.load(f)
        if "features" in json_data and "type" in json_data:
            return "geojson", json_data
    except (json.JSONDecodeError, UnicodeDecodeError):
        pass

    try:
        # Try reading as NetCDF
        dataset = xr.open_dataset(input_file)
        return "netcdf", dataset
    except Exception:
        pass

    return None, None


def convert_to_csv(input_file, output_csv):
    """
    Converts a GeoJSON or NetCDF file to a CSV file.

    Parameters:
        input_file (str): Path to the input file.
        output_csv (str): Path to save the converted CSV file.
    """
    file_type, content = determine_file_type(input_file)
    if file_type == "geojson":
        geojson_to_csv(content, output_csv)
    elif file_type == "netcdf":
        netcdf_to_csv(content, output_csv)
    else:
        logger.error(
            "[ERR] Unsupported file format. Please provide a valid GeoJSON or netCDF file."
        )


def usage():
    """
    Prints the usage instructions for the script.
    """
    print(
        """
Usage: python convert_to_csv.py -i <input_file> -o <output_file>
Options:
  -i, --input    Input file path (GeoJSON or NetCDF file).
  -o, --output   Output file path (CSV file).
  -h, --help     Display this help message and exit.

Example:
  python convert_to_csv.py -i data.geojson -o output.csv
  python convert_to_csv.py -i data.nc -o output.csv
  python convert_to_csv.py -i data.grd -o output.csv
    """
    )


if __name__ == "__main__":
    # Parse command-line arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:h", ["input=", "output=", "help"])
    except getopt.GetoptError as err:
        logger.error(f"[ERR] {err}")
        usage()
        sys.exit(2)

    input_file = None
    output_file = None

    # Process command-line arguments
    for opt, arg in opts:
        if opt in ("-i", "--input"):
            input_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg
        elif opt in ("-h", "--help"):
            usage()
            sys.exit(0)

    # Check if input and output files are provided
    if not input_file or not output_file:
        logger.error("[ERR] Error: Input and output file paths are required.")
        usage()
        sys.exit(2)

    # Perform the conversion
    convert_to_csv(input_file, output_file)
