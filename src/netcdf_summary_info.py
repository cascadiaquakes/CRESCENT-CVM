#!/usr/bin/env python

import netCDF4 as nc
import numpy as np

"""A utility tool designed to read and analyze NetCDF4 classic files. It provides a detailed summary of the file content, including global metadata, dimensions, and variables with their respective attributes and data ranges. This script is essential for verifying the overall structure and integrity of NetCDF files, ensuring that they conform to expected formats and contain the necessary
information for further processing or analysis.

    Call arguments:
        None, code will ask for the input netCDF model file name.
"""


def summarize_netcdf(file_path):
    # Open the NetCDF file using netCDF4
    try:
        dataset = nc.Dataset(file_path, "r")
    except Exception as e:
        print(f"Error opening file: {e}")
        return

    print(f"Summary of NetCDF file: {file_path}")

    # Global metadata
    print("\nGlobal Attributes:")
    for attr in dataset.ncattrs():
        print(f"  - {attr}: {dataset.getncattr(attr)}")

    # Dimensions
    print("\nDimensions:")
    for dim_name, dim in dataset.dimensions.items():
        print(f"  - {dim_name}: {len(dim)}")

    # Variables and their metadata
    print("\nVariables:")
    for var_name, var in dataset.variables.items():
        print(f"  - {var_name} (dimensions: {var.dimensions}, shape: {var.shape})")
        for attr in var.ncattrs():
            print(f"    - {attr}: {var.getncattr(attr)}")

        # Range of variable values
        try:
            data = var[:]
            if data.size > 0:  # Avoid issues with variables that have no data
                non_nan_data = data[~np.isnan(data)]
                if non_nan_data.size > 0:
                    print(
                        f"    Range: {np.min(non_nan_data):.2f} to {np.max(non_nan_data):.2f}"
                    )
                else:
                    print(f"    Range: All values are NaN")
            else:
                print(f"    Range: No data")
        except Exception as e:
            print(f"    Could not read data: {e}")

    dataset.close()


def main():
    # Read the NetCDF file path from user input
    file_path = input("Enter the path to the NetCDF file: ")
    summarize_netcdf(file_path)


if __name__ == "__main__":
    main()
