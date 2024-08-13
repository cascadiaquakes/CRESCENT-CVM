#!/usr/bin/env python

import h5py
import netCDF4 as nc
import os
import sys

def print_structure(name, obj, indent=''):
    """
    Recursively prints the structure of HDF5 groups and datasets.

    Parameters:
    - name (str): The name of the current group or dataset.
    - obj (h5py.Group or h5py.Dataset): The current HDF5 object.
    - indent (str): The indentation level for printing nested structures.
    
    """
    if isinstance(obj, h5py.Group):
        print(f"{indent}├── {name}")
        indent += '│   '
        for key, val in obj.items():
            print_structure(key, val, indent)
    elif isinstance(obj, h5py.Dataset):
        print(f"{indent}├── {name}")

def print_hdf5_structure(hdf5_file):
    """
    Prints the structure of an HDF5 file.

    Parameters:
    - hdf5_file (str): The path to the HDF5 file.
    """
    with h5py.File(hdf5_file, 'r') as f:
        print("HDF5 File Structure:")
        for key, val in f.items():
            print_structure(key, val)

def print_netcdf_structure(netcdf_file):
    """
    Prints the structure of a netCDF file.

    Parameters:
    - netcdf_file (str): The path to the netCDF file.
    """
    with nc.Dataset(netcdf_file, 'r') as ds:
        print("netCDF File Structure:")
        print("Dimensions:")
        for dim_name, dim in ds.dimensions.items():
            print(f"├── {dim_name}: {len(dim)}")
        print("Variables:")
        for var_name, var in ds.variables.items():
            print(f"├── {var_name}: {var.dimensions}")

def main():
    """
    Main function that prompts the user for the input file name, prints the file name, size, and type,
    and then prints the file structure. Supports both HDF5 and netCDF file formats.

    Author: Manochehr Bahavar / EarthScope
    """
    input_file = input("Enter the path to the input file (HDF5 or netCDF): ")

    if not os.path.isfile(input_file):
        print(f"Error: File '{input_file}' does not exist.")
        sys.exit(1)

    file_size = os.path.getsize(input_file)
    file_type = None

    if input_file.endswith('.h5') or input_file.endswith('.hdf5'):
        file_type = "HDF5"
    elif input_file.endswith('.nc'):
        file_type = "netCDF"
    else:
        print(f"Error: Unsupported file type for file '{input_file}'.")
        sys.exit(1)

    print(f"File Name: {os.path.basename(input_file)}")
    print(f"File Size: {file_size} bytes")
    print(f"File Type: {file_type}")
    print()

    if file_type == "HDF5":
        print_hdf5_structure(input_file)
    elif file_type == "netCDF":
        print_netcdf_structure(input_file)

if __name__ == "__main__":
    main()

