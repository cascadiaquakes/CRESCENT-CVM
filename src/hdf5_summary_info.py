#!/usr/bin/env python

import h5py
import numpy as np


def summarize_hdf5(file_path):
    # Open the HDF5 file using h5py
    try:
        with h5py.File(file_path, "r") as hdf:
            print(f"Summary of HDF5 file: {file_path}")

            # Function to recursively print group and dataset information
            def print_group_info(name, obj):
                if isinstance(obj, h5py.Group):
                    print(f"\nGroup: {name}")
                    # Print group attributes
                    for attr in obj.attrs:
                        print(f"  - Attribute {attr}: {obj.attrs[attr]}")
                elif isinstance(obj, h5py.Dataset):
                    print(f"\nDataset: {name} (shape: {obj.shape}, dtype: {obj.dtype})")
                    # Print dataset attributes
                    for attr in obj.attrs:
                        print(f"  - Attribute {attr}: {obj.attrs[attr]}")
                    # Print data range
                    try:
                        data = obj[()]
                        if data.size > 0:
                            non_nan_data = data[~np.isnan(data)]
                            if non_nan_data.size > 0:
                                print(
                                    f"  Range: {np.min(non_nan_data):.2f} to {np.max(non_nan_data):.2f}"
                                )
                            else:
                                print(f"  Range: All values are NaN")
                        else:
                            print(f"  Range: No data")
                    except Exception as e:
                        print(f"  Could not read data: {e}")

            # Walk through the file and print information
            hdf.visititems(print_group_info)

    except Exception as e:
        print(f"Error opening file: {e}")


def main():
    # Read the HDF5 file path from user input
    file_path = input("Enter the path to the HDF5 file: ")
    summarize_hdf5(file_path)


if __name__ == "__main__":
    main()
