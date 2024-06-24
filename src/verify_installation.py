#!/usr/bin/env python

import importlib
import subprocess
import sys
import os

# Path to the requirements.txt file
requirements_path = os.path.join(os.path.dirname(__file__), '..', 'requirements.txt')

# Path to the NetCDF file
netcdf_file_path = os.path.join(os.path.dirname(__file__), '..', 'sample-files', 'Cascadia-ANT+RF-Delph2018', 'Cascadia-ANT+RF-Delph2018.r0.1.nc')

def read_requirements(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip() and not line.startswith('#')]

def check_packages(packages):
    for package in packages:
        package_name = package.split('==')[0]
        try:
            importlib.import_module(package_name)
            print(f"{package_name} is installed correctly.")
        except ImportError:
            print(f"{package_name} is not installed. Attempting to install it now.")
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])

def basic_functionality_test():
    try:
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import netCDF4

        # Basic numpy test
        array = np.array([1, 2, 3])
        assert np.sum(array) == 6

        # Basic pandas test
        df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]})
        assert df.shape == (3, 2)

        # Basic matplotlib test
        plt.plot(array)
        plt.title('Test Plot')
        plt.show()

        print("All basic functionality tests passed.")

    except Exception as e:
        print(f"An error occurred during basic functionality tests: {e}")
        return False
    return True

def netcdf_summary(file_path):
    try:
        from netCDF4 import Dataset
        dataset = Dataset(file_path, 'r')
        print(f"NetCDF file: {file_path}")
        print(f"Dimensions: {dataset.dimensions.keys()}")
        print(f"Variables: {dataset.variables.keys()}")
        dataset.close()
    except Exception as e:
        print(f"An error occurred while reading the NetCDF file: {e}")
        return False
    return True

if __name__ == "__main__":
    print("Reading requirements.txt...")
    required_packages = read_requirements(requirements_path)
    print("Checking required packages...")
    check_packages(required_packages)
    print("Running basic functionality tests...")
    if basic_functionality_test():
        print("Checking NetCDF file...")
        if netcdf_summary(netcdf_file_path):
            print("Installation confirmed. All tests passed successfully.")
        else:
            print("Installation check failed during NetCDF summary.")
    else:
        print("Installation check failed during basic functionality tests.")

