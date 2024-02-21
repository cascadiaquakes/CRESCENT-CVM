#!/usr/bin/env python

import sys
import os
import getopt
import matplotlib.pyplot as plt
import xarray as xr

# Get the directory paths.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LIB_DIR = os.path.join(ROOT_DIR, "lib")
sys.path.append(LIB_DIR)
import shared_lib as lib

sys.path.append(os.getcwd())
import simple_plotter_prop as prop

"""A simple Python code to demonstrate the ease of access Xarray offers
for working the CVM files in netCDF format.
"""

divider = 80 * "="
logger = lib.get_logger()

# Initialize the variables.
filename = prop.filename
depth = prop.depth
variable = prop.variable


def usage():
    """Script usage notes."""
    logger.info(
        f"""
    Script to plot a simple depth slice from a given netCDF model file. It expects to find the 
    plotter property file (simple_plotter_prop.py) under the working directory.
        -m, --model  [default: {prop.filename}] model filename
        -d, --depth  [default: {prop.depth}] slice depth
        -v, --variable [default: {prop.variable}] the variable to plot
        """
    )


# Capture the input parameters.
try:
    argv = sys.argv[1:]
    opts, args = getopt.getopt(
        argv, "hm:d:v:", ["help", "model=", "depth=", "variable="]
    )
except getopt.GetoptError as err:
    # Print the error, and help information and exit:
    logger.error(err)
    sys.exit(1)


for o, a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit(0)
    elif o in ("-m", "--model"):
        filename = a
        if not os.path.isfile(filename):
            usage()
            logger.error(f"[ERR] Could not find the input file: {filename}")
            sys.exit(1)
    elif o in ("-d", "--depth"):
        try:
            depth = float(a)
        except Exception as ex:
            usage()
            logger.error(f"[ERR] Invalid depth (must be float): {depth}\n{ex}")
            sys.exit(1)
    elif o in ("-v", "--variable"):
        variable = a
    else:
        usage()
        logger.error(f"[ERR] Invalid option {o} with value of {a}")
        sys.exit(1)


divider = 80 * "="
vmin = prop.vmin
vmax = prop.vmax
figure_size = prop.figure_size
figure_size_s = prop.figure_size_s
cmap = prop.cmap
x = prop.x
y = prop.y
auxiliary_x = prop.auxiliary_x
auxiliary_y = prop.auxiliary_y
y_min = prop.y_min
x_min = prop.x_min

# Read the netCDF content to an xarray dataset and display its content.
ds = xr.open_dataset(filename)
logger.info(f"\n\n\n{filename} dataset content:\n{divider}\n{ds}")

# Extract the designated variable's dataset and display its content.
ds_var = ds[variable]
logger.info(f"\n\n{variable} variable's dataset content:\n{divider}\n{ds_var}")

# Extract the horizontal (depth) slice for the give depth.
ds_var_depth = ds_var.where(ds_var.depth == depth, drop=True)
ds_var_depth.plot(figsize=figure_size, cmap=cmap, vmin=vmin, vmax=vmax)
plt.show()

# Plot the same slice in the auxiliary coordinates.
ds_var_depth.plot(
    figsize=figure_size, cmap=cmap, vmin=vmin, vmax=vmax, x=auxiliary_x, y=auxiliary_y
)
plt.show()


# extract the uper right portion of the slice.
ds_var_depth_ur = ds_var.where(
    (ds_var.depth == depth) & (ds_var[x] >= x_min) & (ds_var[y] >= y_min), drop=True
)
ds_var_depth_ur.plot(figsize=figure_size_s, cmap=cmap, vmin=vmin, vmax=vmax)
plt.show()
