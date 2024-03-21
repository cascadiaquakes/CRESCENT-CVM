import os
import sys
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import matplotlib.pyplot as plt

# Get the root _directory.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib

# Set up the logger.
logger = lib.get_logger()


def get_point(location):
    f"""
    Get the coordinates for a point.

    Call arguments:
        location: point's location.
    """
    point = None
    while point is None:
        if location != "depth":
            point = input(
                f"Cross-section {location} point as lat,lon ['back', 'exit']? "
            )
        else:
            point = input(f"Cross-section depth range as start, end ['back', 'exit']? ")

        # Done, back!
        if point.strip() == "exit":
            sys.exit()
        elif point.strip() == "back" or not point.strip():
            return "back"
        elif "," not in point:
            if location != "depth":
                logger.error(
                    f"[ERR] invalid {location} coordinates '{point}' input {location} as lat,lon"
                )
            else:
                logger.error(
                    f"[ERR] invalid {location} range '{point}' input {location} as start, end depths"
                )
            point = None
        else:
            try:
                point = point.split(",")
                point = [float(i) for i in point]
                break
            except Exception as ex:
                logger.error(
                    f"[ERR] invalid {location} range '{point}' input {location} as value1,value2\n{ex}"
                )
            point = None
    return point


def display_var_meta(ds, meta, indent, values=True):
    """
    Display metadata fora given variable.

    Call arguments:
        ds - [required] the xarray dataset
        meta - [required] dataset metadata
        indent - [required] display indent
        values - [optional, default: True] display the values
    """
    for var in meta:
        for attr_indx, attr in enumerate(ds[var].attrs):
            if attr_indx == 0:
                logger.info(f"\t{var}: ")
            if attr == "variable":
                ConnectionRefusedError
            else:
                logger.info(f"{indent*' '}{attr}: {ds[var].attrs[attr]}")
        if values:
            logger.info(f"{indent*' '}Values:")
            lib.view_list(ds[var].data, indent=indent)


def slicer(ds, slice_dir, slice_value, slice_limits):
    """
    Slice a dataset along a given variable.

    Call arguments:
        ds - [required] the xarray dataset
        slice_dire - [required] the variable name along which to slice
        slice_value - [required] the slice location
        slice_limits - [required] limits of the slice in other directions.
    """
    slice_limits_keys = list(slice_limits.keys())
    slice_limits_values = list(slice_limits.values())

    sliced_data = ds.where(
        (ds[slice_limits_keys[0]] >= slice_limits_values[0][0])
        & (ds[slice_limits_keys[0]] <= slice_limits_values[0][1])
        & (ds[slice_limits_keys[1]] >= slice_limits_values[1][0])
        & (ds[slice_limits_keys[1]] <= slice_limits_values[1][1]),
        drop=True,
    )
    sliced_data = sliced_data.sel({slice_dir: float(slice_value)})
    return sliced_data


def subsetter(ds, limits):
    """
    Subset a dataset as a volume.

    Call arguments:
        ds - [required] the xarray dataset
        limits - [required] limits of the volume in all directions.
    """
    limits_keys = list(limits.keys())
    limits_values = list(limits.values())

    sliced_data = ds.where(
        (ds[limits_keys[0]] >= limits_values[0][0])
        & (ds[limits_keys[0]] <= limits_values[0][1])
        & (ds[limits_keys[1]] >= limits_values[1][0])
        & (ds[limits_keys[1]] <= limits_values[1][1])
        & (ds[limits_keys[2]] >= limits_values[2][0])
        & (ds[limits_keys[2]] <= limits_values[2][1]),
        drop=True,
    )
    return sliced_data


def gmap(plot_var, cmap, gmap_limits, sliced_data):
    """
    Geographical display of a dataset.

    Call arguments:
        data_var - [required] the variable for which to create the plot.
        cmap - [required] colormap to use.
        gmap_limits - [required] surface geographical limits.
        slice_limits - [required] sliced data to plot.
    """
    # Color from cmap.
    color = cmap

    # Defining the figure
    fig = plt.figure(facecolor="w", edgecolor="k")

    # Axes with Cartopy projection
    ax = plt.axes(projection=ccrs.PlateCarree())
    # and extent
    ax.set_extent(
        [
            gmap_limits["longitude"][0],
            gmap_limits["longitude"][1],
            gmap_limits["latitude"][0],
            gmap_limits["latitude"][1],
        ],
        ccrs.PlateCarree(),
    )

    # Plotting using Matplotlib
    cf = sliced_data[plot_var].plot(
        x="longitude",
        y="latitude",
        transform=ccrs.PlateCarree(),
        cmap=color,
        add_colorbar=True,
    )

    # Plot lat/lon grid
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.1,
        color="k",
        alpha=1,
        linestyle="--",
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Add map features with Cartopy
    ax.add_feature(
        cfeature.NaturalEarthFeature(
            "physical",
            "land",
            "10m",
            edgecolor="face",
            facecolor="none",
        )
    )
    ax.coastlines(linewidth=1)

    plt.tight_layout()
    plt.show()
