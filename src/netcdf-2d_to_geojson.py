#!/usr/bin/env python
import os
import netCDF4 as nc
import geojson
import numpy as np
import json


class NumpyEncoder(json.JSONEncoder):
    """Use a custom JSON encoder that knows how to handle numpy.float32."""

    def default(self, obj):
        if isinstance(obj, np.float32):
            return float(obj)
        return super().default(obj)


def prompt_for_input(prompt, default=None, valid_options=None):
    """Prompts the user for input and provides a default value. It can also check against valid options."""
    while True:
        # If valid_options is provided, create the string for valid options, else make it empty
        options_str = f"[{' / '.join(valid_options)}]" if valid_options else ""
        user_input = input(f"{prompt} {options_str} (default: {default}): ") or default
        if valid_options and user_input not in valid_options:
            print(f"Invalid input. Please choose from {valid_options}.")
        else:
            return user_input


def validate_file(file_path):
    """Validates whether the provided file path exists."""
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    return file_path


def grd_to_geojson_surface(
    grd_file,
    output_file,
    lat_var,
    lon_var,
    z_var,
    depth_unit="km",
    minify=False,
    decimation_factor=2,
    grid_type="rectangular",
):
    """Converts a 2D NetCDF (.nc) or .grd file to a GeoJSON surface."""

    # Open the grd (netCDF) file
    try:
        dataset = nc.Dataset(grd_file, "r")
    except FileNotFoundError:
        raise FileNotFoundError(f"File '{grd_file}' could not be found or opened.")

    # Validate if variables exist in the file
    if (
        lat_var not in dataset.variables
        or lon_var not in dataset.variables
        or z_var not in dataset.variables
    ):
        raise ValueError(
            "One or more variables are not found in the file. Please check variable names."
        )

    # Read the latitude, longitude, and z-variable data
    try:
        lats = dataset.variables[lat_var][::decimation_factor]
        lons = dataset.variables[lon_var][::decimation_factor]
        z_data = dataset.variables[z_var][::decimation_factor, ::decimation_factor]
    except KeyError as e:
        raise KeyError(f"Failed to read variables: {e}")

    # Check if the dimensions are aligned correctly
    if len(lats.shape) != 1 or len(lons.shape) != 1:
        raise ValueError("Latitude and Longitude should be 1D arrays")
    if z_data.shape != (len(lats), len(lons)):
        raise ValueError(
            f"Shape of z-variable data {z_data.shape} does not match the grid dimensions"
        )

    # Prepare a list of GeoJSON features
    features = []

    # Iterate over the decimated grid and create polygons for each surface patch
    for i in range(len(lats) - 1):
        for j in range(len(lons) - 1):
            # Get the four corners of the grid cell
            lat1, lon1, z1 = lats[i], lons[j], z_data[i, j]
            lat2, lon2, z2 = lats[i + 1], lons[j], z_data[i + 1, j]
            lat3, lon3, z3 = lats[i + 1], lons[j + 1], z_data[i + 1, j + 1]
            lat4, lon4, z4 = lats[i], lons[j + 1], z_data[i, j + 1]

            # Skip if any z-value is masked or invalid
            if (
                np.ma.is_masked(z1)
                or np.ma.is_masked(z2)
                or np.ma.is_masked(z3)
                or np.ma.is_masked(z4)
            ):
                continue

            # If z represents depth, convert it to meters and make it negative (depth is positive down)
            if depth_unit == "km":
                z1, z2, z3, z4 = -z1 * 1000, -z2 * 1000, -z3 * 1000, -z4 * 1000

            # Create the grid based on the selected grid type
            if grid_type == "rectangular":
                # Create a single rectangle polygon for the grid cell
                rectangle = geojson.Polygon(
                    [
                        [
                            (lon1, lat1, z1),
                            (lon2, lat2, z2),
                            (lon3, lat3, z3),
                            (lon4, lat4, z4),
                            (lon1, lat1, z1),
                        ]
                    ]
                )
                features.append(
                    geojson.Feature(
                        geometry=rectangle, properties={"elevation": [z1, z2, z3, z4]}
                    )
                )

            elif grid_type == "triangular":
                # Create two triangles for each grid cell (split the rectangle into two triangles)
                triangle1 = geojson.Polygon(
                    [
                        [
                            (lon1, lat1, z1),
                            (lon2, lat2, z2),
                            (lon3, lat3, z3),
                            (lon1, lat1, z1),
                        ]
                    ]
                )
                triangle2 = geojson.Polygon(
                    [
                        [
                            (lon1, lat1, z1),
                            (lon3, lat3, z3),
                            (lon4, lat4, z4),
                            (lon1, lat1, z1),
                        ]
                    ]
                )
                # Add the two triangles as separate features
                features.append(
                    geojson.Feature(
                        geometry=triangle1, properties={"elevation": [z1, z2, z3]}
                    )
                )
                features.append(
                    geojson.Feature(
                        geometry=triangle2, properties={"elevation": [z1, z3, z4]}
                    )
                )

            else:
                raise ValueError(
                    "Invalid grid type. Use 'rectangular' or 'triangular'."
                )

    # Create the FeatureCollection
    feature_collection = geojson.FeatureCollection(features)

    # Write to the output GeoJSON file
    with open(output_file, "w") as f:
        if minify:
            # Minify the output by removing whitespace and line breaks
            geojson.dump(feature_collection, f, separators=(",", ":"), cls=NumpyEncoder)
        else:
            # Pretty print with spaces and line breaks for readability
            geojson.dump(feature_collection, f, indent=2, cls=NumpyEncoder)

    print(
        f"GeoJSON surface saved to {output_file} (Minified: {minify}, Grid Type: {grid_type})"
    )


# Prompt the user for input parameters
print("2D netCDF/.grd file to GeoJSON Converter!")

input_file = prompt_for_input(
    "Enter the path to the input NetCDF/GRD file", default="input.nc"
)
validate_file(input_file)

output_file = prompt_for_input(
    "Enter the path to save the GeoJSON output", default="output.geojson"
)

latitude = prompt_for_input(
    "Enter the name of the latitude variable", default="latitude"
)
longitude = prompt_for_input(
    "Enter the name of the longitude variable", default="longitude"
)
depth = prompt_for_input(
    "Enter the name of the depth/elevation variable", default="depth"
)

depth_unit = prompt_for_input(
    "Enter the depth unit", default="km", valid_options=["km", "m"]
)

decimation_factor = int(
    prompt_for_input("Enter the decimation factor (integer)", default="2")
)

grid_type = prompt_for_input(
    "Enter the grid type",
    default="rectangular",
    valid_options=["rectangular", "triangular"],
)

minify = (
    prompt_for_input(
        "Do you want to minify the GeoJSON file?",
        default="no",
        valid_options=["yes", "no"],
    )
    == "yes"
)

# Call the function with user input
grd_to_geojson_surface(
    input_file,
    output_file,
    latitude,
    longitude,
    depth,
    depth_unit=depth_unit,
    minify=minify,
    decimation_factor=decimation_factor,
    grid_type=grid_type,
)