#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import netCDF4
import traceback

# Get the root _directory.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib

passed_list = ["\t"]
failed_list = ["\t"]
warn_list = ["\t"]
linebreak = "\n"
linebreak_tab = f"{linebreak}\t"

if lib.supports_unicode():
    checkmark = "✔️".encode("utf-8").decode("utf-8")
    failmark = "❌".encode("utf-8").decode("utf-8")
    notemark = "⚠️".encode("utf-8").decode("utf-8")
else:
    checkmark = "[OK]"
    failmark = "[X]"
    notemark = "[!]"

STANDARDS = (
    "https://cascadiaquakes.github.io/cvm-tools-book/standards_and_conventions.html"
)
VAR_STANDARDS = "https://cascadiaquakes.github.io/cvm-tools-book/standards_and_conventions.html#model-variable-naming-and-unit-requirements"
GLOBAL_STANDARDS = "https://cascadiaquakes.github.io/cvm-tools-book/usage/global_metadata.html#global-metadata-file-structure"
VAR_STANDARDS = "https://cascadiaquakes.github.io/cvm-tools-book/usage/variables_metadata.html#variable-s-metadata-file-structure"

# Preferred order of dimensions.
preferred_dimensions_order = {
    "latitude": ["latitude", "longitude"],
    "longitude": ["latitude", "longitude"],
    "x": ["x", "y"],
    "y": ["x", "y"],
}

# List of required global attributes
required_global_attributes = [
    "title",
    "id",
    "model",
    "data_revision",
    "summary",
    "Conventions",
    "Metadata_Conventions",
    "repository_name",
    "repository_institution",
    "repository_pid",
    "author_name",
    "author_institution",
    "reference",
    "reference_pid",
    "grid_ref",
    "utm_zone",
    "ellipsoid",
    "data_layout",
]

# List of optional global attributes
optional_global_attributes = [
    "author_email",
    "author_url",
]

metadata_summary = {}

# Required attributes based on variable type
required_attributes = {
    "coordinates": ["long_name", "units", "standard_name"],
    "auxiliary": ["long_name", "units", "standard_name"],
    "depth": ["long_name", "units", "positive"],
    "model": [
        "long_name",
        "units",
        "display_name",
        "data_source",
    ],
}

accepted_file_type = ["NetCDF-4", "NetCDF-4 Classic Model"]


def check_netcdf_file(file_name):
    try:
        file_size = os.path.getsize(file_name)  # Get the file size

        # Try opening the file with netCDF4
        dataset = netCDF4.Dataset(file_name, "r")

        try:
            # Determine the format
            if dataset.file_format == "NETCDF3_CLASSIC":
                message = f"NetCDF-3 Classic {failmark} (see [1])"
                metadata_summary["File Type"] = message
                failed_list.append(f"Bad file format: NetCDF-3 Classic (see [1])")
            elif dataset.file_format == "NETCDF3_64BIT":
                message = f"NetCDF-3 64-bit offset {failmark} (see [1])"
                metadata_summary["File Type"] = message
                failed_list.append(
                    f"Bad file format: NetCDF-3  64-bit offset  (see [1])"
                )
            elif dataset.file_format == "NETCDF4_CLASSIC":
                metadata_summary["File Type"] = f"NetCDF-4 Classic Model {checkmark}"
                passed_list.append(f"File format: NetCDF-4 Classic Model ")
            elif dataset.file_format == "NETCDF4":
                metadata_summary["File Type"] = f"NetCDF-4 {checkmark}"
                passed_list.append(f"File format: NetCDF-4")
            else:
                metadata_summary["File Type"] = "Unknown"
                file_type_marker = f"{failmark} (see [1])"
                output(f"{file_type_marker} file format (see [1])")
                sys.exit(1)

        except Exception as e:
            output(f"Error determining file format: {e} {failmark}")
            metadata_summary["File Type"] = (
                f"Error determining file format: {e} {failmark}"
            )
            sys.exit(1)

        output(
            f"\n\nNOTE:  Checks based on Standards and Conventions (see [1])\n\n=== NetCDF File Check ===\nFile: {file_name}\nSize: {file_size / (1024 ** 2):.2f} MB\nFile Type: {metadata_summary['File Type']}"
        )

        try:
            # List dimensions and their sizes
            output_header("Dimensions and Sizes")
            list_dimensions(dataset)

            # List primary and auxiliary coordinates
            output_header("Coordinate Information")
            primary_coords, auxiliary_coords = list_coordinates(dataset)

            # Latitude and Longitude Value Check
            output_header("Latitude and Longitude Value Check")
            lat_min, lat_max, lon_min, lon_max = check_lat_lon_values(
                dataset, primary_coords, auxiliary_coords
            )

            # Check geospatial attributes
            output_header("Geospatial Attributes Check")
            check_geospatial_attributes(dataset, lat_min, lat_max, lon_min, lon_max)

            # Check for required and optional global attributes
            output_header("Global Attributes Check")
            check_global_attributes(dataset)

            # Perform CF compliance check
            output_header("CF Compliance Check")
            check_cf_compliance(file_name)

            # List variables with their ranges and units, split into types
            output_header("Variables Summary")
            list_variables(dataset, primary_coords)

        except Exception as e:
            output(f"Error during processing: {e} {failmark}")
            metadata_summary["Processing Error"] = (
                f"Error during processing: {e} {failmark}\nTraceback:{traceback.print_exec()}"
            )

        dataset.close()

    except OSError as e:
        output(f"Error: Unable to open the file '{file_name}'. {e} {failmark}")
        metadata_summary["File Access Error"] = (
            f"Unable to open the file '{file_name}'. {e} {failmark}"
        )

    # Print metadata summary
    try:
        output_header(
            "Metadata Inspection Summary"
        )  # Changed to 'Metadata Inspection Summary'
        # print_metadata_summary()
        output(
            f"{checkmark} Passed list:\n{linebreak_tab.join(passed_list)}"
            f"\n\n{failmark} Failed list:\n{linebreak_tab.join(failed_list)}"
            f"\n\n{notemark}  Warning list:\n{linebreak_tab.join(warn_list)}"
        )
        output(f"\n\n\n")

    except Exception as e:
        output(f"Error while printing metadata summary: {e} {failmark}")


def list_dimensions(dataset):
    try:
        dimensions = {dim: len(dataset.dimensions[dim]) for dim in dataset.dimensions}
        metadata_summary["Dimensions"] = dimensions
        for dim, size in dimensions.items():
            message = f"{dim}: {size}"
            output(f"- {message} {checkmark}", indent=True)
            passed_list.append(message)
    except Exception as e:
        output(f"Error listing dimensions: {e} {failmark}")
        metadata_summary["Dimensions"] = f"Error listing dimensions: {e} {failmark}"


def list_coordinates(dataset):
    try:
        primary_coords = []
        auxiliary_coords = set()  # Use a set to avoid duplicate entries
        third_dimension = None

        # Iterate over variables to identify coordinates
        for var_name, variable in dataset.variables.items():
            if hasattr(variable, "coordinates"):
                auxiliary_coords.update(variable.coordinates.split())
            elif var_name in dataset.dimensions:
                primary_coords.append(var_name)

        # Identify the third dimension
        for coord in primary_coords:
            if coord.lower() not in ["latitude", "longitude", "x", "y"]:
                third_dimension = coord
                break

        metadata_summary["Primary Coordinates"] = (
            primary_coords if primary_coords else None
        )
        metadata_summary["Auxiliary Coordinates"] = (
            list(auxiliary_coords) if auxiliary_coords else None
        )
        metadata_summary["Third Dimension"] = (
            third_dimension if third_dimension else None
        )

        if primary_coords:
            message = f"Primary Coordinates: {', '.join(primary_coords)}"
            output(
                f"{message} {checkmark}",
                indent=True,
            )
            passed_list.append(message)
            if third_dimension:
                message = f"Third Dimension: {third_dimension}"
                output(message, indent=True)
                passed_list.append(message)
        else:
            message = f"No primary coordinates found."
            output(f"{message} {failmark}", indent=True)
            metadata_summary["Primary Coordinates"] = message
            failed_list.append(message)

        if auxiliary_coords:
            message = f"Auxiliary Coordinates: {', '.join(auxiliary_coords)}"
            output(
                f"{message} {checkmark}",
                indent=True,
            )
            passed_list.append(message)
        else:
            output(
                f"No auxiliary coordinates found. {notemark} (Optional)", indent=True
            )
            warn_list.append(f"No optional auxiliary coordinates found.")

        # Return both primary and auxiliary coordinates
        return primary_coords, auxiliary_coords

    except Exception as e:
        output(f"Error listing coordinates: {e} {failmark}")
        metadata_summary["Coordinates"] = f"Error listing coordinates: {e} {failmark}"
        return [], set()  # Return empty lists in case of an error


def check_lat_lon_values(dataset, primary_coords, auxiliary_coords):
    try:
        lat_name = None
        lon_name = None

        # Identify latitude and longitude variable names
        for coord in primary_coords + list(auxiliary_coords):
            if "lat" in coord.lower():
                lat_name = coord
            if "lon" in coord.lower():
                lon_name = coord

        if lat_name and lon_name:
            lat_values = dataset.variables[lat_name][:]
            lon_values = dataset.variables[lon_name][:]

            # Determine longitude format (0-360 or -180/180)
            if lon_values.min() >= 0 and lon_values.max() <= 180:
                lon_format = "Potentially 0-360 degrees or -180/180 degrees"
            elif lon_values.min() >= 0 and lon_values.max() <= 360:
                lon_format = "0-360 degrees"
            elif lon_values.min() >= -180 and lon_values.max() <= 180:
                lon_format = "-180/180 degrees"
            else:
                lon_format = "Unknown"

            metadata_summary["Longitude Format"] = lon_format
            message = f"Longitude Format: {lon_format}"
            output(f"{message} {checkmark}", indent=True)
            passed_list.append(message)

            # Check that latitudes and longitudes are within valid ranges
            lat_min, lat_max = lat_values.min(), lat_values.max()
            lon_min, lon_max = lon_values.min(), lon_values.max()

            if lat_min < -90 or lat_max > 90:
                message = f"Latitude values are out of range: {lat_min:.3f} to {lat_max:.3f} {failmark}"
                output(
                    message,
                    indent=True,
                )
                metadata_summary["Latitude Range"] = message
                failed_list.append(message)
            else:
                message = f"Latitude Range: {lat_min:.3f} to {lat_max:.3f}"
                output(
                    f"{message}  {checkmark}",
                    indent=True,
                )
                passed_list.append(message)
                metadata_summary["Latitude Range"] = f"{lat_min:.3f} to {lat_max:.3f}"

            if lon_format == "0-360 degrees" and (lon_min < 0 or lon_max > 360):
                message = f"Longitude values are out of range: {lon_min:.3f} to {lon_max:.3f} {failmark}"
                output(
                    message,
                    indent=True,
                )
                metadata_summary["Longitude Range"] = message
                failed_list.append(message)
            elif lon_format == "-180/180 degrees" and (lon_min < -180 or lon_max > 180):
                message = f"Longitude values are out of range: {lon_min:.3f} to {lon_max:.3f} {failmark}"
                output(
                    message,
                    indent=True,
                )
                metadata_summary["Longitude Range"] = message
                failed_list.append(message)
            elif lon_format == "Potentially 0-360 degrees or -180/180 degrees":
                message = f"Longitude Range: {lon_min:.3f} to {lon_max:.3f} {checkmark}"
                output(
                    message,
                    indent=True,
                )
                passed_list.append(message)
                metadata_summary["Longitude Range"] = f"{lon_min:.3f} to {lon_max:.3f}"
            else:
                message = f"Longitude Range: {lon_min:.3f} to {lon_max:.3f} {checkmark}"
                output(
                    message,
                    indent=True,
                )
                metadata_summary["Longitude Range"] = f"{lon_min:.3f} to {lon_max:.3f}"

            return lat_min, lat_max, lon_min, lon_max  # Return computed values

        else:
            message = f"Latitude and/or Longitude not found. {failmark}"
            output(message, indent=True)
            metadata_summary["Coordinates"] = message
            failed_list.append(message)
            return None, None, None, None

    except Exception as e:
        message = f"Error checking latitude and longitude values: {e} {failmark}"
        output(message)
        metadata_summary["Coordinates"] = message
        failed_list.append(message)
        return None, None, None, None


def check_geospatial_attributes(dataset, lat_min, lat_max, lon_min, lon_max):
    try:
        required_geospatial_attrs = [
            "geospatial_lat_min",
            "geospatial_lat_max",
            "geospatial_lon_min",
            "geospatial_lon_max",
            "geospatial_vertical_min",
            "geospatial_vertical_max",
        ]

        missing_attrs = [
            attr for attr in required_geospatial_attrs if attr not in dataset.ncattrs()
        ]

        if missing_attrs:
            for attr in missing_attrs:
                if (
                    "vertical" in attr and "depth" in dataset.dimensions
                ):  # Check if depth is present
                    message = (
                        f"Error: {attr} is missing but depth dimension is present."
                    )
                    output(
                        f"{message} {failmark}",
                        indent=True,
                    )
                    failed_list.append(message)
                elif "vertical" in attr:
                    message = f"Warning: {attr} is missing (model does not have a depth dimension). "
                    output(
                        f"{message} {notemark}",
                        indent=True,
                    )
                    warn_list.append(
                        f"{attr} is missing (model does not have a depth dimension)"
                    )
                else:
                    message = f"Error: Required attribute {attr} is missing."
                    output(
                        f"{message} {failmark}",
                        indent=True,
                    )
                    failed_list.append(message)
            return  # Skip further checks if required attributes are missing

        # Check if depth is a dimension and handle accordingly
        if "depth" in dataset.dimensions:
            depth_values = dataset.variables["depth"][:]
            depth_min = depth_values.min()
            depth_max = depth_values.max()
        else:
            depth_min = None
            depth_max = None

        # Use computed latitude and longitude values from earlier
        geospatial_metadata = {
            "geospatial_lat_min": lat_min,
            "geospatial_lat_max": lat_max,
            "geospatial_lon_min": lon_min,
            "geospatial_lon_max": lon_max,
            "geospatial_vertical_min": (
                depth_min if depth_min is not None else "Not applicable"
            ),
            "geospatial_vertical_max": (
                depth_max if depth_max is not None else "Not applicable"
            ),
        }

        # Compare with geospatial metadata
        for key, value in geospatial_metadata.items():
            if key in dataset.ncattrs():
                metadata_value = dataset.getncattr(key)
                if value == "Not applicable":
                    output(
                        f"{key}: {value} (No depth dimension) {notemark}", indent=True
                    )
                    warn_list.append(f"{key}: {value} (No depth dimension)")
                elif is_within_tolerance(metadata_value, value):
                    message = (
                        f"{key}: Metadata = {metadata_value}, Computed = {value:.3f}"
                    )
                    output(
                        f"{message}  {checkmark}",
                        indent=True,
                    )
                    passed_list.append(message)
                else:
                    message = (
                        f"{key}: Metadata = {metadata_value}, Computed = {value:.3f}"
                    )
                    output(
                        f"{message} {failmark}",
                        indent=True,
                    )
                    failed_list.append(
                        f"{key}: Metadata = {metadata_value}, Computed = {value:.3f} {failmark}"
                    )
                    output(
                        f"Error: {key} mismatch. Metadata = {metadata_value}, Computed = {value:.3f}",
                        indent=True,
                        level=2,
                    )
                    message = f"Metadata = {metadata_value}, Computed = {value:.3f} {failmark}"
                    metadata_summary[key] = message
            else:
                output(
                    f"{key} is not found in metadata. {notemark} (Optional)",
                    indent=True,
                )
                warn_list.append(f"Optional {key} is not found in metadata.")
                if value is not None and value != "Not applicable":
                    message = f"Computed {key}: {value:.3f} {checkmark}"
                    output(message, indent=True)
                    passed_list.append(message)
    except Exception as e:
        message = f"Error checking geospatial attributes: {e} {failmark}"
        output(message)
        metadata_summary["Geospatial Attributes"] = message
        failed_list.append(message)


def is_within_tolerance(metadata_value, computed_value, tolerance=0.1):
    try:
        if metadata_value is None or computed_value is None:
            return False
        return abs(metadata_value - computed_value) <= abs(tolerance * metadata_value)
    except Exception as e:
        output(f"Error during tolerance check: {e} {failmark}")
        return False


def check_global_attributes(dataset):
    try:
        missing_required_attributes = []
        missing_optional_attributes = []

        for attribute in required_global_attributes:
            if attribute not in dataset.ncattrs():
                missing_required_attributes.append(attribute)

        for attribute in optional_global_attributes:
            if attribute not in dataset.ncattrs():
                missing_optional_attributes.append(attribute)

        if missing_required_attributes:
            message = f"Missing Required Global Attributes: {failmark} (see [1] & [3])"
            output(
                message,
                indent=True,
            )
            for attr in missing_required_attributes:
                message = f"- {attr}"
                output(message, indent=True, level=2)
                failed_list.append(
                    f"Missing Required Global Attribute: {attr} (see [1] & [3])"
                )
            metadata_summary["Missing Required Global Attributes"] = (
                missing_required_attributes
            )
        else:
            message = f"All required global attributes are present."
            output(f"{message}  {checkmark}", indent=True)
            passed_list.append(message)

        if missing_optional_attributes:
            output(
                f"Missing Global Attributes (optional): {notemark}",
                indent=True,
            )

            for attr in missing_optional_attributes:
                warn_list.append(f"Missing Optional Global Attribute: {attr}")
                output(f"- {attr}", indent=True, level=2)
            metadata_summary["Missing Optional Global Attributes"] = (
                missing_optional_attributes
            )
        else:
            message = f"All optional global attributes are present."
            output(f"message {checkmark}", indent=True)
            passed_list.append(message)
    except Exception as e:
        output(f"Error checking global attributes: {e} {failmark}")
        metadata_summary["Global Attributes"] = (
            f"Error checking global attributes: {e} {failmark}"
        )


def check_cf_compliance(file_name):
    try:
        # Use the cfchecker command-line tool to check CF compliance
        result = os.popen(f"cfchecks -v CF-1.0 {file_name}").read()

        # Initialize status variables
        errors_detected = 0
        warnings_given = 0
        information_messages = 0

        # Parse the output for specific lines indicating status
        for line in result.splitlines():
            if "ERRORS detected" in line:
                errors_detected = int(line.split(":")[-1].strip())
            elif "WARNINGS given" in line:
                warnings_given = int(line.split(":")[-1].strip())
            elif "INFORMATION messages" in line:
                information_messages = int(line.split(":")[-1].strip())

        # Determine the compliance result based on the parsed information
        if errors_detected == 0 and warnings_given == 0:
            message = f"CF Compliance: Compliant"
            output(f"{message} {checkmark}", indent=True)
            metadata_summary["CF Compliance"] = "Compliant"
            passed_list.append(message)
        elif errors_detected == 0 and warnings_given > 0:
            output(f"CF Compliance: with warning {notemark}", indent=True)
            output(f"Warnings: {warnings_given} {notemark}", indent=True, level=2)
            metadata_summary[f"CF Compliance"] = "with warning {notemark}"
            warn_list.append("CF Compliance")
        elif errors_detected > 0:
            message = f"CF Compliance: Non-compliant"
            output(f"{message} {failmark}", indent=True)
            failed_list.append(message)
            output(f"Issues:\n{result}", indent=True, level=2)
            metadata_summary["CF Compliance"] = f"Non-compliant {failmark}"
    except Exception as e:
        output(f"CF compliance check failed for '{file_name}'. Error: {e} {failmark}")
        metadata_summary["CF Compliance"] = f"Check failed: {e} {failmark}"


def list_variables(dataset, primary_coords):
    try:
        coordinate_variables = {}
        auxiliary_coordinates = set()  # Set to store auxiliary coordinates
        model_variables = {}

        # Retrieve the third dimension from the metadata summary
        third_dimension = metadata_summary.get("Third Dimension")

        # Step 1: Identify all auxiliary coordinates from model variable attributes
        for var_name, variable in dataset.variables.items():
            if "coordinates" in variable.ncattrs():
                coordinates = variable.getncattr("coordinates").split()
                auxiliary_coordinates.update(coordinates)

        for var_name, variable in dataset.variables.items():
            var_type = None
            if var_name in dataset.dimensions:
                if third_dimension == "depth" and var_name == "depth":
                    var_type = "depth"
                else:
                    var_type = "coordinates"
            elif var_name in auxiliary_coordinates:
                var_type = "auxiliary"
            else:
                var_type = "model"

            # Check for required attributes
            missing_attributes = []
            for attr in required_attributes[var_type]:
                if attr not in variable.ncattrs():
                    missing_attributes.append(f"{attr} for {var_name}")

            # Recommend positive down for depth.
            if third_dimension == "depth" and "positive" in variable.ncattrs():
                if variable.positive != "down":
                    output(
                        f"Warning: For {third_dimension} the positive direction is set to '{variable.positive}'\n"
                        f"When using 'depth' as the third dimension, setting it to 'down' is highly recommended because:\n"
                        f"- Most oceanographic and geophysical models follow this convention.\n"
                        f"- It ensures compatibility with common analysis tools, which typically assume depth increases. {notemark}",
                        indent=True,
                    )
                warn_list.append(
                    f"For {third_dimension} the positive direction is set to '{variable.positive}'"
                )
            var_min = variable[:].min()
            var_max = variable[:].max()
            var_units = variable.units if "units" in variable.ncattrs() else "No units"
            var_dims = list(variable.dimensions)
            var_summary = {
                "Range": f"{var_min:.3f} to {var_max:.3f}",
                "Units": var_units,
                "Dimensions": ", ".join(var_dims),
            }
            if missing_attributes:
                var_summary["Missing Attributes"] = missing_attributes

            # Determine expected order of dimensions based on coordinates and third dimension
            if var_type == "model" and len(variable.dimensions) == 3:
                # Ensure the third dimension is not already in primary_coords
                for coord in primary_coords:
                    if coord != third_dimension:
                        expected_order = preferred_dimensions_order[coord]
                        break

                expected_order = [third_dimension] + expected_order

                if var_dims != expected_order:
                    output(
                        f"Error: Variable '{var_name}' has incorrect dimension order: {', '.join(var_dims)}. Expected order is '{', '.join(expected_order)}'. {failmark}",
                        indent=True,
                    )
                    metadata_summary[var_name] = (
                        f"Incorrect dimension order: {', '.join(var_dims)} {failmark}"
                    )
                else:
                    message = f"Variable '{var_name}' has correct dimension order: {', '.join(var_dims)}."
                    output(
                        f"{message} {checkmark}",
                        indent=True,
                    )
                    passed_list.append(message)

            # Classify variable
            if var_type == "coordinates" or (var_type == "depth" and third_dimension):
                coordinate_variables[var_name] = var_summary
            elif var_type == "auxiliary" and var_type != "depth":
                coordinate_variables[var_name] = (
                    var_summary  # Auxiliary coordinates are also part of coordinates
                )
            else:
                model_variables[var_name] = var_summary

        output("\nCoordinate Variables:", indent=True)
        for var_name, summary in coordinate_variables.items():
            message = f"{var_name}: Range = {summary['Range']}, Units = {summary['Units']}, Dimensions = {summary['Dimensions']}"
            output(
                f"- {message} {checkmark}",
                indent=True,
                level=2,
            )
            passed_list.append(message)
            if "Missing Attributes" in summary:
                output(
                    f"Missing Attributes: {', '.join(summary['Missing Attributes'])} {failmark}",
                    indent=True,
                    level=3,
                )

        output("\nModel Variables:", indent=True)
        for var_name, summary in model_variables.items():
            mark = checkmark
            if var_name not in prop.valid_variable_list:
                mark = failmark
            message = f"{var_name}: Range = {summary['Range']}, Units = {summary['Units']}, Dimensions = {summary['Dimensions']}"
            output(
                f"- {message} {mark}",
                indent=True,
                level=2,
            )
            if mark == checkmark:
                passed_list.append(message)

            if mark == failmark:
                message = f"Invalid model variable name: '{var_name}' (see [2] & [4])"
                output(
                    message,
                    indent=True,
                    level=3,
                )
                failed_list.append(message)

            if "Missing Attributes" in summary:
                message = f"Missing Attributes: {', '.join(summary['Missing Attributes'])} {failmark}"
                output(
                    message,
                    indent=True,
                    level=3,
                )
                failed_list.append(
                    f"Missing Attributes: {', '.join(summary['Missing Attributes'])} (see [2] & [4])"
                )

        metadata_summary["Coordinate Variables"] = coordinate_variables
        metadata_summary["Model Variables"] = model_variables
    except Exception as e:
        output(f"Error listing variables: {e} {failmark}")
        metadata_summary["Variables"] = f"Error listing variables: {e} {failmark}"

    output(
        f"\n\nNOTES:\n    [1] Standards and Conventions see: {STANDARDS}\n    "
        f"[2] Model Variable Naming and Unit Requirements see: {VAR_STANDARDS}\n    "
        f"[3] Global Metadata File Structure see : {GLOBAL_STANDARDS}\n    "
        f"[4] Variable’s Metadata File Structure see : {VAR_STANDARDS}\n    "
        f"\n\n"
    )


def print_metadata_summary():
    try:
        for key, value in metadata_summary.items():
            if value is not None:
                if isinstance(value, dict):
                    output(f"{key}:", indent=True)
                    for sub_key, sub_value in value.items():
                        if isinstance(sub_value, str) and failmark in sub_value:
                            output(
                                f"- {sub_key}: {sub_value} {failmark}",
                                indent=True,
                                level=2,
                            )
                        else:
                            output(f"- {sub_key}: {sub_value}", indent=True, level=2)
                else:
                    if isinstance(value, str) and failmark in value:
                        output(f"{key}: {value} {failmark}", indent=True)
                    else:
                        output(f"{key}: {value}", indent=True)
    except Exception as e:
        output(f"Error printing metadata summary: {e} {failmark}")


def output(text, indent=False, level=1):
    indentation = "  " * (level - 1) if indent else ""
    print(f"{indentation}{text}")


def output_header(header):
    print(f"\n=== {header} ===\n")


if __name__ == "__main__":
    # Prompt the user to enter the file name
    file_name = input("Please enter the file name: ").strip()
    check_netcdf_file(file_name)
