import os
import sys
import xarray as xr
import numpy as np


# Get the root _directory.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop
import writer_prop as writer_prop

LIB_DIR = os.path.join(ROOT_DIR, prop.LIB_DIR)
sys.path.append(LIB_DIR)
import shared_lib as lib
import writer_lib
from shared_lib import get_module_param, get_param, utc_now

logger = lib.get_logger()


def netcdf_to_geocsv(
    input_file,
    caller="cvm_convert.py",
):
    """Convert a geoCSV file content to GeoCSV.

    Keyword arguments:
    input_file -- [default None] The netCDF file to convert.
    caller -- [default cvm_convert.py] -- The calling script"""

    geocsv_metadata = list()
    geocsv_metadata.append(f"# dataset: {get_module_param(prop, 'dataset_version')}")
    geocsv_metadata.append(f"# delimiter: {prop.delimiter}")
    geocsv_metadata.append(f"# created: {utc_now()['date_time']} UTC ({caller})")
    geocsv_metadata.append(f"# netcdf_file: {os.path.basename(input_file)}")
    with xr.open_dataset(input_file, engine="netcdf4") as ds:
        df = ds.to_dataframe()

        for row in ds.attrs:
            # For some reason history lines are split (investigate).
            if row == "history":
                ds.attrs[row] = ds.attrs[row].replace("\n", ";")
            geocsv_metadata.append(f"# global_{row}: {ds.attrs[row]}")
        for var in ds.variables:
            if "variable" not in ds[var].attrs:
                geocsv_metadata.append(f"# {var}_variable: {var}")
                geocsv_metadata.append(f"# {var}_dimensions: {len(ds[var].dims)}")

            for att in ds[var].attrs:

                geocsv_metadata.append(f"# {var}_{att}: {ds[var].attrs[att]}")
                if att == "missing_value":
                    geocsv_metadata.append(f"# {var}_missing_value: {prop.na_rep}")
                if att == "variable":
                    geocsv_metadata.append(f"# {var}_dimensions: {len(ds[var].dims)}")
                    geocsv_metadata.append(f"# {var}_column: {var}")
    metadata = "\n".join(geocsv_metadata)
    data = df.to_csv(sep=prop.delimiter, na_rep=prop.na_rep)

    return metadata, data


# Convert dataset attributes to JSON, handling non-serializable values
def convert_to_json_safe(obj):
    if isinstance(obj, (np.integer, np.floating, np.bool_)):
        return obj.item()
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj


# Output metadata as JSON.
def json_metadata(input_file):
    with xr.open_dataset(input_file) as ds:
        metadata_json = {
            key: convert_to_json_safe(value) for key, value in ds.attrs.items()
        }

        file_size_bytes = os.path.getsize(input_file)
        file_size_kb = file_size_bytes / 1024
        metadata_dict = dict()
        metadata_dict["vars"] = list(ds.variables)

        if file_size_kb:
            metadata_dict["size_kb"] = round(file_size_kb, 2)
        metadata_dict["data_vars"] = list(ds.data_vars)
        z_var = list(
            set(metadata_dict["vars"])
            - set(metadata_dict["data_vars"])
            - set(["latitude", "longitude", "easting", "northing", "x", "y"])
        )
        metadata_dict["z_var"] = ""
        if len(z_var) >= 1:
            metadata_dict["z_var"] = z_var[0]

        # Include variable attributes
        variables_json = {}
        for var_name, var in ds.variables.items():
            var_attrs = {
                key: convert_to_json_safe(value) for key, value in var.attrs.items()
            }
            variables_json[var_name] = var_attrs
        metadata_dict["variables"] = variables_json

        # Combine dataset and variable attributes
        output_json = lib.merge_dictionaries(metadata_json, [metadata_dict])

    return output_json


# This function is not used and should be removed if no error occurs.
def geocsv_to_netcdf_old(
    input_file,
    caller="cvm_convert.py",
):
    """Convert a GeoCSV file content to netCDF.

    Keyword arguments:
    input_file -- [default None] The GeoCSV file to convert.
    caller -- [default cvm_convert.py] -- The calling script"""

    with open(input_file) as csvfile:

        metadata_dict = dict()
        var_dict = dict()
        data = csvfile.read()
        lines = data.split("\r\n")
        delimiter = None

        for line in lines:
            _key = line.strip()
            if not (_key):
                continue
            # Only the header part (metatdata).
            if not _key[0] == "#":
                # Watch for possible breaks in long history lines.
                if attr == "history":
                    logger.error(
                        f"[ERR] Unexpected line break. Please check the {attr} metadata value."
                    )
                    sys.exit(2)
                else:
                    break
            _key = _key[1:].strip()
            # Global
            if _key.startswith("global"):
                attr, value = _key.replace("global_", "").split(":", 1)
                attr = attr.strip()
                value = value.strip()
                metadata_dict[attr] = value
            # Variables.
            else:
                attrib, value = _key.split(":", 1)
                if attrib == "delimiter":
                    delimiter = value
                    continue
                if "_" not in attrib:
                    continue
                var, attr = attrib.strip().split("_", 1)

                if attr in ("column", "file"):
                    continue
                if var not in var_dict:
                    var_dict[var] = dict()
                var_dict[var][attr] = value.strip()
        metadata_dict["geocsv_file"] = os.path.basename(input_file)
    if delimiter is None:
        logger.error(
            f"[ERR] The geoCSV file {input_file} header does not define the required 'delimiter' parameter"
        )
        sys.exit(2)
    df = lib.read_csv(input_file, delimiter)
    # dataset = pd.read_csv(input_file, sep=delimiter, comment="#")
    ds = df.set_index(["depth", "latitude", "longitude"]).to_xarray()
    ds.attrs = metadata_dict
    for var in var_dict:
        ds[var].attrs = var_dict[var]
    ds.to_netcdf(
        "test.nc",
    )
    coords = writer_lib.get_coords(df, var_dict)

    return metadata_dict, df
