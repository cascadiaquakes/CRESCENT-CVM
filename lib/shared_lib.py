import os
import sys
import io
import logging
import time
from datetime import datetime, timezone
from pprint import pprint

from pyproj import Proj
import xarray as xr
import pandas as pd
from tqdm import tqdm

# Get the directory paths.
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

PROP_DIR = os.path.join(ROOT_DIR, "prop")
sys.path.append(PROP_DIR)
import shared_prop as prop
import writer_prop as writer_prop


def cf_coordinate_names(x,y, aux=False):
    """Produce CF standard coordinate variable names.

    Keyword arguments:
    x -- [required] the current x-coordinate variable name
    y -- [required] the current y-coordinate variable name
    aux -- [default=False] flag to indicate if these variables are for the auxiliary coordinates.
    """    
    if x.lower() == "longitude" and y.lower() == "latitude":
        return {"x": "longitude", "y":"latitude"}
    elif aux:
        return {"x": x, "y": y}
    elif x.lower() != "longitude" and y.lower() != "latitude":
        return {"x": "x", "y": "y"}
    else:
        logging.error(f"[cf_coordinate_names ERR] invalid coordinate name combinations ({x},{y})")
        raise


def get_geocsv_metadata(ds):
    """Compose GeoCSV style metadata for a given Dataset.

    Keyword arguments:
    ds -- [required] the xarray dataset
    """
    geocsv_metadata = list()
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
    metadata = f'{"\n".join(geocsv_metadata)}\n'
    return metadata


def view_list(lst, indent=14, chunk=20):
    """Print a chunk from start and end of a list to view.

    Keyword arguments:
    lst -- [required] list of the numbers
    chunk -- [default: 20] print chunk elements from start and end of the list.
    """
    file = io.StringIO()
    # Not needed, but move pointer to start of file.
    file.seek(0)
    if len(lst) <= 60:
        pprint(
            list(lst),
            file,
            width=80,
            compact=True,
            indent=indent,
        )
    else:
        pprint(
            list(lst[0:chunk]),
            file,
            width=80,
            compact=True,
            indent=indent,
        )
        print(f"{indent*' '}...", file=file)
        pprint(
            list(lst[-1 * chunk :]),
            file,
            width=80,
            compact=True,
            indent=indent,
        )
    print(file.getvalue().replace("]", " ").replace("[", " "))


def closest(lst, value):
    """Find the closest number in a list to the given value

    Keyword arguments:
    lst -- [required] list of the numbers
    value -- [required] value to find the closest list member for.
    """
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - value))]


def utc_now():
    """Return the current UTC time."""
    try:
        _utc = datetime.now(tz=timezone.utc)
        utc = {
            "date_time": _utc.strftime("%Y-%m-%dT%H:%M:%S"),
            "datetime": _utc,
            "epoch": _utc.timestamp(),
        }
        return utc
    except:
        logging.error(f"[UTC ERR] Failed to get the current UTC time")
        raise


def read_csv(
    data_file,
    delimiter,
):
    """Load the CSV data file and show the progress.

    Keyword arguments:
    data_file -- [required] input CSV data file
    delimiter -- [required] delimiter marker
    """
    delimiter = delimiter.strip()
    if len(delimiter) <= 0:
        # df = pandas.read_csv(data_file, delim_whitespace=True)
        df = pd.concat(
            [
                chunk
                for chunk in tqdm(
                    pd.read_csv(
                        data_file,
                        chunksize=writer_prop.read_chunk_size,
                        delim_whitespace=True,
                        comment="#",
                        skip_blank_lines=True,
                        header=0,
                    ),
                    desc="...loading data",
                )
            ]
        )
    else:
        # df = pandas.read_csv(data_file, delimiter=delimiter)
        df = pd.concat(
            [
                chunk
                for chunk in tqdm(
                    pd.read_csv(
                        data_file,
                        chunksize=writer_prop.read_chunk_size,
                        sep=delimiter,
                        comment="#",
                        skip_blank_lines=True,
                        header=0,
                    ),
                    desc="...loading data",
                )
            ]
        )
    return df


def check_file_type(file_path):
    """Check a given file's format

    Keyword arguments:
    file_path -- [default None] File to check."""
    file_type = dict()
    # Check for a valid GeoCSV file.
    try:
        with open(file_path, newline="") as csvfile:
            start = csvfile.read(300).strip()
            if start.startswith(f"# dataset:"):
                file_type["engine"] = "geocsv"
                version = start.split("\n")[0].split(":")[1].strip()
                if version == prop.dataset_version:
                    file_type["valid"] = True
                else:
                    file_type["valid"] = (
                        f"[ERR] Invalid {file_type['engine']} file type. For 'dataset', "
                        f"expected {prop.dataset_version}, but received {version}"
                    )
            else:
                file_type["engine"] = None
                file_type["valid"] = (
                    f"[ERR] For a GeoCSV file, the '# dataset: value' should always be present "
                    f"and should be the first line of a dataset."
                )
    except UnicodeDecodeError:
        # Check for a valid netCDF file.
        try:
            nc_dataset = xr.open_dataset(file_path, engine="netcdf4")
            file_type["engine"] = "netcdf"
            file_type["valid"] = True
        except Exception as ex:
            file_type["engine"] = None
            file_type["valid"] = (
                f"[ERR] dataset File {file_path} type not recognized.\n {ex}"
            )
    except Exception as ex:
        file_type["engine"] = None
        file_type["valid"] = f"[ERR] File {file_path} type not recognized.\n {ex}"

    return file_type


def get_logger(log_file=None, log_stream=None, config=True, log_level=0):
    """Set up the logging.

    Keyword arguments:
    log_file -- [default None] name of the log file to write to
    log_stream -- [default None] the log stream to log to
    config -- [default True] should configuring the logging module
    log_level -- [default 0] minimum priority level of messages to log (NOTSET=0, DEBUG=10,
                INFO=20, WARN=30, ERROR=40, CRITICAL=50)
    """
    root = logging.getLogger()
    logging.getLogger("matplotlib.font_manager").disabled = True
    if root.handlers:
        for handler in root.handlers:
            root.removeHandler(handler)
    if config or log_file is not None or log_stream is not None:
        if log_file is not None:
            logging.basicConfig(
                filename=log_file,
                format="%(message)s",
                encoding="utf-8",
                level=log_level,
            )
        elif log_stream is not None:
            logging.basicConfig(
                stream=log_stream,
                format="%(message)s",
                encoding="utf-8",
                level=log_level,
            )
        else:
            logging.basicConfig(format="%(message)s", encoding="utf-8", level=log_level)

    # Retrieve the logger instance
    logger = logging.getLogger()
    return logger


def get_param(params, var, required=True):
    """Extract a variable from a parameter module..

    Keyword arguments:
    param -- [required] the parameter module
    var -- [required] the variable to extract from the param object.
    required [default True] -- is var required? If True, the code will
    exit if not found.

    """

    # Set up the logger.
    logger = get_logger()

    if var not in dir(params):
        logger.error(
            f"[ERR] parameter '{var}' is required but it is missing from the parameter file."
        )
        if required:
            sys.exit(2)
        else:
            return None

    return getattr(params, var)


def merge_dictionaries(parent_dict, dict_list):
    """Merge a list of dictionaries to the parent dictionary using update() method.

    Keyword arguments:
    parent_dict -- [required] The parent dictionary.
    dict_list -- a list of dictionaries to merge with the parent dictionary.
    """
    for _dic in dict_list:
        parent_dict.update(_dic)

    return parent_dict


def insert_after(parent_dict, dict_item, after_key, after=True):
    """Insert an item into a dictionary before or after a given item.

    Keyword arguments:
    parent_dict -- [required] The parent dictionary
    dict_item -- [required] The dictionary item to insert
    after_key -- [required] The dictionary key to insert the item after
    after -- [default True] if True, will insert after the given key. Otherwise, insert before it
    """
    delta = 0
    if after:
        delta = 1

    pos = list(parent_dict.keys()).index(after_key) + delta
    items = list(parent_dict.items())
    items.insert(pos, dict_item)

    return dict(items)


def set_variable_keys(variable, var):
    updated_var = dict()
    for key in variable:
        updated_var[f"{var}_{key}"] = variable[key]


def project_lonlat_utm(
    longitude, latitude, utm_zone, ellipsoid, xy_to_latlon=False, preserve_units=False
):
    """
    Performs cartographic transformations. Converts from longitude, latitude to UTM x,y coordinates
    and vice versa using PROJ (https://proj.org).

     Keyword arguments:
    longitude (scalar or array) – Input longitude coordinate(s).
    latitude (scalar or array) – Input latitude coordinate(s).
    xy_to_latlon (bool, default=False) – If inverse is True the inverse transformation from x/y to lon/lat is performed.
    preserve_units (bool) – If false, will ensure +units=m.
    """
    P = Proj(
        proj="utm",
        zone=utm_zone,
        ellps=ellipsoid,
    )
    # preserve_units=preserve_units,

    x, y = P(
        longitude,
        latitude,
        inverse=xy_to_latlon,
    )
    return x, y


def time_it(t0):
    """Compute the elapsed time since last call

    Keyword arguments:
    t0 -- [required] Time of the last call
    """
    t1 = time.time()
    dt = t1 - t0
    if dt >= 60.0:
        text = f"Elapsed time: {dt / 60.0:0.1f}m"
    else:
        text = f"Elapsed time: {dt:0.2f}s"

    return time.time(), text
