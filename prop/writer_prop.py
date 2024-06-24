# The GeoCSV version.
dataset_version = "GeoCSV 2.0"
conventions = {
    "Conventions": "CF-1.0",
    "Metadata_Conventions": "Unidata Dataset Discovery v1.0",
}

repository = {"repository_name": "EMC", "repository_institution": "EarthScope"}

json_indent = 4

global_prefix = "global_"

netcdf_format = {
    "NETCDF3C": {"format": "NETCDF3_CLASSIC", "encoding": None},
    "NETCDF4C": {
        "format": "NETCDF4_CLASSIC",
        "encoding": {"zlib": True, "complevel": 9},
    },
}

valid_output_types = ["metadata", "geocsv", "netcdf"]
valid_h5_output_types = ["metadata", "hdf5"]

read_chunk_size = 10000
