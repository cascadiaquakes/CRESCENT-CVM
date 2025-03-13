# Cascadia Community Velocity Model Tools (CVM Tools)

**Last updated:** 2025-03-13 **Release:** r0.9.7

This repository contains Python 3 tools developed to support the [Cascadia Community Velocity Models (CVM)](https://cascadiaquakes.org/cvm/) project. These tools are designed to facilitate the storage, extraction, and visualization of CVMs.

## Special Notes

- **Release r0.9.7:** Enhanced support for depth-positive convention (Issue #10).
- **Release r0.9.6:** Now using a standard variable naming convention and cvm_writer can decimate data.
- **Release r0.9.5:** Added a utility tool for converting 2D netCDF or GeoJSON files to CSV format.
- **Release r0.9.5:** Added support for 2D netCDF models and non-geospatial models.
- **Release r0.9:** The metadata attributes **grid_mapping** and **grid_mapping_name** have been replaced by **grid_ref**. All netCDF files must now include the **grid_ref** attribute.
- **Release r0.8:** Metadata template files have changed from Python data files to text files. The **cvm_write** tool no longer supports Python metadata files.
- **Release r0.8:** Metadata templates have been split into **global** and **variable** files to support both netCDF and HDF5 formats.

### Download

Clone the repository or download and unzip the files.

### Installation

Ensure you have the following:

- [Python](https://www.python.org/) 3
- Required Python modules (listed in `requirements.txt`):
  - Install with:
    ```
    pip install -r requirements.txt
    ```

Tested under Python 3.12.0 on macOS 14.2.1. Compatibility with older versions of Python 3 may vary.

## Documentation

For detailed documentation and tutorials, visit: [CVM-Tools Documentation](https://cascadiaquakes.github.io/cvm-tools-book/).

### Feedback or Issues

If you have any comments, questions, or bug reports, please open an issue ticket or contact [manochehr.bahavar@earthscope.org](mailto:manochehr.bahavar@earthscope.org).
