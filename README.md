# Cascadia Community Velocity Model Tools (CVM Tools)

[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)

Python 3 tools for storing, extracting and visualizing the Cascadia Community Velocity Models ([CVM](https://cascadiaquakes.org/cvm/)).

**Last updated:** 2025-03-13 **Release:** r0.9.7

**Documentation:** https://cascadiaquakes.github.io/cvm-tools-book/

**Viewer:** https://cvm.cascadiaquakes.org

---

## Release Notes

- **Release r0.9.7:** Enhanced support for depth-positive convention (Issue #10).
- **Release r0.9.6:** Now using a standard variable naming convention and cvm_writer can decimate data.
- **Release r0.9.5:** Added a utility tool for converting 2D netCDF or GeoJSON files to CSV format.
- **Release r0.9.5:** Added support for 2D netCDF models and non-geospatial models.
- **Release r0.9:** The metadata attributes **grid_mapping** and **grid_mapping_name** have been replaced by **grid_ref**. All netCDF files must now include the **grid_ref** attribute.
- **Release r0.8:** Metadata template files have changed from Python data files to text files. The **cvm_write** tool no longer supports Python metadata files.
- **Release r0.8:** Metadata templates have been split into **global** and **variable** files to support both netCDF and HDF5 formats.

---

## Installation

Clone the repository or download and unzip the files, then make sure you have:

- [Python](https://www.python.org/) 3
- Required Python modules (listed in `requirements.txt`):
  - Install with:
    ```
    pip install -r requirements.txt
    ```

Tested under Python 3.12.0 on macOS 14.2.1. Compatibility with older versions of Python 3 may vary.

---

## License

This project is licensed under the BSD 3-Clause License. See [LICENSE](LICENSE).

---

## Contacts

- **Manochehr Bahavar** (original author), EarthScope: [manochehr.bahavar@earthscope.org](mailto:manochehr.bahavar@earthscope.org)
- **Amanda M. Thomas**, University of California, Davis: [amthom@ucdavis.edu](mailto:amthom@ucdavis.edu)
- **William Marfo**, University of California, Davis: [wmarfo@ucdavis.edu](mailto:wmarfo@ucdavis.edu)

For bugs or questions, please [open an issue](https://github.com/cascadiaquakes/CRESCENT-CVM/issues).

---

## Acknowledgments

<table>
  <tr>
    <td align="center" valign="middle">
      <img src="https://raw.githubusercontent.com/cascadiaquakes/.github/main/assets/crescent-circle.png" alt="CRESCENT Logo" width="150" />
    </td>
    <td align="center" valign="middle">
      <img src="https://raw.githubusercontent.com/cascadiaquakes/.github/main/assets/nsf-logo.png" alt="NSF Logo" width="450" />
    </td>
  </tr>
</table>

This project is supported by the National Science Foundation through the Cascadia Region Earthquake Science Center (CRESCENT), Cooperative Agreement EAR-2225286.
