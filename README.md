# MCtools
Python scripts for generation of .egsinp and .egsphant files for BEAM/DOSXYZnrc from DICOM data.

**NOTE: These scripts have only been used with a Varian 21iX beam model and are provided with no support or guarantee that they will work out of the box. Tests should be run to verify the geometry of the resulting `.egsinp` files with your beam model.**

## Installation

The tools are not currently packaged for proper python installation. The tools can be imported in a python shell in the directory they are located by

```python
import dicomparse
import ctcreate
```

## Templates

Example templates for BEAM/DOSXYZnrc input files are provided in the `/templates` folder. The BEAMnrc model is only an example that resembles the basic layout of a Varian 21iX linear accelerator, and does not represent the actual specifications of a clinical linear accelerator.

The template examples show where certain variables (preceeded by `!!`, e.g. `!!X` representing the position of the X jaw) need to be placed in the template files.

Templates should be installed in a folder at `EGS_HOME/templates`

## BEAM/DOSXYZnrc Input File Generation

`.egsinp` files for BEAMnrc and DOSXYZnrc can be generated using the function

```
dicomparse.egsinpGenerate(plandir, EGS_HOME, beamind)
```

where `plandir` is the directory containing the DICOM RP file, `EGS_HOME` is the location of the `EGS_HOME` directory on the disk, and `beamind` is a list of beam indexes to generate files for (not providing `beamind` results in input files being generated for all beams).

Resulting `.egsinp` files are found in the `plandir` directory.

### Wedges

Testing has only been performed on wedges using Varian Enhanced Dynamic Wedges (EDWs). `EDW.py` is a script used for generating dynamic jaw input files using Varian Golden Segmented Treatment Tables.

## `.egsphant` Generation

`.egsphant` files can be generated using the function:

```
ctcreate.ctcreate(EGS_HOME, ctdir, structs={})
```

where `EGS_HOME` is as explained previously, `ctdir` is the directory containing the DICOM CT, dose (RD) and structure set (RS) files, and `structs` is a dictionary containing the names and corresponding physical density of any structures to override.

The function uses CT number to physical density information read from text files located in `EGS_HOME/templates`. Text files should be named as the serial number of the CT scanner to which the conversion table corresponds. An example is provided.