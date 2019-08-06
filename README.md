[![Anaconda-Server Badge](https://anaconda.org/conda-forge/foam_ct_phantom/badges/version.svg)](https://anaconda.org/conda-forge/foam_ct_phantom) [![Anaconda-Server Badge](https://anaconda.org/conda-forge/foam_ct_phantom/badges/latest_release_date.svg)](https://anaconda.org/conda-forge/foam_ct_phantom) [![Anaconda-Server Badge](https://anaconda.org/conda-forge/foam_ct_phantom/badges/platforms.svg)](https://anaconda.org/conda-forge/foam_ct_phantom) [![Anaconda-Server Badge](https://anaconda.org/conda-forge/foam_ct_phantom/badges/license.svg)](https://anaconda.org/conda-forge/foam_ct_phantom) [![Anaconda-Server Badge](https://anaconda.org/conda-forge/foam_ct_phantom/badges/downloads.svg)](https://anaconda.org/conda-forge/foam_ct_phantom)

[![Build Status](https://travis-ci.com/dmpelt/foam_ct_phantom.svg?branch=master)](https://travis-ci.com/dmpelt/foam_ct_phantom) [![Build status](https://ci.appveyor.com/api/projects/status/qp6ceia7iu05v7kr/branch/master?svg=true)](https://ci.appveyor.com/project/dmpelt/foam-ct-phantom/branch/master)


# foam_ct_phantom


A Python package for generating foam-like phantoms for CT.

* [\[Latest Release\]](https://github.com/dmpelt/foam_ct_phantom/releases/latest)
* [\[Version history\]](https://github.com/dmpelt/foam_ct_phantom/blob/master/CHANGELOG.md)
* [\[Bug Tracker\]](https://github.com/dmpelt/foam_ct_phantom/issues)
* [\[Documentation\]](https://dmpelt.github.io/foam_ct_phantom/)

A publication about this code is currently in preparation.

Development of this code is supported by Centrum Wiskunde & Informatica (CWI), with financial support provided by The Netherlands Organisation for Scientific Research (NWO), project number 016.Veni.192.235.

# Installation

To install this code in a conda environment, run:

```bash
conda install -c conda-forge foam_ct_phantom
```

In other environments, the code can be installed by running:

```bash
python setup.py install
```

The code requires the following Python modules: numpy, scipy, tqdm, h5py, psutil, numba, sortedcollections, pathlib2.
For compiling the code, the scikit-build module is required.

To run on GPU (only for faster cone-beam projection generation), a CUDA-capable GPU must be present and CUDA drivers must be installed. In addition, please make
sure that the version of the cudatoolkit package installed by conda matches the CUDA version of your drivers. Specific versions
of cudatoolkit can be installed by running (where 'X.X' is the CUDA version, e.g. '10.0'):

```bash
conda install cudatoolkit=X.X
```
