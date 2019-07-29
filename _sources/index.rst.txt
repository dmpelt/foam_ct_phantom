.. FoamCTPhantom documentation master file, created by
   sphinx-quickstart on Mon Jul 29 11:34:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FoamCTPhantom's documentation!
=========================================

A Python package for generating foam-like phantoms for CT.

* `[Latest Release] <https://github.com/dmpelt/foam_ct_phantom/releases/latest>`_
* `[Version history] <https://github.com/dmpelt/foam_ct_phantom/blob/master/CHANGELOG.md>`_
* `[Bug Tracker] <https://github.com/dmpelt/foam_ct_phantom/issues>`_
* `[Documentation] <https://dmpelt.github.io/foam_ct_phantom/>`_

A publication about this code is currently in preparation.

Development of this code is supported by Centrum Wiskunde & Informatica (CWI), with financial support provided by The Netherlands Organisation for Scientific Research (NWO), project number 016.Veni.192.235.

Installation
-------------

To install this code in a conda environment, run:

.. code-block:: bash


    conda install -c conda-forge foam_ct_phantom

In other environments, the code can be installed by running:

.. code-block:: bash


    python setup.py install

The code requires the following Python modules: numpy, scipy, tqdm, h5py, psutil, numba, sortedcollections, pathlib2.
For compiling the code, the scikit-build module is required.

To run on GPU (only for faster cone-beam projection generation), a CUDA-capable GPU must be present and CUDA drivers must be installed. In addition, please make
sure that the version of the cudatoolkit package installed by conda matches the CUDA version of your drivers. Specific versions
of cudatoolkit can be installed by running (where 'X.X' is the CUDA version, e.g. '10.0'):

.. code-block:: bash


    conda install cudatoolkit=X.X


Usage
-----

Please see the example scripts for usage information.

.. toctree::
   :maxdepth: 2

   auto_examples/index
   apiref/modules

* :ref:`genindex`
* :ref:`modindex`
