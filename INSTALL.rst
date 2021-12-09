*******************************
Building and installing Fortnet
*******************************

If you have problems with the build, you can find suggestions for some
frequently occuring scenarios in the `Troubleshooting <#troubleshooting>`_
section.


Requirements
============

In order to compile Fortnet, you need the following software components:

* A Fortran 2003 compliant compiler

* Compatible serial HDF5 with Fortran API and High-Level routines
  (version 1.10.x or newer)

* GNU make

* CMake (version 3.16 or newer)

* Python (version >= 3.2) for the source preprocessor


Requirements for testing Fortnet
--------------------------------

In order to execute the code tests and validate them against precalculated
results, you will additionally need:

* Python (version >= 3.2) with numPy and h5py


Tested build environments
-------------------------

Fortnet is built and tested for both serial and MPI environments on the
following architectures:

+---------------+------------------------+---------------+-------------------+-----+
| Architecture  | Compiler               | MPI           | Serial HDF5       |Notes|
+===============+========================+===============+===================+=====+
| x86_64 /      | | GNU Fortran 8.3.0    | | OpenMPI 4.0 | | 1.10.4          |  /  |
| Linux         | | GNU Fortran 8.4.0    | | OpenMPI 4.0 | | 1.10.4          |     |
+---------------+------------------------+---------------+-------------------+-----+
| x86_64 /      | GNU Fortran  9.3.0     | OpenMPI 4.0   | 1.10.4            |  /  |
| Linux         |                        |               |                   |     |
+---------------+------------------------+---------------+-------------------+-----+
| x86_64 /      | | GNU Fortran 10.1.0   | | OpenMPI 4.0 | | 1.10.4          |  /  |
| Linux         | | GNU Fortran 10.2.0   | | OpenMPI 4.0 | | 1.10.4          |     |
+---------------+------------------------+---------------+-------------------+-----+
| x86_64 /      | GNU Fortran 11.1.0     | OpenMPI 4.0   | 1.10.4            |  /  |
| Linux         |                        |               |                   |     |
+---------------+------------------------+---------------+-------------------+-----+
| x86_64 /      | Intel Fortran 18.0.5   | MPICH 3.2     | 1.12.1 / 1.13.0   |  /  |
| Linux         |                        |               |                   |     |
+---------------+------------------------+---------------+-------------------+-----+
| x86_64 /      | Intel Fortran 19.0.5   | MPICH 3.3     | 1.12.1 / 1.13.0   |  /  |
| Linux         |                        |               |                   |     |
+---------------+------------------------+---------------+-------------------+-----+
| x86_64 /      | | Intel Fortran 20.2.1 | | IMPI 2021.1 | | 1.12.1 / 1.13.0 |  /  |
| Linux         | | Intel Fortran 20.2.3 | | IMPI 2021.1 | | 1.12.1 / 1.13.0 |     |
+---------------+------------------------+---------------+-------------------+-----+


Obtaining the source
====================

Clone the `public git repository <https://github.com/vanderhe/fortnet>`_.
The tagged revisions correspond to stable releases, while the default branch
contains the latest development version. ::

  git clone https://github.com/vanderhe/fortnet.git
  cd fortnet

The project uses git-submodules for some external dependencies, which will be
automatically retrieved during configuration.


Building
========

**Important note:** CMake caches its variables in the `CMakeCache.txt` file in
the build folder (e.g. ``_build/CMakeCache.txt``). Make sure to delete this file
before re-running CMake if you have changed variables in `config.cmake` or in
the toolchain files in the `sys/` folder. (Deleting the `CMakeCache.txt` file is
not necessary if you change a variable via the ``-D`` command line option.)

In order to build Fortnet carry out the following steps:

* | Inspect the `config.cmake` file and customise the global build parameters.
  | (If you are unsure, leave the defaults as they are.)

* Invoke CMake to configure the build. Specify the installation destination
  (e.g. ``$HOME/opt/fnet``) and pass an arbitrary folder (e.g. ``_build``) for
  the build and the directory containing the source files (e.g. ``.``) as
  arguments to CMake. Additionally define your Fortran compiler as
  environment variables, e.g. (in a BASH compatible shell)::

    FC=mpifort cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/fnet -B _build .

  Based on the detected compilers, the build system will read further settings
  from a corresponding toolchain file in the `sys/` folder. Either from a
  compiler specific one (e.g. `gnu.cmake`, `intel.cmake`, etc.) or the generic
  one (`generic.cmake`) if the detected compiler combination does not correspond
  to any of the specific settings. The selected toolchain is indicated in the
  CMake output. (The toolchain file selection can be manually overriden by
  setting the ``TOOLCHAIN`` CMake variable.)

  You may adjust any CMake variable defined in `config.make` or in the
  toolchain files by either modifying the files directly or by setting
  (overriding) the variable via the ``-D`` command line option. For example, in
  order to disable MPI parallelism, you would have to override the ``WITH_MPI``
  variable with the CMake command line argument ``-D``::

    -DWITH_MPI=0

* If the configuration was successful, start the build by::

    cmake --build _build -- -j

  This will compile the code using several threads and showing only the most
  relevant information.

  If, for debugging purposes, you wish to see the exact compiling commands, you
  should execute a serial build with verbosity turned on instead::

    cmake --build _build -- VERBOSE=1

* Note: By default the code is compiled with distributed memory parallelism
  (MPI) enabled. In case you want to disable it, override the corresponding
  variable ``WITH_MPI`` as shown above.


Testing Fortnet
===============

* After successful compilation, change to the build folder and execute the code
  tests::

    pushd _build
    ctest
    popd

  You can also run the tests in parallel in order to speed this up. If you use
  parallel testing, ensure that the number of MPI-processes is reduced
  accordingly. As an example, assuming your workstation has 4 cores and you have
  set up the ``TEST_MPI_PROCS`` variable to ``2`` (in `config.cmake`), issue ::

    ctest -j2

  for an MPI enabled binary running two tests simultaneously, each using 2
  cores.

  The ``TEST_MPI_PROCS`` cache variable can be updated or changed also after
  the compilation by invoking CMake with the appropriate ``-D`` options, e.g.::

    cmake -B _build -DTEST_MPI_PROCS=2 .
    pushd _build
    ctest
    popd


Installing Fortnet
==================

* The compiled executables, libraries, module files etc. can be copied into an
  installation directory by ::

    cmake --install _build

  where the destination directory can be configured by the variable
  ``CMAKE_INSTALL_PREFIX`` (in the `config.cmake` file). The default location is
  the `_install` subdirectory within the build directory.


Generating developer documentation
==================================

Developer documentation can be generated using the FORD source code
documentation generator by issuing ::

  cd doc/fortnet/ford && ford fortnet.md

in the main source directory. The documentation will be created in the
`doc/fortnet/ford/doc` folder.


Developer build instructions
============================

You should avoid customizing the build by directly changing variables in the
CMake config files, as your changes may accidently be checked in into the
repository. Instead, create a customized CMake config file, where you
pre-populate the appropriate cache variables. Then use the `-C` option to load
that file::

  FC=mpifort cmake -C custom.cmake -B _build .

The customized config file is read by CMake before the compiler detection
stage.


Advanced build configuration
============================

Controlling the toolchain file selection
----------------------------------------

You can override the toolchain file, and select a different provided case,
passing the ``-DTOOLCHAIN`` option with the relevant name, e.g. ::

  -DTOOLCHAIN=gnu

or by setting the toolchain name in the ``FNET_TOOLCHAIN`` environment
variable selects it. If you want to load an external toolchain file instead of
the bundled ones, you can specify the file path with the ``-DTOOLCHAIN_FILE``
option ::

  -DTOOLCHAIN_FILE=/some/path/myintel.cmake

or with the ``FNET_TOOLCHAIN_FILE`` environment variable.

Similarly, you can also use an alternative build config file instead of
`config.cmake` by specifying it with the ``-DBUILD_CONFIG_FILE`` option or by
defining the ``FNET_BUILD_CONFIG_FILE`` environment variable.


Preventing the download of external sources
-------------------------------------------

Depending on the value of the ``HYBRID_CONFIG_METHODS`` configuration variable,
some dependencies (e.g. mpifx) are automatically downloaded during the
configuration phase and built during the Fortnet build process. If you want to
ensure that nothing gets downloaded during the build, pass the variable
definition ::

  -DHYBRID_CONFIG_METHODS="Find"

to CMake during the configuration. In this case, CMake will only try to find
those dependencies on the system (by searching in the standard system paths and
in the locations defined in the environment variable ``CMAKE_PREFIX_PATH``) and
stop if some components were not found.


Troubleshooting
===============

* **CMake finds the wrong compiler**

  CMake should be guided with the help of the environment variable ``FC`` to
  make sure it uses the right compilers, e.g. ::

    FC=mpifort cmake [...]

* **CMake does not find HDF5**

  You have to make sure that an HDF5 installation is present, that matches your
  compiler. The rudimentary steps to compile HDF5 from source could look similar
  to this (assumind you already installed an Intel compiler)::

    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.1/src/CMake-hdf5-1.12.1.tar.gz
    tar xfz CMake-hdf5-1.12.1.tar.gz
    cd CMake-hdf5-1.12.1/hdf5-1.12.1/
    CC=icc CXX=icpc FC=ifort F9X=ifort ./configure --prefix=${PWD}/hdf5 --enable-fortran --with-default-api-version=v110 --enable-shared
    make -j -l2
    make install

  For self-compiled HDF5 instances, CMake should be guided with the help of the
  environment variable ``HDF5_ROOT`` to make sure it searches at the right
  location, e.g. ::

    export HDF5_ROOT=/home/user/CMake-hdf5-1.12.1/hdf5-1.12.1/hdf5
