**********************************************************
Fortnet: A Behler-Parrinello-Neural-Network Implementation
**********************************************************

|lgpl badge|

Fortnet is a Behler-Parrinello-Neural-Network implementation, written in modern
Fortran. Using atom-centered symmetry functions to characterize local atomic
environments, Fortnet provides easy access to the famous BPNN neural network
architecture to predict atomic or global properties of your physical system,
featuring powerful but optional MPI parallelism.

|Fortnet logo|


Installation
============

Building from source
--------------------

**Note:** This section describes the building with default settings in a typical
Linux environment. For more detailed information on the build customization and
the build process, consult the **detailed building instructions** in
`INSTALL.rst <INSTALL.rst>`_.

Download the latest stable source code from `GitHub
<https://github.com/vanderhe/fortnet/>`_::

  git clone https://github.com/vanderhe/fortnet.git

You need CMake (>= 3.16) to build Fortnet. If your environment offers no CMake
or only an older one, you can easily install the latest CMake via Python's
``pip`` command::

  pip install cmake

Start CMake by passing your compiler as environment variable ``FC`` and the
location where the code should be installed and the build directory
(``_build``) as options::

  FC=mpifort cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/fnet -B _build .

If the configuration was successful, start the build with::

  cmake --build _build -- -j

After successful build, you should test the code by running::

  pushd _build
  ctest -j
  popd

If the tests were successful, install the package via::

  cmake --install _build

For further details see the `detailed building instructions <INSTALL.rst>`_.


Documentation
=============

Consult following resources for documentation:

* `Step-by-step instructions with selected examples (Fortnet Recipes)
  <https://fortnet.readthedocs.io/>`_


Contributing
============

New features, bug fixes, documentation, tutorial examples and code testing is
welcome during the ongoing Fortnet development!

The project is `hosted on github <https://github.com/vanderhe/fortnet/>`_.
Please check `CONTRIBUTING.rst <CONTRIBUTING.rst>`_ for guide lines.

I am looking forward to your pull request!


License
=======

Fortnet is released under the GNU Lesser General Public License. See the included
`LICENSE <LICENSE>`_ file for the detailed licensing conditions.



.. |Fortnet logo| image:: ./utils/art/logo.svg
    :alt: Fortnet logo
    :width: 90
    :target: https://github.com/vanderhe/fortnet/

.. |lgpl badge| image:: ./utils/art/gnu-lgplv3.svg
    :alt: LGPL v3.0
    :scale: 100%
    :target: https://opensource.org/licenses/LGPL-3.0
