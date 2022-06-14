.. _sec-interfaces-sockets:

####################
Socket-Communication
####################

For calculations that heavily rely on file-IO and situations where
initialization times are not negligible, socket communcation is often preferable
and might lead to a significant reduction of wallclock time. The i-PI protocol
:cite:`ipi` is among the most famous options in this regard and supported by
Fortnet.

Below, compiling and running Fortnet with support for socket communication is
outlined.

.. note::

   Note: At this time, only communication with ASE or i-PI has been tested. In
   general, however, other software which correctly implements the i-PI protocol
   should be supported.

Compiling with Socket Support
=============================

To enable socket-communication in Fortnet the ``WITH_SOCKETS`` flag in the
configuration file ``config.cmake`` must be set to ``TRUE``, before starting the
compilation process. Equivalently, this configuration step is also possible
directly from the command line, e.g. invoke:
::

   FC=mpifort CC=gcc cmake -DTOOLCHAIN=gnu -DWITH_SOCKETS=1 -B _build .

.. _sec_interfaces-sockets-input:

Providing the Input for Fortnet
===============================

[Input: `recipes/interfaces/sockets/input/`]

To establish a connection via socket communication, Fortnet is called with the
appropriate ``Socket{}`` driver option, e.g.:
::

   Driver = Socket {
     File = 'fortnet'
     Protocol = i-PI {
       # Non-periodic water molecule (H2O)
       Periodic = No
       AtomicNumbers = 8 1 1
     }
     MaxSteps = 1000
     Verbosity = 0
   }

It instructs Fortnet to read geometry driving commands via a named temporary
communication file (stored in /tmp/). The code is asked to perform up to
``MaxSteps = 1000`` geometry steps under control of the external driver. For
codes that support the `EXIT` command in the protocol, Fortnet may be told to
continue running until told to stop by using::

  MaxSteps = -1

For initialization purposes, rudimentary geometry specifications, like the
(atomic) number of atoms and boundary conditions, must be provided. Make sure
that the atomic numbers are given in the same order as for the geometries that
will be send via the socket connection later on.

Additionally, a network initialization file (``NetstatFile``) must be specified
and read. Operating in socket mode always requires force analysis, which is
assumed (and overwritten) even if not requested by the user:
::

   Data {
     NetstatFile = 'fortnet.hdf5'
   }

   Options {
     ReadNetStats = Yes
     Mode = 'predict'
   }

   Analysis {
     Forces = Analytical {}
   }

For the full input file assembled from the code snippets shown here, please
consult the archive, whose location is stated at the begin of this section.
