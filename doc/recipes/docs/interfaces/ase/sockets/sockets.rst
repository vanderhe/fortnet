.. highlight:: none
.. _sec-interfaces-ase-sockets:

********************
Socket-Communication
********************

[Input: `recipes/interfaces/ase/sockets/`]

For calculations that heavily rely on file-IO, it is better to use communication
between Fortnet and ASE with the i-PI protocol. In these cases, the reduction of
the wallclock time can be significant.

Below, an exemplary geometry optimization, controlled by ASE via socket
communication, is performed. A water molecule is used for demonstration.

.. _sec_interfaces_ase_sockets-geoopt-input:

Geometry Optimization by ASE
============================

Providing the input for Fortnet
-------------------------------

To establish a connection via socket communication, Fortnet is called with the
appropriate ``Socket{}`` driver option, specifying, among others, the expected
atom elements and boundary conditions:
::

   Driver = Socket {
     File = 'fortnet'
     Protocol = i-PI {
       # Non-periodic water molecule (H2O)
       Periodic = No
       AtomicNumbers = 8 1 1
     }
     MaxSteps = -1
     Verbosity = 0
   }

It instructs Fortnet to read geometry driving commands via a named temporary
communication file (stored in /tmp/). The code is asked to continue running
until told to stop (``MaxSteps = -1``).

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

Main Script
-----------
The geometry used is a simple, unoptimized water molecule, whose file location
is specified in the main script (`GEO_IN`)::

   3

   O     0.0000000000E+00  -1.0000000000E+00   0.0000000000E+00
   H     0.0000000000E+00   0.0000000000E+00   7.8306400000E-01
   H     0.0000000000E+00   0.0000000000E+00  -7.8306400000E-01

Following the necessary imports, the main script defines the type of socket to
be used (UNIX socket) as well as the path to the Fortnet executable
(`FNET`) and the geometry (`GEO_IN`). The main method then first reads in
the geometry to optimize. Subsequently, the ASE driver ``BFGS()`` for geometry
optimization is specified and where to write the trajectory and the logfile.

Finally, the last trajectory image that contains the optimized water molecule is
extracted and the geometry written to disk:

.. code-block:: python

   import sys
   from subprocess import Popen
   from ase.io import read, write
   from ase.optimize import BFGS
   from ase.io.trajectory import Trajectory
   from ase.calculators.socketio import SocketIOCalculator


   UNIXSOCKET = 'fortnet'
   FNET = 'fnet'
   GEO_IN = './H2O_unoptimized.xyz'
   GEO_OUT = './H2O_optimized.xyz'
   TRAJECTORY = 'opt.traj'


   def main():
       '''Main driver routine.'''

       system = read(GEO_IN, format='xyz')

       opt = BFGS(system, trajectory=TRAJECTORY, logfile='opt.log')

       with SocketIOCalculator(log=sys.stdout, unixsocket=UNIXSOCKET) as calc:
	   Popen(FNET)
	   system.set_calculator(calc)
	   opt.run(fmax=1E-04)

       traj = Trajectory(TRAJECTORY)
       atoms = traj[-1]
       write(GEO_OUT, atoms, format='xyz')

   if __name__ == "__main__":
       main()

You may also consider inspecting the trajectory of the run by invoking:
::

  ase gui opt.traj

.. note::

    To correctly close sockets on the ASE side, call `calc.close()` at the end
    or, more elegantly, enclose the class ``SocketIOCalculator`` using the
    `with` statement as done in the example shown here. Nevertheless, in the
    current state of ASE, the socket gets closed without warning missing the
    'EXIT' string of the i-PI protocol, which always leads to an error message
    issued by Fortnet at the end of a calculation driven by
    socket-communication.
