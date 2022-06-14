#!/usr/bin/env python3

'''
Runs regression test for socket communication.
'''


import struct
import socket
import os
import numpy as np
import numpy.linalg as la
from sockettools import readgen, receive_all, a0

# Expecting two geometry steps of communication with Fortnet
NR_STEPS = 2


def connect():

    pid = os.getpid()

    server_address = '/tmp/ipi_fnet%i' % pid

    # Make sure the socket does not already exist
    try:
        os.unlink(server_address)
    except OSError:
        if os.path.exists(server_address):
            raise

    # write file for fortnet_in.hsd to include:
    with open("file.hsd", "w") as fd:
        fd.write('# The externally set filename for this run\n')
        fd.write("+Driver = +Socket {\n")
        fd.write('  !File = "fnet%i"\n' % pid)
        fd.write("}\n")

    # plain text file with the same information
    with open("file.txt", "w") as fd:
        fd.write("/tmp/ipi_fnet%i" % pid)

    serversocket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    serversocket.bind(server_address)
    # become a server socket with at most 1 connection
    serversocket.listen(1)
    connection, _ = serversocket.accept()

    return connection


def main():

    specienames, species, coords, origin, latvecs = readgen("H2O_ice.gen")
    # convert to atomic units (native i-PI unit system)
    coords /= a0
    latvecs /= a0

    connection = connect()

    for istep in range(NR_STEPS):
        print("Step %i" % istep)
        connection.sendall(b'POSDATA     ')
        connection.sendall(latvecs)

        connection.sendall(la.inv(latvecs))
        connection.sendall(np.int32(len(coords)))
        connection.sendall(coords)

        connection.sendall(b'GETFORCE    ')
        # needs work:
        buf = receive_all(connection, 12)
        if buf != b'FORCEREADY  ':
            raise ValueError('Unexpected value of "GETFORCE": "%s"!' % buf)

        # expecting energy and number of atoms
        buf = receive_all(connection, 12)
        unpacked_data = struct.unpack('di', buf)
        print(unpacked_data)

        # forces
        buf = receive_all(connection, len(coords)*3*8)
        frmat = '%i' % (3*len(coords))
        frmat += 'd'
        unpacked_data = struct.unpack(frmat, buf)
        print(unpacked_data)

        # stress x lattive volume
        buf = receive_all(connection, 9*8)
        unpacked_data = struct.unpack('9d', buf)
        print(unpacked_data)

        # dummy '0' at the end of the data
        buf = receive_all(connection, 4)
        unpacked_data = struct.unpack('i', buf)
        print(unpacked_data)

    connection.shutdown(socket.SHUT_RDWR)
    connection.close()


if __name__ == '__main__':
    main()
