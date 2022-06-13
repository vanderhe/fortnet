#!/usr/bin/env bash

#  Script for cleaning up after socket/file communication

if [ -e file.txt ]; then
    # kill anything using the file
    fuser -k $(cat file.txt)

    # Remove any left over socket file in /tmp
    rm -f $(cat file.txt)

    # remove the file itself
    rm file.txt
fi

if [ -e port.txt ]; then
    # kill anything using that port number
    fuser -n tcp -k $(cat port.txt)

    # clean up file
    rm -f port.txt
fi
