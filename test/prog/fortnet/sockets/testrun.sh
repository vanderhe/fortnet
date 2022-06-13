#!/usr/bin/env bash

#  Script for starting up socket/file communication

FORTNET_CMD=$*

# Start a python server to drive the Fortnet instance
sleep 2
./prerun.py &
echo "$!" > subprocess.pid
sleep 2

# run the actual calculation
$FORTNET_CMD

# clean up afterwards
./postrun.sh
