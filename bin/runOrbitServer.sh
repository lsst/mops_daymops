#!/bin/bash

ORBFIT=/home/jmyers/orb401/src/panst/orbit_server.x

for REQUEST in *.in.request
do
   BN=`basename $REQUEST .in.request`
   echo $BN | /usr/bin/time -o $BN.orbit_server.time $ORBFIT > $BN.orbit_server.log
done