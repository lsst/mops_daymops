#!/bin/sh

# --------------------------------------------------------- 
# INPUT PARAMETERS
# To run on a single host, keep nodes set equal to 1 
# Increase nodes and nslices for a larger parallel execution. 
# For example, for two nodes with 16 total cores, set nodes=2 
# and nslices=15 (pipeline itself takes one core) 
nodes=1
nslices=2
# --------------------------------------------------------- 

pwd=`pwd`
PYTHONPATH=${pwd}/../python:${PYTHONPATH}
export PYTHONPATH

# Command line arguments 
if [ "$#" != 2 ]; then
    echo "------------------------------------------"
    echo "Usage:  run.sh <policy-file-name> <runId>"
    echo "------------------------------------------"
    exit 0
fi

pipelinePolicyName=${1}
runId=${2}

# Add 1 to the number of slices to get the universe size 
usize=$(( $nslices + 1 ))

echo "nodes ${nodes}"
echo "nslices ${nslices}"
echo "usize ${usize}"

# Patch $MPDHOSTSFILE
# echo `hostname`:4 > $MPDHOSTSFILE

# MPI commands will be in PATH if mpich2 is in build
echo "Running mpdboot --totalnum=${nodes} --file=$MPDHOSTSFILE --verbose"
mpdboot --totalnum=${nodes} --file=$MPDHOSTSFILE --verbose

sleep 3s
echo "Running mpdtrace -l"
mpdtrace -l
sleep 2s

echo "Running mpiexec -usize ${usize}  -machinefile $MPDHOSTSFILE -np 1 runPipeline.py ${pipelinePolicyName} ${runId}"
mpiexec -usize ${usize}  -machinefile $MPDHOSTSFILE -np 1 runPipeline.py ${pipelinePolicyName} ${runId}

sleep 1s

# echo "Running mpdallexit"
# mpdallexit

