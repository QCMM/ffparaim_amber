#!/bin/bash

export dir=$1

for i in $dir/*/* ; do
    cd $i
    # orca run
    sbatch ../../../submit_orca.sh
    # return to dir
    cd ../../..
done
