#!/bin/bash

if [ -z ${INFILE+x} ]; then
    echo "ERROR: INFILE environment variable is not specified" >&2
    exit 1
fi

if [ ! -e "/var/data/$INFILE" ]; then
    echo "ERROR: Specified INFILE does not exist" >&2
    exit 2
fi

outfile="${OUTFILE:-outfile}"
igrf="${IGRF:-13}"
mode="${MODE:-par}"
steps="${STEPS:-25000}"

./build/Trajectories_IGRF_T04_C "/var/data/$INFILE" "/var/data/$outfile" "$igrf" "$mode" "$steps"