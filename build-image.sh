#!/bin/sh

root_dir=$(dirname "$(realpath "$0")")

docker build -t trajectories_igrf_t04_c:latest "$root_dir"