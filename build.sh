#!/bin/sh

root_dir=$(dirname "$(realpath "$0")")

cmake -S "$root_dir"/src -B "$root_dir"/build && cd "$root_dir"/build && make