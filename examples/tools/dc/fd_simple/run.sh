#!/bin/bash
MY_PWD=$(pwd)

# exe
cd ../
dc="$(pwd)/bin/dc"

# generate input and mesh file
cd $MY_PWD

"$dc" -i input.yaml -k fd_simple
