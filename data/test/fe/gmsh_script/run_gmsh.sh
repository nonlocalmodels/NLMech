#!/bin/bash 
MY_PWD=$(pwd)

filename="mesh"
fileout="mesh"

# create .geo file
octave "$filename"".m" > "$fileout"".geo"

# use gmsh to create mesh
gmsh "$fileout"".geo" -2 

# use gmsh to create mesh
gmsh "$fileout"".geo" -2 -o "$fileout"".vtk"

# clean 
rm "$fileout"".geo"

# # use paraview python script to convert vtk file to vtu
# pvpython convert_vtk.py

# # clean 
# rm mesh.vtk

# # rename file to desired name
# if [[ $# == 1 ]]; then
# 	filename=$1
# 	mv mesh.vtu "$filename"".vtu"
# fi