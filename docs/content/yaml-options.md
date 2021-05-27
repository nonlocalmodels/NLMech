# Configuration file

We use the [YAML](https://de.wikipedia.org/wiki/YAML) file format for the configuration files. For more details on the YAML file format, we refer to this [tutorial](https://gettaurus.org/docs/YAMLTutorial/). Here, we introduce the option of within the YAMl file to control the simulation.  

## Mesh generation

To generate a mesh the tool [mesh]() is provided and it is used as

```sh
mesh -i input_mesh.yaml -d 1
```

where `-d` specify the dimension of the problem and `-i` specifies the YAML fine with the mesh information. Let us look at the [example](https://github.com/nonlocalmodels/NLMech/tree/main/examples/qsModel/1D) for the one-dimensional implicit time integration example:


```yaml
Output:
  Path: ./
  Mesh: mesh
  File_Format: vtu
Domain:
  - 0.
  - 16
Horizon: 2
Horizon_h_Ratio: 4 
Compress_Type: zlib
```

### Output

The tag `Output` describes details about the generate `.vtu` files using the following attributes:

* `Path:` Describes the path where the generated mesh is written. Using the value `./` means the file will be written in the folder where the executable `mesh` was executed. 
* `Mesh:` Describes the name of the generated file. 
* `File_Format`: Describes the file format of the generated file, e.g. `vtu` and `gmsh`. For the `vtu` file format we refer to the [VTK](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) file format documentation and for the `gmsh` file format to the [Gmsh](http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php) documentation.
* `Compress_Type` Describes the compression type of the generated file. 

### Domain

The tag `Domain` describes the simulation domain using the following attributes:

* The dashes `-` list the first point and the last point of the domain
* `Horizon` Describes the length of the horizon as an integer
* `Horizon_h_Ratio` Describes the number of nodes within the horizon

For more complex geometries, the [Gmsh](https://gmsh.info/) tool is recommended. Further examples are available [here](https://nonlocalmodels.github.io/examples/).

#### References

* C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering 79(11), pp. 1309-1331, 2009. 
* W. Schroeder, Ken Martin, and Bill Lorensen, Visualization Toolkit: An Object-Oriented Approach to 3D Graphics, 4th Edition

