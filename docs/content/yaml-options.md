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

## Simulation

Here, all tags and attributes of ther YAML files are listed and described. A detailed set of examples is available [here](https://nonlocalmodels.github.io/examples/).

### Model deck

Example of a model deck:

```yaml
Model:
  Dimension: 2
  Discretization_Type:
    Spatial: finite_difference
    Time: velocity_verlet
  Final_Time: 0.001000
  Time_Steps: 50000
  Horizon: 0.002000
  Horizon_h_Ratio: 4
```

The tag `Model` describes the model with following attributes:

* `Dimension` Defines the dimension of the problem, e.g. one-dimensional 
* `Discretization_Type` Describes the discretization in time and space
  * `Time` Defines if a central difference scheme `central_difference` or a velocity verlet `velocity_verlet` scheme is used for the discretization in time
  * `Spatial` Defines if a finite difference `finite_difference` scheme is used
* `Horizon` Defines the horizon
* `Horizon_h_Ratio` Defines the ratio of nodes with the distance of h in the horizon
* `Mesh_Size` Defines the mesh size h
* `Final_Time` Defines the final time 
* `Time_Steps` Defines the amount of time steps

### Restart

Example of a `Restart` deck:

```yaml
Restart:
  File: out/output_3.vtu
  Step: 1500
```

The tag `Restart` describes a restart of the simulation using following attributes:

* `File` The path and file name of the last sucessfull simulation step
* `Step` The time step to restart from

One example of the restart of one simulation is available [here](https://nonlocalmodels.github.io/examples/restart.html).

### Mesh

Example of a `Mesh` deck:

```yaml
Mesh:
  File: mesh.vtu
  Keep_Element_Conn: true
```

The tag `Mesh` describes the mesh of the simulation using following attributes:

* `File` Path and file name to the mesh either in the gmsh or vtu file format.
* `Keep_Element_Conn` Keep the mesh information available

### Output

Example of a `Output` deck

```yaml
Output:
  Path: out_top_bottom_const/
  Tags:
    - Displacement
    - Velocity
    - Force
    - Damage_Z
  Output_Interval: 500
  Compress_Type: zlib
  Perform_FE_Out: true
```

The tag `Output` describes the mesh of the simulation using following attributes:

* `File_Format` Specifies if the gmsh or vtu file format is used
* `Path` Defines the path were the output is written to
* `Compress_Type` Defines the compression type for the vtu file (`zlib` or `ascii`)
* `Output_Interval` Defines the interval a output file is written
* `Tag` Defines the field appended to the out file. The mesh is always added and all other simulation data needs to be listed.
* `Perform_FE_Out` Store the mesh information in the output

### Boundary conditions

#### Displacement boundary conditions

Example of a `Displacement_BC` deck

```yaml
Displacement_BC:
  Sets: 1
  Set_1:
    Location:
      Rectangle: [0.000000e+00, 9.900000e-02, 1.000000e-01, 1.000000e-01]
    Direction: [1,2]
    Time_Function:
      Type: constant
      Parameters:
        - 0.0
    Spatial_Function:
      Type: constant
      Parameter: [1]
```


The tag `Displacement_BC` describes the displacement boundary conditions and the attribute `Sets` speifies the amount of boundary conditions. If boundary conditions is defined by `Set_1`, `Set_2`, and so on using following attributes:

* `Location` Defines the location of the boundary condition
  *  `Line: [ a, b ]` Applies the boudnary condition along the line from `a` to `b`
  *  `Rectangle: [a, b, c, d]` Applies the boundary condition on all nodes inside the recangle given by the two points `(a,b)` and `c,d)`  
* `Direction: [n]` Defines the direction (x=1, y=2, or z=3) of the boundary condition 
* `Time_Function` Specifies how the boundary condition is applied over time
  * `Type` Defines the type of loading
  * `Parameters` Defines the corresponding parameter of the load
* `Spatial_Function` Specifies how the boundary condition is applied in space
  * `Type` Defines the type of loading
  * `Parameters` Defines the corresponding parameter of the load

#### Force boundary conditions

Example of a `Force_BC` deck:

```yaml
Force_BC:
  Sets: 2
  Set_1:
    Location:
      Rectangle: [0, 0, 0.1, 0.002]
    Direction: [2]
    Time_Function:
      Type: linear
      Parameters:
        - -5000000000.0
    Spatial_Function:
      Type: constant
      Parameters: [1]
  Set_2:
    Location:
      Rectangle: [0, 0.098, 0.1, 0.1]
    Direction: [2]
    Time_Function:
      Type: linear
      Parameters:
        - 5000000000.0
    Spatial_Function:
      Type: constant
      Parameters: [1]
```

The tag `Force_BC` describes the displacement boundary conditions and the attribute `Sets` speifies the amount of boundary conditions. If boundary conditions is defined by `Set_1`, `Set_2`, and so on using following attributes:

* `Location` Defines the location of the boundary condition
  *  `Line: [ a, b ]` Applies the boudnary condition along the line from `a` to `b`
  *  `Rectangle: [a, b, c, d]` Applies the boundary condition on all nodes inside the recangle given by the two points `(a,b)` and `c,d)`  
* `Direction: [n]` Defines the direction (x=1, y=2, or z=3) of the boundary condition 
* `Time_Function` Specifies how the boundary condition is applied over time
  * `Type` Defines the type of loading
  * `Parameters` Defines the corresponding parameter of the load
* `Spatial_Function` Specifies how the boundary condition is applied in space
  * `Type` Defines the type of loading
  * `Parameters` Defines the corresponding parameter of the load

### Fracture


Example of a `Fracture` deck:

```yaml
Fracture:
  Cracks:
    Sets: 1
    Set_1:
      Orientation: -1
      Line: [5.000500e-02, 0.000000e+00, 5.000500e-02, 2.000000e-02]
```

The tag `Fracture` describes the displacement boundary conditions and the attribute `Sets` speifies the amount of boundary conditions. If boundary conditions is defined by `Set_1`, `Set_2`, and so on using following attributes:

* `Orientation` Describes the orientation of the crack
* `Line` Describes a line and all bonds interesting this line are initially broken

#### References

* C. Geuzaine and J.-F. Remacle. Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities. International Journal for Numerical Methods in Engineering 79(11), pp. 1309-1331, 2009. 
* W. Schroeder, Ken Martin, and Bill Lorensen, Visualization Toolkit: An Object-Oriented Approach to 3D Graphics, 4th Edition

