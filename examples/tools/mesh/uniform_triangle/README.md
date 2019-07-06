### Mesh
We describe the input file required to run [Mesh](../../../../../tools/mesh/mesh.cpp). 

- We consider rectangle domain of size `0.1 m x 0.1 m`. 

- Horizon is specified to be `horizon = 0.02 m`. 

- Ratio of horizon to mesh size `r = 4`. Thus mesh size is `horizon / 4`.

- Type of mesh is `uniform_tri` which produces uniform mesh of triangles. Other options are `uniform_square` and `uniform_tri_sym`.

### Results
Mesh looks like

<p id="result" align="center">
	<img src="result.png" alt="setup" width="400" height="400" />
</p>
