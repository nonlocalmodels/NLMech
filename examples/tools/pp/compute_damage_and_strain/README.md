### Postprocessing simulation results
We provide various methods in [Compute](../../../../../tools/pp/src/compute.cpp) to perform postprocessing on simulation results. In this example, we will compute damage at nodes and strain and stress.

We assume that input `.yaml` file used in obtaining simulatio is available and is in the same path where simulation results `.vtu` files are.

### Input file for postprocessing

1. Provide path where simulation results and input file is located.

1. Specifiy the prefix the simulation filenames have. E.g. if files are `output_1.vtu, output_2.vtu,...` then prefix is `output`.

1. If simulation input file does not have elastic material properties, then provide them.

1. Specify output path.

1. We can specify multiply compute sets. In this example we have two compute sets. One for damage, and another for strain and stress. Thus `Sets: 2`.

1. For each set give the required instruction. 

#### Damage 
Damage calculation is performed in `Set_1`.  

- Provide filename for this set. We have `Tag_Filename: damage_Z`.

- Specify (optional) start and end index of simulation file, and specify (optional) interval between processing each file. 

	- For example, we want to skip first 10 files, we will have `Dt_Start: 11`. 

	- If we want to skip 1 file each time, i.e. consider `output_11.vtu, output_13.vtu, ...` then we will have `Dt_Interval: 2`.

- Set `Damage_Z: true` which is responsible for computing damage. 

#### Strain and stress
This calculation is performed in `Set_2`.  

- Description of some fields are same as before. However, we can have different values for them in different sets. 
	
	- We take `Tag_Filename: strain`, `Dt_Start: 20`, `Dt_Interval: 4`.

- `Output_Only_Nodes: true` means we do not write element-node connectivity in the output `.vtu` file. As strain/stress are cell data, we set it to `false`.

- Set `Strain_Stress: true` which is responsible for computing strain and stress. 

- Additionally, we can compute the magnitude of strain tensor. 

	- If we want to compute magnitude of tensor, we only need to have 'Magnitude_Strain_Tensor:' with no entry. 

	- If we want to compute magnitude of `yy` component, we can specify `Component: yy`. 
