# RobustBiodieselSCM
Code for the Robust procument fro biodiesel production

## Folders and files

The **src** folder contains the following files with julia code.

- *main.jl*: file runs the program
- *create_instances.jl*: the data used in the optimization problem is created here 
- *data_instances.jl*: the data instance structs are defined here 
- *generate_data.jl*: contains functions for creating supply, capillarity, quality, and the robust data (functions used in *create_instances.jl*)
- *run_models.jl*: runs the different optimization models: NMDT linearized, non-linear and robust

The original code files:
- *original_code.jl*  
- *BioDieselOpt.jl*: contains the small example problem