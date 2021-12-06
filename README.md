# RobustBiodieselSCM

## The src folder 
This folder contains the following files with Julia code.

### *main.jl*
This file runs the program.  

There is a section for the smaller example problem and a section for the large problem.
Only deterministic optimisation can be used for the smaller example problem.  
The number of suppliers and the number of data entries to use in the SVC model (in robust data generation, file: *generate_data.jl*) are set here.

### *create_instances.jl*
The data used in the optimisation problem is created here.  

Important functions:

- **create_toy_instance()**: creates the smaller example problem data.  
- **create_general_instance()**: creates the large problem data. This includes robust data generation.

### *data_instances.jl*
This file contains the definitions for two data instance structures.  

**DetInstance** holds all the data used for deterministic optimisation.   
**RobInstance** holds the additional data used for robust optimisation.   

These data instances are created in *create_instances.jl*.

### *generate_data.jl*
This file contains functions for generating data. 

Important functions:

- **create_supply()**: creates synthetic supply data for the large problem.  
- **create_capillarity_expansion()**: implements capillarity expansion for cooking oil suppliers and updates the supply data accordingly.  
- **create_qualities()**: creates property values for biomasses of suppliers.  
- **create_robust_data()**: creates all the data for the robust problem. This includes optimising an SVC model. Returns data found in **RobInstance** in *data_instances.jl*.

These functions are used in *create_instances.jl*.

### *run_models.jl*
This file contains functions for creating or updating models and functions for optimising these models.

Model creation functions:

- **create_linear_model()**: creates the NMDT linearized model. Uses data from **DetInstance**.
- **turn_model_robust()**: adds robustness to a model. Uses data from both **DetInstance** and **RobInstance**.
- **create_nonlinear_model()**: creates a nonlinear model. Uses data from **DetInstance**. This model is used for validation of results.

Model optimisation functions:

- **run_linear_model()**: optimises the NMDT linearised model. Used for both deterministic and robust models.
- **run_nonlinear_model()**: optimises the (deterministic) nonlinear model. This is used for validation of results.

These functions are called directly from the *main.jl* file.

## Other files
The original code files (not in use):
- *original_code.jl*  
- *BioDieselOpt.jl*: contains the small example problem

## Run the program 
From the src folder run command:  
```
julia main.jl
```