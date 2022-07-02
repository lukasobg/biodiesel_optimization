# RobustBiodieselSCM

## Run the program 
Run the following command from the src folder:

For the main program
```
julia main.jl
```
For the toy problem
```
julia toy_problem.jl
```
For the benchmarks
```
julia run_benchmarks.jl
```

----

## The src folder 
This folder contains the following files with Julia code.


**Files that run the program.**
&nbsp;
### *main.jl*
This file runs the main program.  
&nbsp;
### *toy_problem.jl*
This file runs the toy problem.  
&nbsp;
### *run_benchmarks.jl*
This file runs benchmarks.  

----

**Other files**

&nbsp;
### *create_instances.jl*
The data used in the optimisation problem is created here.  

Important functions:

- **create_toy_instance()**: creates the smaller example problem data.  
- **create_deterministic_instance()**: creates the large problem data.

&nbsp;
### *data_instances.jl*
This file contains the definitions for three data instance structures.  

**ToyInstance** holds all the data used in the toy problem.   
**DetInstance** holds all the data used in deterministic optimisation.   
**RobInstance** holds the additional data used in robust optimisation.   
**SVCInstance** holds the additional data used in robust optimisation.   

These data instances are created in *create_instances.jl*.

&nbsp;
### *generate_data.jl*
This file contains functions for generating data. 

Important functions:

- **create_supply()**: creates synthetic supply data for the large problem.  
- **create_capillarity_expansion()**: implements capillarity expansion for cooking oil suppliers and updates the supply data accordingly.  
- **create_qualities()**: creates property values for biomasses of suppliers.  
- **create_robust_data()**: creates all the data for the robust problem. This includes optimising an SVC model. Returns data found in **RobInstance** in *data_instances.jl*.

These functions are used in *create_instances.jl*.

&nbsp;
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

&nbsp;
### *benchmarks.jl*
This file contains functions for benchmarking. Some changes required when benchmarking a deterministic v robust model. 
- **benchmark_model()**: runs 11 iterations of benchmarking. 11 data instances, 4 models each: nlm, nmdt with p=2, nmdt with p=3, nmdt with p=4.
- **save_benchmark()**: prints out benchmark times and stores in csv files.

&nbsp;
### *model_results.jl*
- **get_results()**: gets a solved model as parameter, prints out and returns the relevant data.

----

## The data folder 

This folder contains the model results and benchmark data as csv-files. Separate folders for deterministic and robust model benchmarks. 

----

## Other files
The original code files (not in use):
- *original_code.jl*  
- *BioDieselOpt.jl*: contains the small example problem