# DAMPF - Dissipation-Assisted Matrix Product Factorization
Development and optimizations of the Dissipation-Assisted Matrix Product Factorization (DAMPF) algorithm for the simulation of open quantum systems

### Suggestions to run the code
To run the code, I suggest to open the whole folder of the project with VSCode, 
set the number of threads of the julia configuration file and run the dampf.jl file.
Otherwise just run the dampf.jl with the terminal command

    julia -t 4 dampf.jl

where "-t nthreads" sets the number of threads. The plots are saved in the plots folder.