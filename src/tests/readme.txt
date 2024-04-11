##### Readme #####
To run the code, I suggest to open the whole folder of the project with VSCode, 
set the number of threads of the julia configuration file and run the dampf.jl file.
Otherwise just run the dampf.jl with the terminal command

    julia -t 4 dampf.jl

Where "-t nthreads" sets the number of threads. The plots will be saved in the plots folder. 
If you run the code inside this folder the plots will be saved under the location src/tests/plots.