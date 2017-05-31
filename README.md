# Nelder Mead simplex
## Optimization school project

To compile the code, in linux command line, type:
```
gcc -g nelder_mead.c functions.c -lm -o nelder_mead
```
To run the algorithm, in linux command line, type:
```
nelder_mead dim func_num min max num_res

where, 
  - dim: dimension of the objective function
  - func_num: objective function id 
  - min, max: the search interval
  - num_res: number of restarts
```
#### For results and more check the report
