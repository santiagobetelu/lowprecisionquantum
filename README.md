## A low precision quantum simulation library

- For simulations of quantum circuits with Schrodinger's formulation (unoptimized)
- Uses low precision arithmetic to save memory and bandwidth
- The paper "The limits of quantum circuit simulation with low precision arithmetic" contains the mathematical analysis of the error. Between 16 and 32 bits per coefficient are sufficient for almost any simulation
- Most precise when the states are random and maximally entangled.

## Paper:
https://arxiv.org/abs/2005.13392

## How to compile and run
A test program is provided. It runs the random circuit test of Table VII. Compile with mpicc as
```bash
mpicc -march=native -Ofast test-smallcomplex.c -o test -lm
```
In a single node computer, run as:
```bash
mpirun ./test
```
On a multinode system run with Slurm 
```bash
sbatch slurmscript.batch
```
Suggested Slurm script:
```bash
#SBATCH –o output-%J.txt
#SBATCH --nodes=8     # 8 nodes (must be a power of 2)
#SBATCH –n 64         # 8 ranks per node (must be a power of 2)
#SBATCH –p normal     # what queue
#SBATCH –t 02:00:00   
srun ./test
```
The number of ranks must be a power of two.

If you want to test your own quantum circuits, do not edit the library, but write the quantum program in a separate file containing the function "qcprogram()" as shown in the example.
