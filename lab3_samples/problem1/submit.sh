#!/bin/bash

# SBATCH --ntasks=384
#SBATCH --nodes=8
#SBATCH --partition=xeon-p8 #(xeon-p8: 48 cores/nodes, xeon-g6-volta: 40 cores/nodes)

python 1C_run_input.py