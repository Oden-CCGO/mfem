#!/bin/bash
#----------------------------------------------------
# Example Slurm job script
# for TACC Stampede2 SKX nodes
#
#   *** Hybrid Job on SKX Normal Queue ***
#
#       This sample script specifies:
#         10 nodes (capital N)
#         40 total MPI tasks (lower case n); this is 4 tasks/node
#         12 OpenMP threads per MPI task (48 threads per node)
#
# Last revised: June 2021
#
# Notes:
#
#   -- Launch this script by executing
#      "sbatch fiber.slurm" on Stampede2 login node.
#
#   -- Check current queue with: "squeue -u sh43394"
#
#   -- Use ibrun to launch MPI codes on TACC systems.
#      Do not use mpirun or mpiexec.
#
#   -- In most cases it's best to keep
#      ( MPI ranks per node ) x ( threads per rank )
#      to a number no more than 48 (total cores).
#
#   -- If you're running out of memory, try running
#      fewer tasks and/or threads per node to give each
#      process access to more memory.
#
#   -- IMPI and MVAPICH2 both do sensible process pinning by default.
#
#----------------------------------------------------

#SBATCH -J casc            # Job name
#SBATCH -o casc.o%j        # Name of stdout output file
#SBATCH -e casc.e%j        # Name of stderr error file
#SBATCH -p skx-normal      # Queue (partition) name
#SBATCH -N 1               # Total # of nodes
#SBATCH -n 1               # Total # of mpi tasks
#SBATCH -t 01:59:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=stefan@oden.utexas.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A FTA-SUB-Ghattas # Allocation name (req'd if you have more than 1)

# Other commands must follow all #SBATCH directives...

module list
pwd
date

#!/bin/bash

#
# Cascadia Codes
# Run script: cascadia_seq
#

# Mesh file (if any)
mesh_file=false
m='mesh_file'

# Problem instance (forward problem)
#   1   - manufactured stationary solution (sine-wave)
#   2   - manufactured stationary solution (linear polynomial)
#   3   - manufactured stationary solution (quadratic polynomial)
#   10  - manufactured time-dependent solution (exp decaying sine-wave)
#   20  - manufactured time-dependent solution (exp decaying linear polynomial)
#   30  - manufactured time-dependent solution (exp decaying quadratic polynomial)
#   40  - manufactured time-dependent solution (linear growing quadratic polynomial)
#   50  - manufactured time-dependent solution (quadratically growing quadratic polynomial)
#   100 - unknown solution: forcing with single (x,y) Gaussian deformation)
#   200 - unknown solution: forcing with superimposed Gaussian deformations)
p=100

# Configure forward solve
# 0: disable forward operator
# 1: enable forward operator (one solve)
# 2: use forward operator to export p2o map
#    (number of solves ~ number of parameters)
fwd=1

# Configure adjoint solve
# 0: disable adjoint operator
# 1: enable adjoint operator (one solve)
# 2: use adjoint operator to export adjoint p2o map
#    (number of solves ~ number of sensors)
adj=0

# Configure the prior (regularization)
# 0: Do not assemble prior
# 1: Laplacian prior (assemble + write to file)
# 2: Bi-Laplacian prior (assemble + write to file)
prior=0

# Regularization parameters
# alpha1 ~ |m|
# alpha2 ~ |grad m|
# alpha3 ~ |dm/dt|
alpha1=1.0
alpha2=1.0
alpha3=0.0

# Polynomial order of approximation
#   order_p = o   (scalar-valued H1 pressure space)
#   order_u = o-1 (vector-valued L2 velocity space)
#   order_m = 1   (scalar-valued H1 parameter space)
o=2

# ODE solver type
#   1: Forward Euler
#   2: RK2
#   3: RK3 SSP
#   4: RK4
#   6: RK6
#  11: Backward Euler
ode=4

# Enable/disable mass lumping
lump=true

# Number of uniform h-refinements
ref=0

# Final time
tf=0.2

# Number of time steps
nt=20

# Parameter is defined for every n-th time step
# - only used for unknown solution
# - must evenly divide the number of time steps
param_rate=1

# Observations are defined for every n-th time step (sensor frequency)
# - only used for unknown solution
# - must evenly divide the number of time steps
# - must be a multiple of param_rate
obs_rate=1

# Number of observers in x,y direction (uniformly placed)
nx_obs=2
ny_obs=3

# Specify format for output data
# hdf = true  --> binary (HDF5)
# hdf = false --> text
hdf=true

# Enable/disable writing observations (fwd/adj output)
obs=true

# Enable/disable writing paraview vis files
vis=true

# If vis files enabled, every n-th step they are written
vs=5

# Configure program arguments
args=" -p ${p} -fwd ${fwd} -adj ${adj}"
args+=" -prior ${prior} -alpha1 ${alpha1} -alpha2 ${alpha2} -alpha3 ${alpha3}"
args+=" -o ${o} -ode ${ode} -ref ${ref}"
args+=" -tf ${tf} -nt ${nt}"
args+=" -pr ${param_rate} -or ${obs_rate}"
args+=" -nxo ${nx_obs} -nyo ${ny_obs}"
if [ "$lump" = true ] ; then
   args+=" -lump"
else
   args+=" -no-lump"
fi
if [ "$obs" = true ] ; then
   args+=" -obs"
else
   args+=" -no-obs"
fi
if [ "$vis" = true ] ; then
   args+=" -vis -vs ${vs}"
else
   args+=" -no-vis"
fi
if [ "$hdf" = true ] ; then
   args+=" -hdf"
else
   args+=" -no-hdf"
fi
if [ "$mesh_file" = true ] ; then
   args+=" -m ${m}"
fi

# Run program
ibrun ./cascadia_seq ${args}


#
# set output
mkdir ${SLURM_JOB_ID}
cd ${SLURM_JOB_ID}
mkdir outputs
cd ..
dir_output="${SLURM_JOB_ID}/outputs/"

date
# ---------------------------------------------------
