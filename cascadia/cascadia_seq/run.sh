#!/bin/bash

#
# Cascadia Codes
# Run script: cascadia_seq
#

# Mesh file (if any)
mesh_file=false
m='mesh_file'

# Output directory
d='output'

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
# 0: write adjoint vectors in standard ordering
# 1: write adjoint vectors in block-reverse ordering
adj_reverse_order=0

# Configure the prior (regularization)
#  0: Do not assemble prior
#  1: Laplacian prior (assemble)
# 11: Laplacian prior (assemble + write to file)
#  2: Bi-Laplacian prior (assemble)
# 22: Bi-Laplacian prior (assemble + write to file)
prior=0
#  0: Write prior in standard indexing (time(outer), space(inner))
#  1: Write prior in re-indexed dofs (space(outer), time(inner))
prior_reindex=0

# Regularization parameters
# alpha1 ~ |m|
# alpha2 ~ |grad m|
# alpha3 ~ |dm/dt|
alpha1=1.0e-3
alpha2=1.0e-3
alpha3=1.0e-5

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
tf=10.0

# Number of time steps
nt=1000

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
# choose: 2,3,5,9,17 [0.25  0.75 ]
# choose: 2,3,7,13   [0.125 0.875]
nx_obs=3
ny_obs=5

# Specify format for output data
# hdf = true  --> binary (HDF5)
# hdf = false --> text
hdf=true

# Enable/disable writing observations (fwd/adj output)
obs=true

# Enable/disable writing paraview vis files
vis=true

# If vis files enabled, every n-th step they are written
vs=10

# Configure program arguments
args=" -d ${d} -p ${p} -fwd ${fwd} -adj ${adj} -revadj ${adj_reverse_order}"
args+=" -prior ${prior} -indprior ${prior_reindex}"
args+=" -alpha1 ${alpha1} -alpha2 ${alpha2} -alpha3 ${alpha3}"
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
./cascadia_seq ${args}
