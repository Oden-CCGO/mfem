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
fwd=2

# Configure adjoint solve
# 0: disable adjoint operator
# 1: enable adjoint operator (one solve)
# 2: use adjoint operator to export adjoint p2o map
#    (number of solves ~ number of sensors)
adj=0

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
ref=3

# Final time
tf=0.01

# Number of time steps
nt=8

# Parameter is defined for every n-th time step
# - only used for unknown solution
# - must evenly divide the number of time steps
param_rate=2

# Observations are defined for every n-th time step (sensor frequency)
# - only used for unknown solution
# - must evenly divide the number of time steps
# - must be a multiple of param_rate
obs_rate=4

# Number of observers in x,y direction (uniformly placed)
nx_obs=3
ny_obs=3

# Enable/disable writing observations (fwd/adj output)
obs=false

# Enable/disable writing paraview vis files
vis=false

# If vis files enabled, every n-th step they are written
vs=5

# Configure program arguments
args=" -p ${p} -fwd ${fwd} -adj ${adj}"
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
if [ "$mesh_file" = true ] ; then
   args+=" -m ${m}"
fi

# Run program
./cascadia_seq ${args}
