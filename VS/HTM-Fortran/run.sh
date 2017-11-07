#!/bin/sh
export LD_LIBRARY_PATH=/tmp
export KNP_AFFINITY=compact
export OMP_NUM_THREADS=228
exec ./main_spike_fortran.out
