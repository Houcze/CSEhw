#!/bin/sh
#PJM -N "gor"
#PJM -L rscgrp=lecture-o
#PJM -L node=1
#PJM --omp thread=12
#PJM -L elapse=00:15:00
#PJM -g gt00
#PJM -j
#PJM -e err
#PJM -o testr.lst

module load fj
export OMP_NUM_THREADS=12
export XOS_MMM_L_PAGING_POLICY=demand:demand:demand

numactl -l ./L3-sol0
numactl -l ./L3-sol0
numactl -l ./L3-sol0
numactl -l ./L3-sol0
numactl -l ./L3-sol0
numactl -C 12-23 -m 4 ./L3-rsol0
numactl -C 12-23 -m 4 ./L3-rsol0
numactl -C 12-23 -m 4 ./L3-rsol0
numactl -C 12-23 -m 4 ./L3-rsol0
numactl -C 12-23 -m 4 ./L3-rsol0
