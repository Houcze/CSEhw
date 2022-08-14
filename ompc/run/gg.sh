#!/bin/sh
#PJM -N "gor"
#PJM -L rscgrp=lecture-o
#PJM -L node=1
#PJM --omp thread=48
#PJM -L elapse=00:15:00
#PJM -g gt00
#PJM -j
#PJM -e err
#PJM -o testr.lst

module load fj
export OMP_NUM_THREADS=48
export XOS_MMM_L_PAGING_POLICY=demand:demand:demand

numactl  ./L3-rsol0
numactl  ./L3-rsol0
numactl  ./L3-rsol0
numactl  ./L3-rsol0
numactl  ./L3-rsol0
