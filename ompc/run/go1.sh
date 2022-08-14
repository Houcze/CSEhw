#!/bin/sh
#PJM -N "go1"
#PJM -L rscgrp=lecture-o
#PJM -L node=1
#PJM --omp thread=48
#PJM -L elapse=00:15:00
#PJM -g gt80
#PJM -j
#PJM -e err
#PJM -o report.lst

module load fj
export OMP_NUM_THREADS=48
export XOS_MMM_L_PAGING_POLICY=demand:demand:demand

numactl -l ./L3-rsol0
numactl -C 12-59 -m 4-7 ./L3-rsol0