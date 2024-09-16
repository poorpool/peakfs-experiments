#!/bin/bash
# 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 是 Numa 0
# 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 是 Numa 1
echo "remember to modify .h mount point"
mpirun -np 36 \
  --hostfile mpi_hosts \
  --bind-to cpu-list:ordered --cpu-list 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 \
  --report-bindings \
  -mca pml ucx -mca btl ^vader,tcp,openib,uct \
  -x FILETYPE=SHARED \
  -x IOMETHOD=POSIX \
  -x BWEXP=1.0 \
  /home/cyx/chores/madbench2/MADbench2.x 24576 16 1 8 8 6 6 
