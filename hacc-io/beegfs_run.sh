#!/bin/bash
mpirun -np 36 \
  --hostfile mpi_hosts \
  --bind-to cpu-list:ordered --cpu-list 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 \
  --report-bindings \
  -mca pml ucx -mca btl ^vader,tcp,openib,uct \
  /home/cyx/chores/HACC_IO_KERNEL/HACC_IO 32768000 /mnt/beegfs
