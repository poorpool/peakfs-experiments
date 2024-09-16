#!/bin/bash
mpirun -np 36 \
  --hostfile mpi_hosts \
  --bind-to cpu-list:ordered --cpu-list 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 \
  --report-bindings \
  -mca pml ucx -mca btl ^vader,tcp,openib,uct \
  -x RDMAV_HUGEPAGES_SAFE=1 \
  -x LD_PRELOAD=/home/cyx/peakfs/install/libpeakfs_preload.so \
  -x PEAKFS_CFG_PATH=/home/cyx/peakfs/cfg/cfg-pdsl.json \
  -x PEAKFS_VIEW_PATH=/home/cyx/peakfs/cfg/view.json \
  -x PEAKFS_HOOK_FILTER=/home/cyx/chores/madbench2/MADbench2.x \
  -x LD_LIBRARY_PATH=/home/cyx/spdk/install/lib:/home/cyx/peakfs/deps/install/lib:/home/cyx/peakfs/install \
  -x FILETYPE=SHARED \
  -x IOMETHOD=POSIX \
  -x BWEXP=1.0 \
  /home/cyx/chores/madbench2/MADbench2.x 24576 16 1 8 8 6 6 
