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
  -x PEAKFS_HOOK_FILTER=/home/cyx/chores/HACC_IO_KERNEL/HACC_IO \
  -x LD_LIBRARY_PATH=/home/cyx/spdk/install/lib:/home/cyx/peakfs/deps/install/lib:/home/cyx/peakfs/install \
  /home/cyx/chores/HACC_IO_KERNEL/HACC_IO 32768000 /home/cyx/gekkofs_mount
