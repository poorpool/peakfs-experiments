#!/bin/bash
# 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 是 Numa 0
# 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 是 Numa 1
mpirun -np 36 \
  --hostfile mpi_hosts \
  --bind-to cpu-list:ordered --cpu-list 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 \
  --report-bindings \
  -mca pml ucx -mca btl ^vader,tcp,openib,uct \
  -x LD_PRELOAD=/home/cyx/gekkofs/build/src/libgekkofs_interceptor.so \
  -x GEKKOFS_CONFIG_PATH=/home/cyx/gekkofs/config/interceptor_test.yaml \
  -x GEKKOFS_INTERCEPT_PROCESS=/home/cyx/chores/HACC_IO_KERNEL/HACC_IO \
  -x GEKKOFS_LOG_NAME=/tmp/gekkofs-interceptor.log \
  -x LD_LIBRARY_PATH=/home/cyx/gekkofs_deps/install/lib \
  /home/cyx/chores/HACC_IO_KERNEL/HACC_IO 32768000 /home/cyx/gekkofs_mount
