#!/bin/bash
mpirun -np 36 \
  --hostfile mpi_hosts \
  --bind-to cpu-list:ordered --cpu-list 22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 \
  --report-bindings \
  -mca pml ucx -mca btl ^vader,tcp,openib,uct \
  -x LD_PRELOAD=/home/cyx/gekkofs/build/src/libgekkofs_interceptor.so \
  -x GEKKOFS_CONFIG_PATH=/home/cyx/gekkofs/config/interceptor_test.yaml \
  -x GEKKOFS_INTERCEPT_PROCESS=/home/cyx/chores/madbench2/MADbench2.x \
  -x GEKKOFS_LOG_NAME=/tmp/gekkofs-interceptor.log \
  -x LD_LIBRARY_PATH=/home/cyx/gekkofs_deps/install/lib \
  -x FILETYPE=SHARED \
  -x IOMETHOD=POSIX \
  -x BWEXP=1.0 \
  /home/cyx/chores/madbench2/MADbench2.x 24576 16 1 8 8 6 6  
