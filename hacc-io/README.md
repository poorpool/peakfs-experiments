# HACC-IO

Source code from https://asc.llnl.gov/coral-benchmarks

Simulating N-N checkpoint workloads.

## Configuration and Compilation

The transfer size and file header size are defined in `RestartIO_GLEAN.h`. Make sure the transfer size is what you want, and make sure the file header size is a multiple of the transfer size.

```bash
./build.sh
```

We slightly modified the code so that HACC-IO uses N-N mode and doesn't preallocate files.

## Run the program

The first argument to HACC-IO is the number of particles. A checkpoint file consists of several read/write parts: 1 file header + 7 float arrays (4 bytes) + 1 int64_t array (8 bytes) + 1 uint16_t array (2 bytes), so ensure that **the particle count times two is a multiple of the transfer size**.

```bash
./beegfs_run.sh
./maplefs_run.sh
./mcatfs_run.sh
```

There are 32768000 particles per process, 36 test processes on the client node.