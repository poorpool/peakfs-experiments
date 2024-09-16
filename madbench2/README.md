# MADbench2

Source code from https://crd.lbl.gov/divisions/scidata/c3/c3-research/madbench2/

Simulating N-1 matrix saving/loading workloads.

## Configuration and Compilation

During compilation, we set the mode to IO (replacing actual computation with busy work) and specify the file path in `MADbench2.h``.

```bash
./build.sh
```

We made some simple modifications to the test code: replaced `fread` with `read`, and split the `fwrite` into multiple `write` operations with a maximum size of the transfer size, followed by one `fsync`.

## Run the program

We set NO_PIX to 24576, NO_BIN to 16, NO_GANG to 1, SBLOCKSIZE and FBLOCKSIZE to 8, and RMOD and WMOD to 6. Their meanings can be found at https://crd.lbl.gov/divisions/scidata/c3/c3-research/madbench2/

MADbench2 has three phases: S, W, and C. After each phase, we manually execute `echo 3 > /proc/sys/vm/drop_caches` to clear the cache. We use the official Excel spreadsheet to calculate the write bandwidth for the W phase and the read bandwidth for the C phase. See `madbench2.xls` for details.

There are 36 test processes on the client node.

```bash
./beegfs_run.sh
./gekkofs_run.sh
./peakfs_run.sh
```