# PeakFS Experiments

## Metadata and data performance: mdtest and IOR

[IO500](https://github.com/IO500/io500) is a comprehensive benchmark suite for testing parallel file systems. It uses mdtest to measure metadata throughput (KIOPS) and IOR (GiB/s) to measure data bandwidth. The IO500 test suite includes many sub-tests, with fixed parameters controlling the I/O patterns. Users only need to specify parameters representing the test data volume. The PeakFS paper involves the following sub-tests:

- mdtest-easy-write/stat/delete: Metadata write/stat/delete throughputs on empty files.
  - Each test process creates its own directory and create empty files within it.
  - Test processes are shuffled to perform stat and delete operations on files created by other processes.
  - Corresponding to Fig.10, Fig.14, and Fig.16(a) in the PeakFS paper.
- mdtest-hard-write/stat/delete: Metadata write/stat/delete throughputs on small files.
  - Every test process creates and writes files in the same shared directory.
  - The file size is 3901 bytes.
  - Test processes are shuffled to perform stat and delete operations on files created by other processes.
  - Corresponding to Fig.11 in the PeakFS paper.
  - Read throughput was also tested, but due to space limitations and its similar trend to other sub-tests, the results were not included in the paper.
- ior-easy-write/read: Data write/read bandwidths under N-N I/O pattern
  - Each test process creates its own file and writes it (known as N-N I/O).
  - Test processes are shuffled to read the data written by other processes.
  - Corresponding to Fig.12, Fig.15, and Fig.16(b) in the PeakFS paper.
- ior-hard-write/read: Data write/read bandwidths under N-1 I/O pattern
  - Every test process writes the same shared file without overlap (known as N-1 I/O).
  - The I/O transfer size is 47008 bytes.
  - Test processes are shuffled to read the data written by other processes.
  - Corresponding to Fig.13 in the PeakFS paper.

Due to confidentiality requirements of our paper collaborators, we regret that we cannot publish the code of PeakFS and our own implementation of GekkoFS. However, we are pleased to show the configurations and scripts we used in `mdtest_and_ior_config/` during testing for reference. We adjust the data volume in the configuration to ensure that the tests run for a minimum duration (e.g., 40 seconds). 

To test the metadata and data performance of a parallel file system like BeeGFS, simply run `./beegfs_io500.sh config-beegfs.ini`.

## Comprehensive Experiments

### Mixed metadata-intensive workloads: md-workbench

md-workbench generates a semi-synchronous metadata-intensive workload that was designed to mimic what a parallel compilation may look like to a file system. The extended mode of the IO500 test suite includes this test. md-workbench consists of three stages: precreate, benchmark, and cleanup. We use the metadata throughput from the second stage, benchmark.

The configurations and scripts are in `md-workbench_config/`, corresponding to Fig.17 in the PeakFS paper. We adjust the data volume of mdtest-easy to ensure that the second stage of md-workbench runs for at least 30 seconds.

To test the mixed metadata performance of a parallel file system like BeeGFS, simply run `./beegfs_mdwb.sh mdwb-beegfs.ini`.

### Checkpoint/restart workloads: HACC-IO

[HACC-IO](ttps://asc.llnl.gov/coral-benchmarks) is the I/O kernel of the cosmology N-body simulation application HACC. It simulates the process of saving and loading particle information between HACC's computing stages, i.e., check-point/restart workloads. The source code, configurations, and scripts can be found in `hacc-io`.