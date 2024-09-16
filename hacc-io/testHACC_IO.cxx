#include "mpi.h"
#include <errno.h>
#include <fcntl.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "RestartIO_GLEAN.h"

using namespace std;

int main(int argc, char *argv[]) {
  char *fname = 0;
  char *buf = 0;
  int numtasks, myrank, status;
  MPI_File fh;

  status = MPI_Init(&argc, &argv);
  if (MPI_SUCCESS != status) {
    printf(" Error Starting the MPI Program \n");
    MPI_Abort(MPI_COMM_WORLD, status);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (argc != 3) {
    printf(" USAGE <exec> <particles/rank>  < Full file path>  ");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  int64_t num_particles = atoi(argv[1]);

  fname = (char *)malloc(strlen(argv[2]) + 1);
  strncpy(fname, argv[2], strlen(argv[2]));
  fname[strlen(argv[2])] = '\0';

  // Let's Populate Some Dummy Data
  float *xx, *yy, *zz, *vx, *vy, *vz, *phi;
  int64_t *pid;
  uint16_t *mask;

  xx = (float *)malloc(num_particles * sizeof(float));
  yy = (float *)malloc(num_particles * sizeof(float));
  zz = (float *)malloc(num_particles * sizeof(float));
  vx = (float *)malloc(num_particles * sizeof(float));
  vy = (float *)malloc(num_particles * sizeof(float));
  vz = (float *)malloc(num_particles * sizeof(float));
  phi = (float *)malloc(num_particles * sizeof(float));
  pid = (int64_t *)malloc(num_particles * sizeof(int64_t));
  mask = (uint16_t *)malloc(num_particles * sizeof(uint16_t));

  for (uint64_t i = 0; i < num_particles; i++) {
    xx[i] = (float)i;
    yy[i] = (float)i;
    zz[i] = (float)i;
    vx[i] = (float)i;
    vy[i] = (float)i;
    vz[i] = (float)i;
    phi[i] = (float)i;
    pid[i] = (int64_t)i;
    mask[i] = (uint16_t)myrank;
  }

  RestartIO_GLEAN *rst = new RestartIO_GLEAN();

  rst->Initialize(MPI_COMM_WORLD);

  // rst->SetPOSIX_IO_Interface(1);
  rst->SetPOSIX_IO_Interface();

  rst->CreateCheckpoint(fname, num_particles);

  rst->Write(xx, yy, zz, vx, vy, vz, phi, pid, mask);

  rst->Close();

  // Let's Read Restart File Now

  float *xx_r, *yy_r, *zz_r, *vx_r, *vy_r, *vz_r, *phi_r;
  int64_t *pid_r;
  uint16_t *mask_r;
  int64_t my_particles;

  my_particles = rst->OpenRestart(fname);

  if (my_particles != num_particles) {
    cout << " Particles Counts Do NOT MATCH " << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  rst->Read(xx_r, yy_r, zz_r, vx_r, vy_r, vz_r, phi_r, pid_r, mask_r);

  rst->Close();

  // Verify The contents
  for (uint64_t i = 0; i < num_particles; i++) {
    if ((xx[i] != xx_r[i]) || (yy[i] != yy_r[i]) || (zz[i] != zz_r[i]) ||
        (vx[i] != vx_r[i]) || (vy[i] != vy_r[i]) || (vz[i] != vz_r[i]) ||
        (phi[i] != phi_r[i]) || (pid[i] != pid_r[i]) ||
        (mask[i] != mask_r[i])) {
      cout << " Values Don't Match Index:" << i << endl;
      cout << "XX " << xx[i] << " " << xx_r[i] << " YY " << yy[i] << " "
           << yy_r[i] << endl;
      cout << "ZZ " << zz[i] << " " << zz_r[i] << " VX " << vx[i] << " "
           << vx_r[i] << endl;
      cout << "VY " << vy[i] << " " << vy_r[i] << " VZ " << vz[i] << " "
           << vz_r[i] << endl;
      cout << "PHI " << phi[i] << " " << phi_r[i] << " PID " << pid[i] << " "
           << pid_r[i] << endl;
      cout << "Mask: " << mask[i] << " " << mask_r[i] << endl;

      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (0 == myrank)
    cout << " CONTENTS VERIFIED... Success " << endl;

  rst->Finalize();

  delete rst;
  rst = 0;

  // Delete the Arrays
  free(xx);
  free(xx_r);
  free(yy);
  free(yy_r);
  free(zz);
  free(zz_r);
  free(vx);
  free(vx_r);
  free(vy);
  free(vy_r);
  free(vz);
  free(vz_r);

  free(phi);
  free(phi_r);

  free(pid);
  free(pid_r);

  free(mask);
  free(mask_r);

  MPI_Finalize();

  return 0;
}
