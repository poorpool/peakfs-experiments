/*****************************/
/* MADbench2.c      Oct 2006 */
/* -----------               */
/*                           */
/* Julian Borrill            */
/* jdborrill@lbl.gov         */
/*****************************/

#include "MADbench2.h"
#include "mpi.h"
#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#define PCOUNT 8
#define TCOUNT 5
#define SLENGTH 64

/*******************/
/* Data Structures */
/*******************/

typedef struct {
  MPI_Comm comm;
  int context;
  int id;
  int no_pe;
  int no_pe_row;
  int no_pe_col;
  int my_pe;
  int my_pe_row;
  int my_pe_col;
} GANG;

typedef struct {
  int my_no_row;
  int my_no_col;
  int my_no_elm;
  int offcount;
  int desc[9];
} MATRIX;

typedef struct {
  int my_no_row;
  int desc[9];
} VECTOR;

/*************/
/* Functions */
/*************/

void initialize(int, char **);
void define_gang(int, GANG *, MATRIX *, VECTOR *);

void build_S();
void build_dSdC(double *, double *, double *, double *, double *, double *, int,
                int);

void invert_D();

void calc_W();
void rere_dSdC(double *, int);
void sowr_W(double *, double *, double *, double *, int);
void remap(double *, double *, int);
void fill(double *, int, int, MATRIX, double *);

void calc_dC();
void get_W(int, double *);
void dbydC(int, double *, int, double *, double *, double *);

void finalize();

void io_distmatrix(double *, GANG, MATRIX, int, char *);

void io_resync(char *);
void error_check(const char *, const char *, int);
void report_time();
int busy_work(int *, int, GANG);
double checksum(double *, int);

/********************/
/* Global Variables */
/********************/

int no_pix, no_bin, no_gang, sblocksize, fblocksize, r_mod, w_mod;
const char *IOMETHOD, *IOMODE, *FILETYPE, *REMAP;
double BWEXP = -1.0;

int no_pe, my_pe;
GANG gang1, gang2;
MPI_Comm g2pe_comm;
MPI_Request mpi_rrequest, mpi_wrequest;
MPI_Status mpi_status;
MPI_Info mpi_info;
MATRIX pp_matrix1, pp_matrix2;
VECTOR p_vector1, p_vector2;
double *Vspace, *Mspace;
char filename[SLENGTH], string[SLENGTH], fn_string[SLENGTH];
// FILE *df;
int fd;
MPI_File dfh;

double timer;
double myt[TCOUNT];

static int i0 = 0, i1 = 1;
static int binwidth = 10;
static double binweight = 1.0;
static int null_desc[] = {0, -1, 0, 0, 0, 0, 0, 0, 0};
static double d0 = 0.0, d1 = 1.0;
static char lo = 'L', no = 'N', tr = 'T';

/**********************************************************************************************************************************/

int main(int argc, char **argv) {

  initialize(argc, argv);
  PMPI_Barrier(MPI_COMM_WORLD);

  build_S();
  PMPI_Barrier(MPI_COMM_WORLD);

#ifndef IO
  invert_D();
  PMPI_Barrier(MPI_COMM_WORLD);
#endif

  calc_W();
  PMPI_Barrier(MPI_COMM_WORLD);

  calc_dC();
  PMPI_Barrier(MPI_COMM_WORLD);

  finalize();
  return 0;
}

/**********************************************************************************************************************************/

void initialize(int argc, char **argv) {
  int b, r, n, mpeak, gfb;
  int p, parameter[PCOUNT];
  char *BWEXPstring;
  struct stat buf;
  // FILE *f;

  /* Initialize MPI */

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &no_pe);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_pe);

  timer = MPI_Wtime();
  for (n = 0; n < TCOUNT; n++)
    myt[n] = 0.0;
  MPI_Pcontrol(1, "infin");

  /* Initialize run parameters */

  if (my_pe == 0) {
    error_check("insuffient arguments", "command line", argc >= PCOUNT);
    for (p = 0; p < PCOUNT; p++)
      parameter[p] = atoi(argv[p]);
  }
  MPI_Bcast(parameter, PCOUNT, MPI_INT, 0, MPI_COMM_WORLD);

  no_pix = parameter[1];
  no_bin = parameter[2];
  no_gang = parameter[3];
  sblocksize = parameter[4];
  fblocksize = parameter[5];
  r_mod = parameter[6];
  w_mod = parameter[7];

  if (my_pe == 0) {
#ifdef IO
    fprintf(stdout, "\nMADbench 2.0 IO-mode\n");
#else
    fprintf(stdout, "\nMADbench 2.0\n");
#endif
  }

  if (my_pe == 0)
    fprintf(stdout,
            "no_pe = %d  no_pix = %d  no_bin = %d  no_gang = %d  sblocksize = "
            "%d  fblocksize = %d  r_mod = %d  w_mod = %d\n",
            no_pe, no_pix, no_bin, no_gang, sblocksize, fblocksize, r_mod,
            w_mod);

  /* Validate run parameters */

  error_check("non-square no_pe", "mpi_init",
              (int)pow(floor(sqrt((double)no_pe)), 2.0) == no_pe);
  error_check("no_pe indivisible by no_gang", "command line",
              no_pe % no_gang == 0);
  error_check("no_bin indivisible by no_gang", "command line",
              no_bin % no_gang == 0);
  error_check("non-square no_gang", "command line",
              (int)pow(floor(sqrt((double)no_gang)), 2.0) == no_gang);
  error_check("sblocksize too large for no_pix & no_pe", "command line",
              ceil((double)no_pix / (double)sblocksize) >= sqrt((double)no_pe));
  error_check("fblocksize not an integer number of doubles", "command line",
              fblocksize % sizeof(double) == 0);
  error_check("no_pe/no_gang indivisible by r_mod", "command line",
              (no_pe / no_gang) % r_mod == 0);
  error_check("no_pe/no_gang indivisible by w_mod", "command line",
              (no_pe / no_gang) % w_mod == 0);

  /* Initialize & validate run environment */

  IOMETHOD = getenv("IOMETHOD");
  if (IOMETHOD == NULL) {
    IOMETHOD = (char *)malloc(SLENGTH * sizeof(char));
    IOMETHOD = "POSIX";
  }
  IOMODE = getenv("IOMODE");
  if (IOMODE == NULL) {
    IOMODE = (char *)malloc(SLENGTH * sizeof(char));
    IOMODE = "SYNC";
  }
  FILETYPE = getenv("FILETYPE");
  if (FILETYPE == NULL) {
    FILETYPE = (char *)malloc(SLENGTH * sizeof(char));
    FILETYPE = "UNIQUE";
  }
  REMAP = getenv("REMAP");
  if (REMAP == NULL) {
    REMAP = (char *)malloc(SLENGTH * sizeof(char));
    REMAP = "CUSTOM";
  }
#ifdef IO
  BWEXPstring = getenv("BWEXP");
  if (BWEXPstring != NULL)
    BWEXP = atof(BWEXPstring);
#endif

  if (my_pe == 0) {
    fprintf(stdout, "IOMETHOD = %s  IOMODE = %s  FILETYPE = %s  REMAP = %s",
            IOMETHOD, IOMODE, FILETYPE, REMAP);
#ifdef IO
    if (!(BWEXP < 0.0))
      fprintf(stdout, "  BWEXP = %.2f", BWEXP);
#endif
    fprintf(stdout, "\n\n");
    fflush(stdout);
  }
  PMPI_Barrier(MPI_COMM_WORLD);

  /* Initialize gang1 */

  gang1.id = 0;
  gang1.comm = MPI_COMM_WORLD;
  MPI_Comm_size(gang1.comm, &gang1.no_pe);
  MPI_Comm_rank(gang1.comm, &gang1.my_pe);
  define_gang(1, &gang1, &pp_matrix1, &p_vector1);

  /* Initialize gang2 */

  gang2.no_pe = gang1.no_pe / no_gang;
  gang2.id = gang1.my_pe / gang2.no_pe;
  gang2.my_pe = gang1.my_pe % gang2.no_pe;
  MPI_Comm_split(gang1.comm, gang2.id, gang2.my_pe, &gang2.comm);
  define_gang(no_gang, &gang2, &pp_matrix2, &p_vector2);

  /* Define gang2 processor number communicator */

  MPI_Comm_split(gang1.comm, gang2.my_pe, gang2.id, &g2pe_comm);

  /* Set matrix1/matrix2 data file double counts */

  if (my_pe == 0) {
    mpeak = (no_gang * pp_matrix1.my_no_elm > pp_matrix2.my_no_elm)
                ? no_gang * pp_matrix1.my_no_elm
                : pp_matrix2.my_no_elm;
    fblocksize /= sizeof(double);
    gfb = mpeak / (no_gang * fblocksize);
    if (mpeak % (no_gang * fblocksize) > 0)
      gfb++;
    pp_matrix2.offcount = gfb * no_gang * fblocksize;
    pp_matrix1.offcount = gfb * fblocksize;
  }
  MPI_Bcast(&pp_matrix1.offcount, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&pp_matrix2.offcount, 1, MPI_INT, 0, MPI_COMM_WORLD);

  /* Allocate & initialize persistent matrix1 & vector2 dataspaces */

  error_check("malloc", "Mspace",
              (Mspace = (double *)malloc(pp_matrix1.my_no_elm *
                                         sizeof(double))) != NULL);
  error_check("malloc", "Vspace",
              (Vspace = (double *)malloc(3 * p_vector2.my_no_row *
                                         sizeof(double))) != NULL);
  for (n = 0; n < pp_matrix1.my_no_elm; n++)
    Mspace[n] = 0.0;

  /* Initialize files subdirectory & open unique or shared file(s) */

  if (strcmp(IOMODE, "NONE") != 0) {

    for (n = 0; n < no_pe; n++) {
      if (my_pe == n && stat(FILE_PREFIX, &buf) != 0)
        mkdir(FILE_PREFIX, S_IRWXU);
      PMPI_Barrier(MPI_COMM_WORLD);
    }

    if (strcmp(FILETYPE, "UNIQUE") == 0)
      sprintf(filename, FILE_PREFIX "/data_%d", my_pe);
    else
      sprintf(filename, FILE_PREFIX "/data");

    if (strcmp(IOMETHOD, "POSIX") == 0)
      error_check("open", filename,
                  (fd = open(filename, O_RDWR | O_CREAT,
                             S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) != 0);
    else if (strcmp(IOMETHOD, "MPI") == 0) {
#ifdef IBM
      MPI_Info_create(&mpi_info);
      MPI_Info_set(mpi_info, "IBM_largeblock_io", "true");
#else
      mpi_info = MPI_INFO_NULL;
#endif
      error_check("MPI_File_open", filename,
                  MPI_File_open((strcmp(FILETYPE, "UNIQUE") == 0)
                                    ? MPI_COMM_SELF
                                    : MPI_COMM_WORLD,
                                filename, MPI_MODE_CREATE | MPI_MODE_RDWR,
                                mpi_info, &dfh) == 0);
    }
  }

  /* Finish timing */

  MPI_Pcontrol(-1, "infin");
}

void define_gang(int no_gang, GANG *gang_ptr, MATRIX *pp_matrix_ptr,
                 VECTOR *p_vector_ptr) {
  int ng, np, context, info;
  int *map;

  /* Define square processor array */

  gang_ptr->no_pe_col = gang_ptr->no_pe_row =
      (int)sqrt((double)gang_ptr->no_pe);

  /* Build blacs grids gang by gang */

  map = (int *)malloc(gang_ptr->no_pe * sizeof(int));
  for (ng = 0; ng < no_gang; ng++) {
    for (np = 0; np < gang_ptr->no_pe; np++)
      map[np] = ng * gang_ptr->no_pe + np;
    blacs_get(&i0, &i0, &context);
    blacs_gridmap(&context, map, &gang_ptr->no_pe_row, &gang_ptr->no_pe_row,
                  &gang_ptr->no_pe_col);
    if (gang_ptr->id == ng)
      gang_ptr->context = context;
  }
  free(map);
  gang_ptr->my_pe_row = (gang_ptr->my_pe) % (gang_ptr->no_pe_row);
  gang_ptr->my_pe_col = (gang_ptr->my_pe) / (gang_ptr->no_pe_row);

  /* Build distributed matrix/vector descriptors */

  p_vector_ptr->my_no_row = pp_matrix_ptr->my_no_row = numroc(
      &no_pix, &sblocksize, &gang_ptr->my_pe_row, &i0, &gang_ptr->no_pe_row);
  pp_matrix_ptr->my_no_col = numroc(&no_pix, &sblocksize, &gang_ptr->my_pe_col,
                                    &i0, &gang_ptr->no_pe_col);
  pp_matrix_ptr->my_no_elm =
      pp_matrix_ptr->my_no_row * pp_matrix_ptr->my_no_col;
  descinit(pp_matrix_ptr->desc, &no_pix, &no_pix, &sblocksize, &sblocksize, &i0,
           &i0, &gang_ptr->context, &pp_matrix_ptr->my_no_row, &info);
  descinit(p_vector_ptr->desc, &no_pix, &i1, &sblocksize, &sblocksize, &i0, &i0,
           &gang_ptr->context, &p_vector_ptr->my_no_row, &info);
}

/**********************************************************************************************************************************/

void build_S() {
  int n, b, lmin, lmax;
  double pi, two_pi;
  double *S, *dSdCb, *LP_lminus1, *LP_l, *wbuffer, *temp, *ra, *dec;

  strcpy(fn_string, "S");

  /* Timing */

  timer = MPI_Wtime();
  for (n = 0; n < TCOUNT; n++)
    myt[n] = 0.0;
  MPI_Pcontrol(1, fn_string);

  /* Assign persistent matrix dataspaces */

  S = Mspace;

  /* Allocate matrix dataspaces */

  error_check("malloc", "dSdCb",
              (dSdCb = (double *)malloc(pp_matrix1.my_no_elm *
                                        sizeof(double))) != NULL);
  error_check("malloc", "LP_lminus1",
              (LP_lminus1 = (double *)malloc(pp_matrix1.my_no_elm *
                                             sizeof(double))) != NULL);
  error_check(
      "malloc", "LP_l",
      (LP_l = (double *)malloc(pp_matrix1.my_no_elm * sizeof(double))) != NULL);
  error_check("malloc", "wbuffer",
              (wbuffer = (double *)malloc(pp_matrix1.my_no_elm *
                                          sizeof(double))) != NULL);

  /* Initialize legendre polynomials, ra & dec */

  for (n = 0; n < pp_matrix1.my_no_elm; n++) {
    LP_lminus1[n] = 0.0;
    LP_l[n] = 1.0;
  }
  pi = 4.0 * atan(1.0);
  two_pi = 2.0 * pi;
  srand48(0);
  ra = (double *)malloc(no_pix * sizeof(double));
  dec = (double *)malloc(no_pix * sizeof(double));
  for (n = 0; n < no_pix; n++) {
    ra[n] = two_pi * drand48();
    dec[n] = pi * (drand48() - 0.5);
  }

  /* Build & write dSdCb matrices, accumulating S */

  for (b = 0, lmax = 1; b <= no_bin; b++) {

    /* Build */

    if (b < no_bin) {
      lmin = lmax;
      lmax += binwidth;
#ifdef IO
      busy_work(&no_pix, 2, gang1);
#else
      build_dSdC(S, dSdCb, LP_lminus1, LP_l, ra, dec, lmin, lmax);
#endif
    }

    /* Resynchronize */

    if (b > 0) {
      if (strcmp(IOMODE, "ASYNC") == 0)
        io_resync("w");
    }

    /* Write */

    if (b < no_bin) {
      memcpy(wbuffer, dSdCb, pp_matrix1.my_no_elm * sizeof(double));
      io_distmatrix(wbuffer, gang1, pp_matrix1, b, "w");
    }
  }

  /* Tidy up */

  free(dec);
  free(ra);
  free(wbuffer);
  free(LP_l);
  free(LP_lminus1);
  free(dSdCb);

  MPI_Pcontrol(-1, fn_string);
  report_time();
}

void build_dSdC(double *S, double *dSdCb, double *LP_lminus1, double *LP_l,
                double *ra, double *dec, int lmin, int lmax) {
  int n, d, c, r, pc, pr, l;
  double inv_four_pi, norm, cos_chi, p0, p1, p2;
  double vc[3], vr[3];
  double *one_over_l, *two_l_plus_1;

  /* Initialize constants & dSdCb */

  inv_four_pi = 1.0 / (16.0 * atan(1.0));
  one_over_l = (double *)malloc(lmax * sizeof(double));
  two_l_plus_1 = (double *)malloc(lmax * sizeof(double));
  for (l = lmin; l < lmax; l++) {
    one_over_l[l] = 1.0 / (double)l;
    two_l_plus_1[l] = (double)(2 * l + 1);
  }
  for (n = 0; n < pp_matrix1.my_no_elm; n++)
    dSdCb[n] = 0.0;

  /* Build dSdCb */

  for (c = 0, n = 0; c < pp_matrix1.my_no_col; c++) {
    pc = ((c / sblocksize) * gang1.no_pe_col + gang1.my_pe_col) * sblocksize +
         c % sblocksize;
    vc[0] = vc[1] = cos(dec[pc]);
    vc[0] *= cos(ra[pc]);
    vc[1] *= sin(ra[pc]);
    vc[2] = sin(dec[pc]);
    for (d = 0, norm = 0.0; d < 3; d++)
      norm += vc[d] * vc[d];
    for (d = 0, norm = 1.0 / sqrt(norm); d < 3; d++)
      vc[d] *= norm;
    for (r = 0; r < pp_matrix1.my_no_row; r++, n++) {
      pr = ((r / sblocksize) * gang1.no_pe_row + gang1.my_pe_row) * sblocksize +
           r % sblocksize;
      vr[0] = vr[1] = cos(dec[pr]);
      vr[0] *= cos(ra[pr]);
      vr[1] *= sin(ra[pr]);
      vr[2] = sin(dec[pr]);
      for (d = 0, norm = 0.0; d < 3; d++)
        norm += vr[d] * vr[d];
      for (d = 0, norm = 1.0 / sqrt(norm); d < 3; d++)
        vr[d] *= norm;
      for (d = 0, cos_chi = 0.0; d < 3; d++)
        cos_chi += vc[d] * vr[d];
      p0 = LP_lminus1[n];
      p1 = LP_l[n];
      for (l = lmin; l < lmax; l++) {
        p2 = 2.0 * cos_chi * p1 - p0 - (cos_chi * p1 - p0) * one_over_l[l];
        dSdCb[n] += two_l_plus_1[l] * p2;
        p0 = p1;
        p1 = p2;
      }
      dSdCb[n] *= inv_four_pi;
      LP_lminus1[n] = p0;
      LP_l[n] = p1;
    }
  }

  /* Accumulate S */

  for (n = 0; n < pp_matrix1.my_no_elm; n++)
    S[n] += binweight * dSdCb[n];

  /* Tidy up */

  free(two_l_plus_1);
  free(one_over_l);
}

/**********************************************************************************************************************************/

void invert_D() {
  int c, r, pc, pr, n, info;
  double *D, *invD, *invDT;

  strcpy(fn_string, "D");

  /* Timing */

  timer = MPI_Wtime();
  for (n = 0; n < TCOUNT; n++)
    myt[n] = 0.0;
  MPI_Pcontrol(1, fn_string);

  /* Assign persistent matrix1 dataspace */

  invD = D = Mspace;

  /* Allocate/initialize matrix1 dataspaces */

  error_check("malloc", "invDT",
              (invDT = (double *)malloc(pp_matrix1.my_no_elm *
                                        sizeof(double))) != NULL);
  memset(invDT, 0, pp_matrix1.my_no_elm * sizeof(double));

  /* Add white noise */

  for (c = 0, n = 0; c < pp_matrix1.my_no_col; c++) {
    pc = ((c / sblocksize) * gang1.no_pe_col + gang1.my_pe_col) * sblocksize +
         c % sblocksize;
    for (r = 0; r < pp_matrix1.my_no_row; r++, n++) {
      pr = ((r / sblocksize) * gang1.no_pe_row + gang1.my_pe_row) * sblocksize +
           r % sblocksize;
      if (pr == pc)
        D[n] += 1.0;
    }
  }

  /* Invert D */

  pdpotrf(&lo, &no_pix, D, &i1, &i1, pp_matrix1.desc, &info);
  pdpotri(&lo, &no_pix, D, &i1, &i1, pp_matrix1.desc, &info);

  /* Complete upper triangle of invD */

  pdtran(&no_pix, &no_pix, &d1, invD, &i1, &i1, pp_matrix1.desc, &d0, invDT,
         &i1, &i1, pp_matrix1.desc);
  for (c = 0, n = 0; c < pp_matrix1.my_no_col; c++) {
    pc = ((c / sblocksize) * gang1.no_pe_col + gang1.my_pe_col) * sblocksize +
         c % sblocksize;
    for (r = 0; r < pp_matrix1.my_no_row; r++, n++) {
      pr = ((r / sblocksize) * gang1.no_pe_row + gang1.my_pe_row) * sblocksize +
           r % sblocksize;
      if (pr < pc)
        invD[n] = invDT[n];
    }
  }

  /* Tidy up */

  free(invDT);

  /* Timing */

  MPI_Pcontrol(-1, fn_string);
  report_time();
}

/**********************************************************************************************************************************/

void calc_W() {
  int b, g, n;
  double *invD1;
  double *invD2, *dSdCb, *Wb, *rbuffer, *wbuffer, *grbuffer, *gmbuffer;
  double *d, *z;
  double t;

  strcpy(fn_string, "W");

  /* Timing */

  timer = MPI_Wtime();
  for (n = 0; n < TCOUNT; n++)
    myt[n] = 0.0;
  MPI_Pcontrol(1, fn_string);

  /* Assign persistent matrix1 dataspace */

  invD1 = Mspace;

  /* Allocate matrix2 dataspaces */

  if (no_gang > 1) {
    error_check("malloc", "invD2",
                (invD2 = (double *)malloc(pp_matrix2.my_no_elm *
                                          sizeof(double))) != NULL);
    error_check("malloc", "grbuffer",
                (grbuffer = (double *)malloc(pp_matrix1.my_no_elm *
                                             sizeof(double))) != NULL);
    error_check("malloc", "gmbuffer",
                (gmbuffer = (double *)malloc(pp_matrix2.my_no_elm *
                                             sizeof(double))) != NULL);
  } else
    invD2 = invD1;

  error_check("malloc", "dSdCb",
              (dSdCb = (double *)malloc(pp_matrix2.my_no_elm *
                                        sizeof(double))) != NULL);
  error_check("malloc", "Wb",
              (Wb = (double *)malloc(pp_matrix2.my_no_elm * sizeof(double))) !=
                  NULL);
  error_check("malloc", "rbuffer",
              (rbuffer = (double *)malloc(pp_matrix2.my_no_elm *
                                          sizeof(double))) != NULL);
  error_check("malloc", "wbuffer",
              (wbuffer = (double *)malloc(pp_matrix2.my_no_elm *
                                          sizeof(double))) != NULL);

#ifndef IO
  /* Remap invD1 to all gangs */

  if (no_gang > 1) {
    remap(invD1, invD2, 0);
    MPI_Bcast(invD2, pp_matrix2.my_no_elm, MPI_DOUBLE, 0, g2pe_comm);
  }

  /* Assign vector dataspaces */

  d = Vspace;
  z = Vspace + pp_matrix2.my_no_row;

  /* Initialize data vector & calculate z = D^{-1} d */

  for (n = 0; n < p_vector2.my_no_row; n++)
    d[n] = 1.0;
  pdsymv(&lo, &no_pix, &d1, invD2, &i1, &i1, pp_matrix2.desc, d, &i1, &i1,
         p_vector2.desc, &i1, &d0, z, &i1, &i1, p_vector2.desc, &i1);
#endif

  for (b = 0; b <= no_bin + no_gang; b += no_gang) {

    /* Resynchronize reading */

    if (b > 0 && b <= no_bin) {
      if (no_gang == 1 && strcmp(IOMODE, "ASYNC") == 0)
        io_resync("r");
      memcpy(dSdCb, rbuffer, pp_matrix2.my_no_elm * sizeof(double));
    }

    /* Read */

    if (b < no_bin) {

      if (no_gang == 1) {

        io_distmatrix(rbuffer, gang1, pp_matrix1, b, "r");

      } else {

        for (g = 0; g <= no_gang; g++) {

          /* Resynchronize reading */

          if (g > 0) {
            if (strcmp(IOMODE, "ASYNC") == 0)
              io_resync("r");
            memcpy(gmbuffer, grbuffer, pp_matrix1.my_no_elm * sizeof(double));
          }

          /* Read */
          // printf("w dist map read\n");
          if (g < no_gang)
            io_distmatrix(grbuffer, gang1, pp_matrix1, b + g, "r");

          /* Remap */

          if (g > 0) {
#ifdef IO
            busy_work(&no_pix, 2, gang1);
#else
            remap(gmbuffer, rbuffer, g - 1);
#endif
          }
        }
      }
    }

    /* Solve */

    if (b > 0 && b <= no_bin)
      pdgemm(&no, &no, &no_pix, &no_pix, &no_pix, &d1, invD2, &i1, &i1,
             pp_matrix2.desc, dSdCb, &i1, &i1, pp_matrix2.desc, &d0, Wb, &i1,
             &i1, pp_matrix2.desc);

    /* Resynchronize writing */

    if (b > no_gang)
      if (strcmp(IOMODE, "ASYNC") == 0)
        io_resync("w");

    /* Write W */

    // printf("w dist map write\n");
    if (b > 0 && b <= no_bin) {
      memcpy(wbuffer, Wb, pp_matrix2.my_no_elm * sizeof(double));
      io_distmatrix(wbuffer, gang2, pp_matrix2, (b - no_gang) / no_gang, "w");
    }
  }

  /* Tidy up */

  free(wbuffer);
  free(rbuffer);
  free(Wb);
  free(dSdCb);
  if (no_gang > 1) {
    free(gmbuffer);
    free(grbuffer);
    free(invD2);
  }

  /* Timing */

  MPI_Pcontrol(-1, fn_string);
  report_time();
}

void remap(double *data1, double *data2, int tg) {
  int max, spr, spc, sp, *tp, pe, np;
  MATRIX spm1;
  double *recv_data1;

  /* Timing */

  PMPI_Barrier(gang1.comm);
  MPI_Pcontrol(-1, fn_string);
  MPI_Pcontrol(1, "Remap");

  /* Remap */

  if (strcmp(REMAP, "SCALAPACK") == 0) {
    pdgemr2d(&no_pix, &no_pix, data1, &i1, &i1, pp_matrix1.desc, data2, &i1,
             &i1, (gang2.id == tg) ? pp_matrix2.desc : null_desc,
             pp_matrix1.desc + 1);
  } else {

    /* Allocate data1 receive buffer */

    max = numroc(&no_pix, &sblocksize, &i0, &i0, &(gang1.no_pe_row)) *
          numroc(&no_pix, &sblocksize, &i0, &i0, &(gang1.no_pe_col));
    recv_data1 = (double *)malloc(max * sizeof(double));

    /* Determine target processor for each source processor */

    tp = (int *)malloc(gang1.no_pe * sizeof(int));
    for (spc = 0, sp = 0; spc < gang1.no_pe_col; spc++) {
      for (spr = 0; spr < gang1.no_pe_row; spr++, sp++) {
        tp[sp] = tg * gang2.no_pe + (spr % gang2.no_pe_row) +
                 (spc % gang2.no_pe_col) * gang2.no_pe_row;
      }
    }

    /* Send processor info & data */

    if (tp[gang1.my_pe] != gang1.my_pe) {
      MPI_Send(&pp_matrix1, 3, MPI_INT, tp[gang1.my_pe], 0, MPI_COMM_WORLD);
      MPI_Send(data1, pp_matrix1.my_no_elm, MPI_DOUBLE, tp[gang1.my_pe], 1,
               MPI_COMM_WORLD);
    }

    /* If I am in the target gang, receive data & distribute */

    if (gang2.id == tg) {

      /* (i) From processors : whose target I am AND which are in the target
       * group AND are not me */

      for (pe = 0, np = 0; pe < gang1.no_pe; pe++) {
        if (tp[pe] == gang1.my_pe) {
          if (pe / gang2.no_pe == tg && pe != gang1.my_pe) {
            MPI_Recv(&spm1, 3, MPI_INT, pe, 0, MPI_COMM_WORLD, &mpi_status);
            MPI_Recv(recv_data1, spm1.my_no_elm, MPI_DOUBLE, pe, 1,
                     MPI_COMM_WORLD, &mpi_status);
            fill(recv_data1, pe, np, spm1, data2);
          }
          np++;
        }
      }

      /* (ii) From processors : whose target I am AND which are not in the
       * target group */

      for (pe = 0, np = 0; pe < gang1.no_pe; pe++) {
        if (tp[pe] == gang1.my_pe) {
          if (pe / gang2.no_pe != tg) {
            MPI_Recv(&spm1, 3, MPI_INT, pe, 0, MPI_COMM_WORLD, &mpi_status);
            MPI_Recv(recv_data1, spm1.my_no_elm, MPI_DOUBLE, pe, 1,
                     MPI_COMM_WORLD, &mpi_status);
            fill(recv_data1, pe, np, spm1, data2);
          }
          np++;
        }
      }

      /* (iii) From processors : whose target I am AND which are me */

      for (pe = 0, np = 0; pe < gang1.no_pe; pe++) {
        if (tp[pe] == gang1.my_pe) {
          if (pe == gang1.my_pe) {
            fill(data1, pe, np, pp_matrix1, data2);
          }
          np++;
        }
      }
    }

    free(tp);
    free(recv_data1);
  }

  /* Timing */

  PMPI_Barrier(gang1.comm);
  MPI_Pcontrol(-1, "Remap");
  MPI_Pcontrol(1, fn_string);
}

void fill(double *d1, int pe, int np, MATRIX spm1, double *d2) {
  int r12, c12, roff, coff, nrb, ncb, rb, cb, nr, nc, r, c, offset, n;

  r12 = gang1.no_pe_row / gang2.no_pe_row;
  c12 = gang1.no_pe_col / gang2.no_pe_col;
  roff = np % r12;
  coff = np / r12;
  nrb = spm1.my_no_row / sblocksize;
  ncb = spm1.my_no_col / sblocksize;
  for (cb = 0, n = 0; cb <= ncb; cb++) {
    nc = (cb == ncb) ? spm1.my_no_col % sblocksize : sblocksize;
    for (c = 0; c < nc; c++) {
      offset = ((cb * c12 + coff) * sblocksize + c) * pp_matrix2.my_no_row;
      for (rb = 0, offset += roff * sblocksize; rb <= nrb;
           rb++, offset += r12 * sblocksize) {
        nr = (rb == nrb) ? spm1.my_no_row % sblocksize : sblocksize;
        memcpy(d2 + offset, d1 + n, nr * sizeof(double));
        n += nr;
      }
    }
  }
}

/**********************************************************************************************************************************/

void calc_dC() {

  int n, br, bc, xbr, xbc, b, bmax, no_bin2, info;
  double *Wb, *WTb, *rbuffer;
  double *dC, *dLdC, *d2LdC2, *tmp;
  double x = 0.0;

  strcpy(fn_string, "C");

  /* Timing */

  timer = MPI_Wtime();
  for (n = 0; n < TCOUNT; n++)
    myt[n] = 0.0;
  MPI_Pcontrol(1, fn_string);

  /* Allocate matrix2 dataspaces */

  error_check("malloc", "Wb",
              (Wb = (double *)malloc(pp_matrix2.my_no_elm * sizeof(double))) !=
                  NULL);
  error_check("malloc", "WTb",
              (WTb = (double *)malloc(pp_matrix2.my_no_elm * sizeof(double))) !=
                  NULL);
  error_check("malloc", "rbuffer",
              (rbuffer = (double *)malloc(pp_matrix2.my_no_elm *
                                          sizeof(double))) != NULL);

  /* Allocate/initialize bin dataspaces */

  no_bin2 = no_bin * no_bin;
  error_check("malloc", "dC = dLdC",
              (dC = dLdC = (double *)malloc(no_bin * sizeof(double))) != NULL);
  error_check("malloc", "d2LdC2",
              (d2LdC2 = (double *)malloc(no_bin2 * sizeof(double))) != NULL);
  for (b = 0; b < no_bin; b++)
    dLdC[b] = 0.0;
  for (b = 0; b < no_bin2; b++)
    d2LdC2[b] = 0.0;

  /* Calculate derivatives*/

  for (b = 0, bmax = no_bin2 / no_gang, xbr = xbc = 0; b <= bmax; b++) {

    bc = no_gang * (b / no_bin);
    br = b % no_bin;

#ifdef IO
    if (b < bmax && br != bc)
      continue;
#else
    if (b < bmax && br < bc)
      continue;
#endif

    /* Resynchronize */

    if (b > 0) {
      if (strcmp(IOMODE, "ASYNC") == 0)
        if (xbr % no_gang == 0)
          io_resync("r");
      memcpy(Wb, rbuffer, pp_matrix2.my_no_elm * sizeof(double));
    }

    /* Read/receive */

    if (b < bmax) {
      get_W(br, rbuffer);
    }

    /* Calculate derivative elements */

    if (b > 0) {
#ifdef IO
      busy_work(&no_pix, 2, gang1);
#else
      dbydC(no_gang * (xbr / no_gang) + (xbr + gang2.id) % no_gang, Wb,
            xbc + gang2.id, WTb, dLdC, d2LdC2);
#endif
      xbc = bc;
      xbr = br;
    }
  }

  /* Gather dLdC & d2LdC2 from the gangs */

  if (no_gang > 1) {
    if (gang2.my_pe == 0) {
      tmp = (double *)malloc(no_bin2 * sizeof(double));
      MPI_Reduce(dLdC, tmp, no_bin, MPI_DOUBLE, MPI_SUM, 0, g2pe_comm);
      if (gang2.id == 0)
        for (b = 0; b < no_bin; b++)
          dLdC[b] = tmp[b];
      MPI_Reduce(d2LdC2, tmp, no_bin2, MPI_DOUBLE, MPI_SUM, 0, g2pe_comm);
      if (gang2.id == 0)
        for (b = 0; b < no_bin2; b++)
          d2LdC2[b] = tmp[b];
      free(tmp);
    }
  }

  /* Calculate dC */

  if (gang1.my_pe == 0) {
    dposv(&lo, &no_bin, &i1, d2LdC2, &no_bin, dLdC, &no_bin, &info);
    x = dC[0];
  }

  /* Tidy up */

  free(d2LdC2);
  free(dLdC);
  free(rbuffer);
  free(WTb);
  free(Wb);

  /* Timing */

  MPI_Pcontrol(-1, fn_string);
  report_time();
  if (my_pe == 0)
    fprintf(stdout, "\ndC[0] = %.5e\n", x);
}

void get_W(int br, double *Wb) {
  int n, to_gang, from_gang;
  double *mpi_buffer;

  if (br % no_gang == 0) {
    io_distmatrix(Wb, gang2, pp_matrix2, br / no_gang, "r");
  } else {
    error_check("malloc", "mpi_buffer",
                (mpi_buffer = (double *)malloc(pp_matrix2.my_no_elm *
                                               sizeof(double))) != NULL);
    to_gang = (gang2.id - 1 + no_gang) % no_gang;
    from_gang = (gang2.id + 1) % no_gang;
    if (gang2.id % 2 == 0) {
      MPI_Ssend(Wb, pp_matrix2.my_no_elm, MPI_DOUBLE, to_gang, br, g2pe_comm);
      MPI_Recv(mpi_buffer, pp_matrix2.my_no_elm, MPI_DOUBLE, from_gang, br,
               g2pe_comm, &mpi_status);
    } else {
      MPI_Recv(mpi_buffer, pp_matrix2.my_no_elm, MPI_DOUBLE, from_gang, br,
               g2pe_comm, &mpi_status);
      MPI_Ssend(Wb, pp_matrix2.my_no_elm, MPI_DOUBLE, to_gang, br, g2pe_comm);
    }
    memcpy(Wb, mpi_buffer, pp_matrix2.my_no_elm * sizeof(double));
    free(mpi_buffer);
  }
}

void dbydC(int br, double *Wb, int bc, double *WTb, double *dLdC,
           double *d2LdC2) {
  double dT_Wb_z, trWb, my_trWbrWbc, trWbrWbc;
  double *d, *z, *Wb_z;

  /* Assign vector dataspaces */

  d = Vspace;
  z = Vspace + p_vector2.my_no_row;
  Wb_z = Vspace + 2 * p_vector2.my_no_row;

  /* If on diagonal, calculate dLdC & transpose Wb */

  if (br == bc) {
    pdgemv(&no, &no_pix, &no_pix, &d1, Wb, &i1, &i1, pp_matrix2.desc, z, &i1,
           &i1, p_vector2.desc, &i1, &d0, Wb_z, &i1, &i1, p_vector2.desc, &i1);
    pddot(&no_pix, &dT_Wb_z, d, &i1, &i1, p_vector2.desc, &i1, Wb_z, &i1, &i1,
          p_vector2.desc, &i1);
    trWb = pdlatra(&no_pix, Wb, &i1, &i1, pp_matrix2.desc);
    if (gang2.my_pe == 0)
      dLdC[br] = 0.5 * (dT_Wb_z - trWb);
    pdtran(&no_pix, &no_pix, &d1, Wb, &i1, &i1, pp_matrix2.desc, &d0, WTb, &i1,
           &i1, pp_matrix2.desc);
  }

  /* If lower triangular, calculate d2LdC2 */

  if (br >= bc) {
    my_trWbrWbc = ddot(&pp_matrix2.my_no_elm, Wb, &i1, WTb, &i1);
    MPI_Reduce(&my_trWbrWbc, &trWbrWbc, 1, MPI_DOUBLE, MPI_SUM, 0, gang2.comm);
    if (gang2.my_pe == 0)
      d2LdC2[br + bc * no_bin] = 0.5 * trWbrWbc;
  }
}

/**********************************************************************************************************************************/

void finalize() {
  struct stat buf;
  int n;

  MPI_Pcontrol(1, "infin");

  /* Free dataspaces */

  free(Vspace);
  free(Mspace);

  /* Close & delete unique or shared file(s) & remove files subdirectory */

  if (strcmp(IOMODE, "NONE") != 0.0) {

    if (strcmp(IOMETHOD, "POSIX") == 0)
      error_check("fclose", filename, close(fd) == 0);
    else if (strcmp(IOMETHOD, "MPI") == 0)
      error_check("MPI_File_close", filename, MPI_File_close(&dfh) == 0);

    if (strcmp(FILETYPE, "SHARED") == 0) {
      if (my_pe == 0)
        error_check("unlink", filename, unlink(filename) == 0);
    } else {
      error_check("unlink", filename, unlink(filename) == 0);
    }
    PMPI_Barrier(MPI_COMM_WORLD);

    for (n = 0; n < no_pe; n++) {
      if (my_pe == n && stat(FILE_PREFIX, &buf) == 0)
        rmdir(FILE_PREFIX);
      PMPI_Barrier(MPI_COMM_WORLD);
    }
  }

  /* Timing */

  MPI_Pcontrol(-1, "infin");

  /* Exit */

  MPI_Finalize();
}

/**********************************************************************************************************************************/

void io_distmatrix(double *data, GANG gang, MATRIX matrix, int rank, char *rw) {
  int token = 0, io_mod;
  int to_pe, from_pe;
  double dt;
  off64_t offset;
  char io_string[SLENGTH];

  /* Check for IOMODE = NONE */

  if (strcmp(IOMODE, "NONE") == 0)
    return;

  /* Synchronize & start timers */

  PMPI_Barrier(gang.comm);
  MPI_Pcontrol(-1, fn_string);
  myt[0] += MPI_Wtime() - timer;
  timer = MPI_Wtime();
  sprintf(io_string, "%s_%s", fn_string, rw);
  MPI_Pcontrol(1, io_string);

  /* Start io_mod token passing */

  io_mod = (*rw == 'r') ? r_mod : w_mod;
  to_pe = gang.my_pe + 1;
  from_pe = gang.my_pe - 1;
  if (gang.my_pe % io_mod > 0)
    PMPI_Recv(&token, 1, MPI_INT, from_pe, token, gang.comm, &mpi_status);

  /* Set offset */

  offset = rank * matrix.offcount * sizeof(double);
  if (strcmp(FILETYPE, "SHARED") == 0)
    offset += (off64_t)my_pe * no_bin * pp_matrix1.offcount * sizeof(double);

  /* POSIX IO */

  if (strcmp(IOMETHOD, "POSIX") == 0) {
    // 使用 lseek 代替 fseeko64
    error_check("lseek", filename, lseek(fd, offset, SEEK_SET) != -1);

    if (*rw == 'r') {
      if (strcmp(IOMODE, "SYNC") == 0) {
        // 使用 read 代替 fread
        ssize_t read_goal = sizeof(double) * matrix.my_no_elm;
        error_check("io block size and transfer size", filename,
                    read_goal % TRANSFER_SIZE == 0);
        ssize_t bytes_read = 0;
        while (bytes_read < read_goal) {
          bytes_read += read(fd, (char *)data + bytes_read, TRANSFER_SIZE);
        }
        error_check("read", filename,
                    bytes_read == sizeof(double) * matrix.my_no_elm);
      }
    } else {
      if (strcmp(IOMODE, "SYNC") == 0) {
        // 使用 write 代替 fwrite
        ssize_t write_goal = sizeof(double) * matrix.my_no_elm;
        error_check("io block size and transfer size", filename,
                    write_goal % TRANSFER_SIZE == 0);
        ssize_t bytes_write = 0;
        while (bytes_write < write_goal) {
          bytes_write += write(fd, (char *)data + bytes_write, TRANSFER_SIZE);
        }
        fsync(fd);
        error_check("write", filename,
                    bytes_write == sizeof(double) * matrix.my_no_elm);
      }
    }
  }

  /* MPI IO */

  else if (strcmp(IOMETHOD, "MPI") == 0) {
    error_check("MPI_File_seek", filename,
                MPI_File_seek(dfh, offset, MPI_SEEK_SET) == 0);
    if (*rw == 'r') {
      if (strcmp(IOMODE, "SYNC") == 0)
        error_check("MPI_File_read", filename,
                    MPI_File_read(dfh, data, matrix.my_no_elm, MPI_DOUBLE,
                                  &mpi_status) == 0);
      else if (strcmp(IOMODE, "ASYNC") == 0)
        MPI_File_iread(dfh, data, matrix.my_no_elm, MPI_DOUBLE, &mpi_rrequest);
    } else {
      if (strcmp(IOMODE, "SYNC") == 0)
        error_check("MPI_File_write", filename,
                    MPI_File_write(dfh, data, matrix.my_no_elm, MPI_DOUBLE,
                                   &mpi_status) == 0);
      else if (strcmp(IOMODE, "ASYNC") == 0)
        MPI_File_iwrite(dfh, data, matrix.my_no_elm, MPI_DOUBLE, &mpi_wrequest);
    }
  }

  /* End io_mod token passing */

  if (to_pe % io_mod > 0)
    PMPI_Send(&token, 1, MPI_INT, to_pe, token, gang.comm);

  /* Finish timing */

  PMPI_Barrier(gang.comm);
  MPI_Pcontrol(-1, io_string);
  dt = MPI_Wtime() - timer;
  if (*rw == 'r')
    myt[2] += dt;
  else
    myt[3] += dt;
  timer = MPI_Wtime();
  MPI_Pcontrol(1, fn_string);
}

/**********************************************************************************************************************************/

void io_resync(char *rw) {
  double dt;
  char io_string[SLENGTH], mpi_error_string[SLENGTH];
  int mpi_error_code, slength = SLENGTH;

  /* Start timing */

  PMPI_Barrier(MPI_COMM_WORLD);
  MPI_Pcontrol(-1, fn_string);
  myt[0] += MPI_Wtime() - timer;
  timer = MPI_Wtime();
  sprintf(io_string, "%s_%s", fn_string, rw);
  MPI_Pcontrol(1, io_string);

  /* Re-synchronize & error check */

  if (strcmp(IOMETHOD, "MPI") == 0) {
    error_check("MPI_Wait", rw,
                MPI_Wait((*rw == 'r') ? &mpi_rrequest : &mpi_wrequest,
                         &mpi_status) == MPI_SUCCESS);
  }

  /* Finish timing */

  PMPI_Barrier(MPI_COMM_WORLD);
  MPI_Pcontrol(-1, io_string);
  dt = MPI_Wtime() - timer;
  if (*rw == 'r')
    myt[2] += dt;
  else
    myt[3] += dt;
  timer = MPI_Wtime();
  MPI_Pcontrol(1, fn_string);
}

void error_check(const char *op, const char *name, int ok) {
  if (!ok) {
    fprintf(stderr, "PE %d failed %s on %s\n", my_pe, op, name);
    MPI_Abort(MPI_COMM_WORLD, my_pe);
  }
}

void report_time() {
  double tsum[TCOUNT], tmin[TCOUNT], tmax[TCOUNT];

  myt[0] += MPI_Wtime() - timer;
  myt[4] = myt[0] + myt[1] + myt[2] + myt[3];
  MPI_Reduce(myt, tsum, TCOUNT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(myt, tmin, TCOUNT, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(myt, tmax, TCOUNT, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (my_pe == 0) {
    if (tsum[0] > 0.0)
      fprintf(stdout, "%s_cc   %10.2f   [%10.2f:%10.2f]\n", fn_string,
              tsum[0] / (double)no_pe, tmin[0], tmax[0]);
    if (tsum[1] > 0.0)
      fprintf(stdout, "%s_bw   %10.2f   [%10.2f:%10.2f]\n", fn_string,
              tsum[1] / (double)no_pe, tmin[1], tmax[1]);
    if (tsum[2] > 0.0)
      fprintf(stdout, "%s_r    %10.2f   [%10.2f:%10.2f]\n", fn_string,
              tsum[2] / (double)no_pe, tmin[2], tmax[2]);
    if (tsum[3] > 0.0)
      fprintf(stdout, "%s_w    %10.2f   [%10.2f:%10.2f]\n", fn_string,
              tsum[3] / (double)no_pe, tmin[3], tmax[3]);
    fprintf(stdout, "          -------\n%s_total%10.2f   [%10.2f:%10.2f]\n\n",
            fn_string, tsum[4] / (double)no_pe, tmin[4], tmax[4]);
    fflush(stdout);

    { // only for beegfs
      // printf("Hurry up! sudo bash -c \"echo 3 > /proc/sys/vm/drop_caches\" on
      // "
      //        "server!!!\n");
      // sleep(10);
    }
  }
}

int busy_work(int *nptr, int scaling_exponent, GANG gang) {
  long long int n, nmax;
  int m;
  double dcount, sexp, *a, b = 1.2, c = 3.4;
  char bw_string[SLENGTH];

  PMPI_Barrier(gang.comm);
  MPI_Pcontrol(-1, fn_string);
  myt[0] += MPI_Wtime() - timer;
  timer = MPI_Wtime();
  sprintf(bw_string, "%s_bw", fn_string);
  MPI_Pcontrol(1, bw_string);

  dcount = (BWEXP < 0.0) ? (double)(*nptr) : (double)(*nptr) * (double)(*nptr);
  sexp = (BWEXP < 0.0) ? (double)scaling_exponent : BWEXP;
  nmax = (long long int)(pow(dcount, sexp) /
                         (double)(gang.no_pe * 2 * sblocksize));

  error_check("malloc", "bw",
              (a = (double *)malloc(sblocksize * sizeof(double))) != NULL);
  for (m = 0; m < sblocksize; m++)
    a[m] = 0.0;
  for (n = 0; n < nmax; n++) {
    for (m = 0; m < sblocksize; m++)
      a[m] += b * c;
  }
  free(a);

  PMPI_Barrier(gang.comm);
  MPI_Pcontrol(-1, bw_string);
  myt[1] += MPI_Wtime() - timer;
  timer = MPI_Wtime();
  MPI_Pcontrol(1, fn_string);

  return (int)(a[0] / a[1]);
}

double checksum(double *M, int mmax) {
  int m;
  double x = 0.0;

  for (m = 0; m < mmax; m++)
    x += M[m];
  fprintf(stdout, "PE%d  %s  %.3e  %.3e\n", my_pe, fn_string, M[0], x);
  fflush(stdout);

  return x;
}
