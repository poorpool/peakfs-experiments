/***********************************/
/*                                 */
/* MADbench2.h            Oct 2006 */
/* -----------                     */
/*                                 */
/* Busy-work re-definitions        */
/* (compile with -D IO)            */
/*                                 */
/* System-specific re-definitions  */
/* & declarations                  */
/* (compile with -D SYSTEM)        */
/*                                 */
/* Julian Borrill                  */
/* jdborrill@lbl.gov               */
/***********************************/

/* Busy work redefinitions */

#ifdef IO

int null() { return 0; }
#define blacs_get(a, b, c)                                                     \
  null();                                                                      \
  *c = 0
#define blacs_gridmap(a, b, c, d, e) null()
#define numroc(a, b, c, d, e)                                                  \
  (((*a) / (*b)) / (*e) + ((*c) < (((*a) / (*b)) % (*e)))) * (*b) +            \
      ((*c) == (((*a) / (*b)) % (*e))) * ((*a) % (*b))
#define descinit(a, b, c, d, e, f, g, h, i, j) null()
#define pdpotrf(a, b, c, d, e, f, g) busy_work(b, 3, gang1)
#define pdpotri(a, b, c, d, e, f, g) busy_work(b, 3, gang1)
#define pdtran(a, b, c, d, e, f, g, h, i, j, k, l) busy_work(a, 2, gang1)
#define pdsymv(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r)           \
  busy_work(b, 2, gang1)
#define pdgemr2d(a, b, c, d, e, f, g, h, i, j, k) busy_work(a, 2, gang1)
#define pdgemm(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s)        \
  busy_work(c, 3, gang2)
#define pdgemv(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s)        \
  busy_work(b, 2, gang2)
#define pddot(a, b, c, d, e, f, g, h, i, j, k, l)                              \
  busy_work(a, 1, gang2);                                                      \
  *b = 0.0
#define dposv(a, b, c, d, e, f, g, h) null()
#define pdlatra(a, b, c, d, e) (double)busy_work(a, 2, gang2)
#define ddot(a, b, c, d, e) (double)null()

#endif

/* System specific definitions */

#if defined BASSI

#ifndef IO
double pdlatra();
double ddot();
#endif

#elif defined JACQUARD

#define fopen64 fopen
#define fseeko64 fseeko
#define off64_t off_t
#define FILE_OFFSET_BITS 64
#define MPI_Request MPIO_Request
#define MPI_Wait MPIO_Wait

#ifndef IO
#define blacs_get blacs_get_
#define blacs_gridmap blacs_gridmap_
#define numroc numroc_
#define descinit descinit_
#define pdpotrf pdpotrf_
#define pdpotri pdpotri_
#define pdtran pdtran_
#define pdsymv pdsymv_
#define pdgemr2d pdgemr2d_
#define pdgemm pdgemm_
#define pdgemv pdgemv_
#define pddot pddot_
#define pdlatra pdlatra_
#define ddot ddot_
#define dposv dposv_
#endif

#elif defined COLUMBIA

#define fopen64 fopen
#define fseeko64 fseeko
#define off64_t off_t

#ifndef IO

#endif

#endif

#define FILE_PREFIX "/home/cyx/gekkofs_mount"
#define TRANSFER_SIZE 131072