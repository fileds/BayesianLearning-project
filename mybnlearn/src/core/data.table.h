#ifndef DATA_TABLE_HEADER
#define DATA_TABLE_HEADER

/* flags for the columns (variables) of a data table, stored in the metadata. */
typedef struct {

  unsigned int own      : 1;  /* this column belongs with the struct. */
  unsigned int discrete : 1;  /* this column is discrete data. */
  unsigned int gaussian : 1;  /* this column is Gaussian data. */
  unsigned int complete : 1;  /* this column contains no missing data points. */
  unsigned int fixed    : 1;  /* this column should never be (re)moved. */
  unsigned int drop     : 1;  /* this column is to be removed. */
  unsigned int padding  : 2;  /* pad to 1 byte. */

} flags;

/* columns (variables) metadata common to all types of data tables. */
typedef struct {

  int nobs;             /* number of observations. */
  int ncols;            /* number of columns. */
  const char **names;   /* column names (optional). */
  flags *flag;

} meta;

/* data table for discrete data. */
typedef struct {

  meta m;               /* metadata. */
  int **col;            /* pointers to the discrete columns. */
  int *nlvl;            /* number of levels of the discrete columns. */

} ddata;

/* data table for gaussian data. */
typedef struct {

  meta m;               /* metadata. */
  double **col;         /* pointers to the continuous columns. */
  double *mean;         /* means of the continuous columns (optional). */

} gdata;

/* data table for conditional gaussian data. */
typedef struct {

  meta m;               /* metadata. */
  int **dcol;           /* pointers to the discrete columns. */
  double **gcol;        /* pointers to the continuous columns. */
  int *nlvl;            /* number of levels of the discrete columns. */
  int ndcols;           /* number of discrete columns. */
  int ngcols;           /* number of continuous columns. */
  int *map;             /* mapping between the original column position and the
                         * positions of the columns in the data structure. */

} cgdata;

void meta_init_flags(meta *m, int offset, SEXP complete, SEXP fixed);
void meta_copy_names(meta *m, int offset, SEXP df);
void meta_copy(meta *src, meta *dest);
void FreeMETA(meta *m);

ddata ddata_from_SEXP(SEXP df, int offset);
ddata empty_ddata(int nobs, int ncols);
void print_ddata(ddata dt);
void ddata_drop_flagged(ddata *dt, ddata *copy);
void ddata_subset_columns(ddata *dt, ddata *copy, int *ids, int nids);
void FreeDDT(ddata dt);

gdata gdata_from_SEXP(SEXP df, int offset);
gdata new_gdata(int nobs, int ncols);
gdata empty_gdata(int nobs, int ncols);
void gdata_cache_means(gdata *dt, int offset);
void gdata_move_column(gdata *dt, gdata *copy, int i, int j);
void print_gdata(gdata dt);
void gdata_drop_flagged(gdata *dt, gdata *copy);
void gdata_subset_columns(gdata *dt, gdata *copy, int *ids, int nids);
void gdata_incomplete_cases_range(gdata *dt, bool *indicator, int col_start,
    int col_end);
#define gdata_incomplete_cases(gdata, indicator, offset) \
  gdata_incomplete_cases_range(gdata, indicator, offset, (*(gdata)).m.ncols - 1)
void gdata_subsample_by_logical(gdata *dt, gdata *copy, bool *indicators,
    int offset);
void FreeGDT(gdata dt);

cgdata cgdata_from_SEXP(SEXP df, int doffset, int goffset);
cgdata new_cgdata(int nobs, int dcols, int gcols);
cgdata empty_cgdata(int nobs, int dcols, int gcols);
void print_cgdata(cgdata dt);
void cgdata_drop_flagged(cgdata *dt, cgdata *copy);
void cgdata_subset_columns(cgdata *dt, cgdata *copy, int *ids, int nids);
void cgdata_incomplete_cases(cgdata *dt, bool *indicator, int doffset,
    int goffset);
void cgdata_subsample_by_logical(cgdata *dt, cgdata *copy, bool *indicators,
    int doffset, int goffset);
void FreeCGDT(cgdata);

#endif
