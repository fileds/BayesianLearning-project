#ifndef SETS_HEADER
#define SETS_HEADER

void cfg(SEXP parents, int *configurations, int *nlevels);
void c_fast_config(int **columns, int nrow, int ncol, int *levels,
    int *configurations, int *nlevels, int offset);
SEXP c_configurations(SEXP parents, int factor, int all_levels);

void first_subset(int *work, int n, int offset);
int next_subset(int *work, int n, int max, int offset);

#endif
