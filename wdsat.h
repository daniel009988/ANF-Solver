//
//  wdsat.h
//

#ifndef wdsat_h
#define wdsat_h

#include <stdio.h>
#include <stdbool.h>

#include "wdsat_utils.h"

bool wdsat_infer(const int_t l);
bool wdsat_solve(int_t n, char mvc_graph[1000], char thread[1000], char cubes[10000], int mode);
bool wdsat_solve_rest_XG(int_t l, int_t set_end, int_t conf[]);
#endif
