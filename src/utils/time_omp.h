#ifndef _TIME_OMP_H_
#define _TIME_OMP_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "_time_struct.h"

void start_omp_timing(int time_id){
    double * st = (double *) malloc(sizeof(double));
    * st = omp_get_wtime();
    GTIMES_BUF[time_id] = (timestruct_t) st;
}

void stop_omp_timing(int time_id){
    double ed = omp_get_wtime();
    double st = *( (double *) GTIMES_BUF[time_id] );
    double elapse = ed - st;
    GTIMES[time_id] = elapse * 1e3F;
}

#endif // _TIME_OMP_H_