#ifndef _TIME_UNIX_H_
#define _TIME_UNIX_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "_time_struct.h"

void start_unix_timing(int time_id){
    struct timespec * st = (struct timespec *) malloc(sizeof(struct timespec));
    clock_gettime(CLOCK_MONOTONIC, st);
    GTIMES_BUF[time_id] = (timestruct_t) st;
}

void stop_unix_timing(int time_id){
    struct timespec * st = (struct timespec *) GTIMES_BUF[time_id]; // Retrieve start timing structure
    struct timespec ed;
    clock_gettime(CLOCK_MONOTONIC, &ed);
    double nanosec = (ed.tv_sec - st->tv_sec) * 1e9 + (ed.tv_nsec - st->tv_nsec);
    free(st);
    GTIMES[time_id] = nanosec / UNIT_TIME_MS;
}

#endif // _TIME_UNIX_H_