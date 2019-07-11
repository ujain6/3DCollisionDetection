#ifndef __TIME_STRUCT_H_
#define __TIME_STRUCT_H_

/* Timing method credited to Colin */
// Define _POSIX_C_SOURCE before including time.h so we can access the timers
#ifndef _POSIX_C_SOURCE
    #define _POSIX_C_SOURCE 201902L
#endif

#ifndef UNIT_TIME_GROUP
    #define UNIT_TIME_NS 1.0
    #define UNIT_TIME_US 1e3F
    #define UNIT_TIME_MS 1e6F
    #define UNIT_TIME_S  1e9F

    // Default unit of time set to us
    #define UNIT_TIME UNIT_TIME_MS
    // Timing Function
    #define MAXTIMESTRUCT 10
#endif

typedef char* timestruct_t;

// Buffering Global Variables
timestruct_t GTIMES_BUF [MAXTIMESTRUCT];
float        GTIMES     [MAXTIMESTRUCT];

float get_time_elapse(int time_id){
    return GTIMES[time_id];
}

#endif // __TIME_STRUCT_H_