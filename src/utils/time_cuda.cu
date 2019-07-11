
#include "time_cuda.h"
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>


void start_cuda_timing(int time_id){
    cudaError_t error;

    cudaEvent_t * start = (cudaEvent_t *) malloc(sizeof(cudaEvent_t));
    error = cudaEventCreate(start);
    if (error != cudaSuccess){
        fprintf(stderr, "Failed to create start event (error code %s)!\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
    GTIMES_BUF[time_id] = (timestruct_t) start;
    // Record the start event
    error = cudaEventRecord(*start, NULL);
    if (error != cudaSuccess){
        fprintf(stderr, "Failed to record start event (error code %s)!\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }   
}

void stop_cuda_timing(int time_id){
    cudaEvent_t * start = (cudaEvent_t *) GTIMES_BUF[time_id]; // Retrieve start event
    cudaEvent_t stop;
    
    cudaError_t error;
    error = cudaEventCreate(&stop);
    if (error != cudaSuccess){
        fprintf(stderr, "Failed to create stop event (error code %s)!\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }   
    // Record the stop event
    error = cudaEventRecord(stop, NULL);
    if (error != cudaSuccess){
        fprintf(stderr, "Failed to record stop event (error code %s)!\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }

    // Wait for the stop event to complete
    error = cudaEventSynchronize(stop);

    if (error != cudaSuccess){
        fprintf(stderr, "Failed to synchronize on the stop event (error code %s)!\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }

    float msecTotal = 0.0f;
    error = cudaEventElapsedTime(&msecTotal, *start, stop);

    if (error != cudaSuccess){
        fprintf(stderr, "Failed to get time elapsed between events (error code %s)!\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }

    GTIMES[time_id] = msecTotal;
    free(start);
}