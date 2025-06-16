// Copyright 2013-2022 Edgar Costa
// See LICENSE file for license details.
// timing.h: header file for timing and profiling routines

#ifndef TIMING_H
#define TIMING_H

#include <stdio.h>
#include <time.h>
#include <sys/resource.h>
#include <sys/time.h>

typedef struct timeval timeval;

static double timevaldiff_in_seconds(timeval start, timeval end)
{
  /* Perform the carry for the later subtraction by updating start. */
  if (end.tv_usec < start.tv_usec) {
    int nsec = (start.tv_usec - end.tv_usec) / 1000000 + 1;
    start.tv_usec -= 1000000 * nsec;
    start.tv_sec += nsec;
  }
  if (end.tv_usec - start.tv_usec > 1000000) {
    int nsec = (end.tv_usec - start.tv_usec) / 1000000;
    start.tv_usec += 1000000 * nsec;
    start.tv_sec -= nsec;
  }

  return end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec)*1e-6;
}



#ifdef __APPLE__

#include <sys/time.h>

typedef timeval timestamp_type;

static void get_timestamp(timestamp_type *t)
{
  gettimeofday(t, NULL);
}

static double timestamp_diff_in_seconds(timestamp_type start,
timestamp_type end) {
    return timevaldiff_in_seconds(start, end);
}

#else

#include <time.h>

typedef struct timespec timestamp_type;

static void get_timestamp(timestamp_type *t)
{
  clock_gettime(CLOCK_REALTIME, t);
}

static double timestamp_diff_in_seconds(timestamp_type start, timestamp_type end)
{
  timestamp_type temp;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    temp.tv_sec = end.tv_sec-start.tv_sec-1;
    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  } else {
    temp.tv_sec = end.tv_sec-start.tv_sec;
    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return temp.tv_sec + 1e-9*temp.tv_nsec;
}

#endif


static inline double get_cpu_time(void){
    return (double)clock() / CLOCKS_PER_SEC;
}


typedef struct timestamp_pair {
    struct rusage usage;
    timestamp_type wall_time;
} time_pair;

static inline void timestamp_mark(time_pair *pair) {
  getrusage(RUSAGE_SELF, &pair->usage);
  get_timestamp(&pair->wall_time);
}

static inline void timestamp_report(double *utime, double *wtime, time_pair *stamp) {
  time_pair now;
  timestamp_mark(&now);
  *utime = timevaldiff_in_seconds(stamp->usage.ru_utime, now.usage.ru_utime);
  *wtime = timestamp_diff_in_seconds(stamp->wall_time, now.wall_time);
}
#endif
