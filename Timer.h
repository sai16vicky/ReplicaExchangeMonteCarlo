#ifndef TIMER_H
#define TIMER_H
#include <sys/times.h>
#include <sys/param.h>
#include <limits.h>
#include <iostream>
#include "Const.h"
#include <time.h>
#include <unistd.h>
#include <math.h>

using namespace std;

// Timer class for measuring cpu-time of program fragments */
class Timer {
  clock_t t_start_u, t_start_s;
  clock_t t_stop_u, t_stop_s;
  enum state {RUNNING,STOPPED};
  state st;
  
public:
  Timer(void);
  void start(void);
  void stop(void);
  void reset(void);
  float resolution(void); 
  // returns resolution of Timer in CPU-seconds */
  float elapsed(void);  // returns stopped time in CPU-seconds */
  float elapsedDetail(void);
};

#ifndef HZ
#define HZ 60
#endif

Timer::Timer(void) {
    t_start_u=0;
    t_start_s=0;
    t_stop_u=0;
    t_stop_s=0;
    st = STOPPED;
}

void Timer::start(void) {
    struct tms buf;
    clock_t t;
    t = times(&buf);
    t_start_u = buf.tms_utime;
    t_start_s = buf.tms_stime;
    st = RUNNING;
}

void Timer::stop(void) {
    struct tms buf;
    clock_t t;
    t = times(&buf);
    t_stop_u = buf.tms_utime;
    t_stop_s = buf.tms_stime;
    st = STOPPED;
}

void Timer::reset(void) {
    struct tms buf;
    clock_t t;
    t = times(&buf);
    t_start_u = buf.tms_utime;
    t_start_s = buf.tms_stime;
}

float Timer::elapsed(void) {
    struct tms buf;
    if (st==STOPPED)
        return(1.0*((t_stop_u+t_stop_s)-(t_start_u+t_start_s))/CLOCKS_PER_SEC);
    else {
        clock_t t;
        t = times(&buf);
        //cout << "user time: " << (1.0*((buf.tms_utime)-(t_start_u))/CLK_TCK) << endl;
        //cout << "sys  time: " << (1.0*((buf.tms_stime)-(t_start_s))/CLK_TCK) << endl;
        //cout << "wall time: " << (1.0*((buf.tms_utime+buf.tms_stime)-(t_start_u+t_start_s))/CLK_TCK) << endl << endl;
        return(1.0*((buf.tms_utime+buf.tms_stime)-(t_start_u+t_start_s))/CLOCKS_PER_SEC);
    }
}

float Timer::elapsedDetail(void) {
    struct tms buf;
    if (st==STOPPED)
        return(1.0*((t_stop_u+t_stop_s)-(t_start_u+t_start_s))/CLOCKS_PER_SEC);
    else {
        clock_t t;
        t = times(&buf);
        cout << "user time: " << (1.0*((buf.tms_utime)-(t_start_u))/CLOCKS_PER_SEC) << endl;
        cout << "sys  time: " << (1.0*((buf.tms_stime)-(t_start_s))/CLOCKS_PER_SEC) << endl;
        cout << "real time: " << (1.0*((buf.tms_utime+buf.tms_stime)-(t_start_u+t_start_s))/CLOCKS_PER_SEC) << endl << endl;
        return(1.0*((buf.tms_utime+buf.tms_stime)-(t_start_u+t_start_s))/CLOCKS_PER_SEC);
    }
}



float Timer::resolution(void) {
    return(1.0/CLOCKS_PER_SEC);
}

#endif