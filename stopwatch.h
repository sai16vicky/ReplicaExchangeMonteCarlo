
#ifndef TBCI_STOPWATCH_H
#define TBCI_STOPWATCH_H
#include <time.h>
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#if defined (__linux__) && !defined(CLK_TCK)
# define CLK_TCK ((clock_t) sysconf (_SC_CLK_TCK))
#endif

#ifdef HAVE_LIMITS_H
# include <limits.h>
#else
# ifndef LONG_MIN
#  if defined(__WORDSIZE) && __WORDSIZE == 64
#   define LONG_MIN (-9223372036854775807L - 1L)
#  else
#   define LONG_MIN (-2147483647L - 1L)
#  endif
# endif
#endif
#if defined(__linux__) && defined(BROKEN_HZ)
# include <asm/param.h>
# define CPS (CLOCKS_PER_SEC*HZ/CLK_TCK)
#else
# define CPS CLOCKS_PER_SEC
#endif

class stopw_base
{
protected:
    double last_time;
    double total;
    const double secs_per_tick;
    const double overflow_secs;
    int running;
    
    virtual double seconds() const
    {
        return ((double) clock()) * secs_per_tick;
    }
    double adv_stopwatch()
    {
        const double secs = seconds();
        double diff = secs - last_time;
#ifdef __i386__ /* Extra precision is killing us on i386 */
        if (diff < 0.0 && -diff < secs_per_tick/4)
            diff = 0.0;
#endif
        if (overflow_secs != 0.0 && diff < 0.0)
            diff += overflow_secs;
        
        last_time = secs;
        total += diff;
        return diff;
    }
    
public:
    stopw_base(const double tick, const double over = -1.0)
    : last_time(0.0)
    , total(0.0)
    , secs_per_tick (tick)
    , overflow_secs (over == -1.0 ? -tick * (2.0 * LONG_MIN): over)
    , running(0)
    {}
    virtual ~stopw_base() {}
    double reset()
    {
        const double t = total;
        last_time = 0.0;
        total = 0.0; running = 0;
        return t;
    }
    void start()
    {
        if (!running) {
            running = 1; last_time = seconds();
        } else
            adv_stopwatch();
    }
    double stop()
    {
        if (running) {
            running = 0;
            adv_stopwatch();
        }
        return total;
    }
    double stop_d()
    {
        if (running) {
            running = 0;
            return adv_stopwatch();
        } else
            return 0.0;
    }
    double read()
    {
        if (running)
            adv_stopwatch();
        return total;
    }
    double read_d()
    {
        if (running)
            return adv_stopwatch();
        else
            return 0.0;
    }
};

class stopwatch : public stopw_base
{
public:
    stopwatch() : stopw_base(1.0/CPS)  {}
};

#ifndef unix
typedef stopwatch stopwatch_u;
typedef stopwatch stopwatch_us;
#else
//#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
//#endif
class stopwatch_u : public stopw_base
{
protected:
    virtual double seconds() const
    {
        struct tms tims;
        /*if (*/times (&tims)/* < 0) STD__ cerr << "Error reading time!" << STD__ endl*/;
        return (double)(tims.tms_utime+tims.tms_cutime) * secs_per_tick;
    };
    
public:
    stopwatch_u() : stopw_base(1000000.0 / (CPS * CLK_TCK))  {};
};

class stopwatch_us : public stopw_base
{
protected:
    virtual double seconds() const
    {
        struct tms tims;
        /*if (*/times (&tims)/* < 0) STD__ cerr << "Error reading time!" << STD__ endl*/;
        return (double)(tims.tms_utime+tims.tms_cutime+tims.tms_stime+tims.tms_cstime)
        * secs_per_tick;
    }
    
public:
    stopwatch_us() : stopw_base(1000000.0 / (CPS * CLK_TCK))  {}
};

#endif /* unix */

#ifdef unix
//#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
//#endif

//#ifdef HAVE_GETTIMEOFDAY
class stopwatch_e : public stopw_base
{
protected:
    double seconds() const
    {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return (double)tv.tv_sec + (double)tv.tv_usec * 0.000001;
    };
    
public:
    stopwatch_e() : stopw_base(0.000001, 0)  {}
};
//#endif /* HAVE_GETTIMEOFDAY */

#else /* unix */
# ifdef _MSC_VER

#include <sys/timeb.h>

class stopwatch_e : public stopw_base
{
protected:
    double seconds() const
    {
        struct _timeb tb;
        _ftime(&tb);
        return (double)tb.time + (double)tb.millitm * 0.001;
    };
    
public:
    stopwatch_e() : stopw_base(0.001)  {}
};

# else /* _MSC_VER */

class stopwatch_e : public stopw_base
{
protected:
    double seconds() const
    {
        //time_t tm;
        return (double) time(0);
    };
    
public:
    stopwatch_e() : stopw_base(1.0)  {}
};

# endif /* _MSC_VER */

#endif /* unix */

#endif /* TBCI_STOPWATCH_H */
