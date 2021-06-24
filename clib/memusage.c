#include <stddef.h>
#ifdef __linux__
#include <sys/sysinfo.h>
#include <sys/resource.h>
#endif

/* returns the available memory on node in KiB */
size_t c_getMemAvail()
{
#ifdef __linux__
  struct sysinfo si;
  sysinfo(&si);
  si.freeram += si.bufferram;
  return si.freeram >> 10;
#else
  return 0;
#endif
}

/* returns the heap memory usage of the process in KiB */
size_t c_getMemUsage()
{
#ifdef __linux__
  struct rusage RU; /* heap memory usage */
  getrusage(RUSAGE_SELF, &RU);
  return RU.ru_maxrss;
#else
  return 0;
#endif
}
