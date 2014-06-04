
#if defined(__TRAP_SIGUSR1) || defined(__TERMINATE_GRACEFULLY)
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>

int init_signal(int signum, void (*new_handler)(int))
{
  static struct sigaction action;

  action.sa_handler   = new_handler;
  // Don't block anything.
  // Not sure if it's the correct behavior (or even if there is one)
  sigemptyset(&action.sa_mask);
  // This will probably make MPI happy
  action.sa_flags     = SA_RESTART;
  
  return sigaction(signum, &action, NULL);
}


int init_signal_USR1(void (*new_handler)(int))
{
  return init_signal(SIGUSR1, new_handler);
}

int init_TERMINATE_GRACEFULLY(void (*new_handler)(int))
{
  int fail;
  fail = init_signal(SIGTERM, new_handler);
  if(fail) return -SIGTERM;
  fail = init_signal(SIGINT,  new_handler);
  if(fail) return -SIGINT;
  fail = init_signal(SIGUSR1,  new_handler);
  if(fail) return -SIGUSR1;
  fail = init_signal(SIGUSR2,  new_handler);
  if(fail) return -SIGUSR2;
  fail = init_signal(SIGXCPU,  new_handler);
  if(fail) return -SIGXCPU;
  return 0;
}

#else
  void dummy ( ) { }
#endif
