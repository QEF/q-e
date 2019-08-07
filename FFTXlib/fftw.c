#include "fftw.h"

 /* die when fatal errors occur */
void fftw_die(char* s)
{
	fprintf(stderr, "%s", s);
	exit(1);
}

