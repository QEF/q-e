#ifndef pbecor_h
#define pbecor_h 1

#pragma acc routine seq
extern void corpbe(double, double,
	int, int, double *, double *,
	double *, double *, double *);

#pragma acc routine seq
void corpbespin(double, double, double,
	int, int, double *, double *, double *,
	double *, double *, double *, double *);

#define invpi075tothird 0.620350490899400016668006812047780
#define r2k 1.56318528359354405905402280568002

#endif
