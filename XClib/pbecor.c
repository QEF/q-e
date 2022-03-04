/* taken from K. Burke's homepage: http://dft.uci.edu/node/187
   translated by f2c (version 20090411).
   code for deriv. wrt. absolute gradient divided by abs. grad. added
*/

#include <math.h>

/* Table of constant values */

#define c_b3  .0310907
#define c_b4  .2137
#define c_b5  7.5957
#define c_b6  3.5876
#define c_b7  1.6382
#define c_b8  .49294
#define c_b9  .01554535
#define c_b10  .20548
#define c_b11  14.1189
#define c_b12  6.1977
#define c_b13  3.3662
#define c_b14  .62517
#define c_b15  .0168869
#define c_b16  .11125
#define c_b17  10.357
#define c_b18  3.6231
#define c_b19  .88026
#define c_b20  .49671

/* ---------------------------------------------------------------------- */
/* ###################################################################### */
/* ---------------------------------------------------------------------- */
/* Subroutine */ 
#pragma acc routine seq
static void gcor2(double a, double a1, double b1, double b2, double b3, double b4, double rtrs, double *gg, double *ggrs)
{
    /* Local variables */
    double q0, q1, q2, q3;

/* slimmed down version of GCOR used in PW91 routines, to interpolate */
/* LSD correlation energy, as given by (10) of */
/* J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992). */
/* K. Burke, May 11, 1996. */
    q0 = a * -2. * (a1 * rtrs * rtrs + 1.);
    q1 = a * 2. * rtrs * (b1 + rtrs * (b2 + rtrs * (b3 + b4 * rtrs)))
	    ;
    q2 = log(1. / q1 + 1.);
    *gg = q0 * q2;
    q3 = a * (b1 / rtrs + b2 * 2. + rtrs * (b3 * 3. + b4 * 4. * rtrs))
	    ;
    *ggrs = a * -2. * a1 * q2 - q0 * q3 / (q1 * (q1 + 1.));
} /* gcor2_ */


#define gamma 0.03109069086965489503494086371273
#define beta 0.06672455060314922
#define delta (beta/gamma)
#define phi 0.409240950261429630406974844266531
#define GAM 0.5198420997897463295344212145565
#define fzz (8./(9.*GAM))


/* ###################################################################### */
/* ---------------------------------------------------------------------- */
/* Subroutine */ 
#pragma acc routine seq
void corpbe(double rs, double t, int lgga, int lpot, double *ec, double *vc, double *h__, double *dvc, double *dv2rho)
{
    /* Local variables */
    double b, b2, q4, t2, q5, t4,
	    eu, pon,
	    eurs, rtrs;
    double tmp1,tmp2,tmp3;

/* ---------------------------------------------------------------------- */
/*  Official PBE correlation code. K. Burke, May 14, 1996. */
/*  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3) */
/*       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho */
/*       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE */
/*       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3) */
/*       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2) */
/*       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2 */
/*       :  UU,VV,WW, only needed for PBE potential */
/*       : lgga=flag to do gga (0=>LSD only) */
/*       : lpot=flag to do potential (0=>energy only) */
/*  output: ec=lsd correlation energy from [a] */
/*        : vcup=lsd up correlation potential */
/*        : vcdn=lsd dn correlation potential */
/*        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON */
/*        : dvcup=nonlocal correction to vcup */
/*        : dvcdn=nonlocal correction to vcdn */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* References: */
/* [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, */
/*     {\sl Generalized gradient approximation made simple}, sub. */
/*     to Phys. Rev.Lett. May 1996. */
/* [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff */
/*     construction of a generalized gradient approximation:  The PW91 */
/*     density functional}, submitted to Phys. Rev. B, Feb. 1996. */
/* [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992). */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* thrd*=various multiples of 1/3 */
/* numbers for use in LSD energy spin-interpolation formula, [c](9). */
/*      GAM= 2^(4/3)-2 */
/*      FZZ=f''(0)= 8/(9*GAM) */
/* numbers for construction of PBE */
/*      gamma=(1-log(2))/pi^2 */
/*      bet=coefficient in gradient expansion for correlation, [a](4). */
/*      eta=small number to stop d phi/ dzeta from blowing up at */
/*          |zeta|=1. */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* find LSD energy contributions, using [c](10) and Table I[c]. */
/* EU=unpolarized LSD correlation energy */
/* EURS=dEU/drs */
/* EP=fully polarized LSD correlation energy */
/* EPRS=dEP/drs */
/* ALFM=-spin stiffness, [c](3). */
/* ALFRSM=-dalpha/drs */
/* F=spin-scaling factor from [c](9). */
/* construct ec, using [c](8) */
    rtrs = sqrt(rs);
    gcor2(c_b3, c_b4, c_b5, c_b6, c_b7, c_b8, rtrs, &eu, &eurs);
/* Computing 4th power */
    *ec = eu;
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* LSD potential from [c](A1) */
/* ECRS = dEc/drs [c](A2) */
/* ECZET=dEc/dzeta [c](A3) */
/* FZ = dF/dzeta [c](A4) */
/* Computing 3rd power */
    *vc = eu - rs * eurs / 3.;
    if (lgga == 0) {
	return;
    }
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* PBE correlation energy */
/* G=phi(zeta), given after [a](3) */
/* DELT=bet/gamma */
/* B=A of [a](8) */
/* Computing 3rd power */
    pon = -eu / gamma;
    b = delta / (exp(pon) - 1.);
    b2 = b * b;
    t2 = t * t;
    t4 = t2 * t2;
    q4 = b * t2 + 1.;
    q5 = b * t2 + 1. + b2 * t4;
    *h__ = gamma * log(q4 * delta * t2 / q5 + 1.);
    if (lpot == 0) {
	return;
    }
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ENERGY DONE. NOW THE POTENTIAL. */
    tmp1 = q4/q5;
    tmp2 = b2*t4*(1.+q4)/(q5*q5);
    tmp3 = 1./(1.+delta*t2*tmp1);

    *dvc = *h__ - beta*t2 * (7./3.*tmp1 + tmp2*((b+delta)*(*vc-eu)/beta-7./3.)) * tmp3;
    *dv2rho = 0.5*beta*phi*rs*(tmp1 - tmp2) * tmp3;
} /* corpbe_ */


/* ###################################################################### */
/* ---------------------------------------------------------------------- */
/* Subroutine */ 
#pragma acc routine seq
void corpbespin(double rs, double t, double zet,
	int lgga, 
	int lpot, double *ec, double *vcup, double *vcdown,
	double *h__, double *dvcup, double *dvcdown, double *dv2rho)
{
    /* Local variables */
    double b, b2, q4, t2, q5, t4,
	    ep, eu, pon,
	    alfm, eprs, eurs, rtrs, f, z4, ecrs, fz, eczet, comm, g, g2, g3, dg;
    double tmp1,tmp2,tmp3,tmp4;
    double alfrsm;

/* ---------------------------------------------------------------------- */
/*  Official PBE correlation code. K. Burke, May 14, 1996. */
/*  INPUT: RS=SEITZ RADIUS=(3/4pi rho)^(1/3) */
/*       : ZET=RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho */
/*       : t=ABS(GRAD rho)/(rho*2.*KS*G)  -- only needed for PBE */
/*       : UU=(GRAD rho)*GRAD(ABS(GRAD rho))/(rho**2 * (2*KS*G)**3) */
/*       : VV=(LAPLACIAN rho)/(rho * (2*KS*G)**2) */
/*       : WW=(GRAD rho)*(GRAD ZET)/(rho * (2*KS*G)**2 */
/*       :  UU,VV,WW, only needed for PBE potential */
/*       : lgga=flag to do gga (0=>LSD only) */
/*       : lpot=flag to do potential (0=>energy only) */
/*  output: ec=lsd correlation energy from [a] */
/*        : vcup=lsd up correlation potential */
/*        : vcdn=lsd dn correlation potential */
/*        : h=NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON */
/*        : dvcup=nonlocal correction to vcup */
/*        : dvcdn=nonlocal correction to vcdn */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* References: */
/* [a] J.P.~Perdew, K.~Burke, and M.~Ernzerhof, */
/*     {\sl Generalized gradient approximation made simple}, sub. */
/*     to Phys. Rev.Lett. May 1996. */
/* [b] J. P. Perdew, K. Burke, and Y. Wang, {\sl Real-space cutoff */
/*     construction of a generalized gradient approximation:  The PW91 */
/*     density functional}, submitted to Phys. Rev. B, Feb. 1996. */
/* [c] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992). */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* thrd*=various multiples of 1/3 */
/* numbers for use in LSD energy spin-interpolation formula, [c](9). */
/*      GAM= 2^(4/3)-2 */
/*      FZZ=f''(0)= 8/(9*GAM) */
/* numbers for construction of PBE */
/*      gamma=(1-log(2))/pi^2 */
/*      bet=coefficient in gradient expansion for correlation, [a](4). */
/*      eta=small number to stop d phi/ dzeta from blowing up at */
/*          |zeta|=1. */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* find LSD energy contributions, using [c](10) and Table I[c]. */
/* EU=unpolarized LSD correlation energy */
/* EURS=dEU/drs */
/* EP=fully polarized LSD correlation energy */
/* EPRS=dEP/drs */
/* ALFM=-spin stiffness, [c](3). */
/* ALFRSM=-dalpha/drs */
/* F=spin-scaling factor from [c](9). */
/* construct ec, using [c](8) */
    rtrs = sqrt(rs);
    gcor2(c_b3, c_b4, c_b5, c_b6, c_b7, c_b8, rtrs, &eu, &eurs);
    gcor2(c_b9, c_b10, c_b11, c_b12, c_b13, c_b14, rtrs, &ep, &eprs);
    gcor2(c_b15, c_b16, c_b17, c_b18, c_b19, c_b20, rtrs, &alfm, &alfrsm);
/* Computing 4th power */
    z4 = zet*zet*zet*zet;
    f = (pow(1.+zet,4./3.)+pow(1.-zet,4./3.)-2.)/GAM;
    *ec = eu*(1.-f*z4) + ep*f*z4 - alfm*f*(1.-z4)/fzz;
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* LSD potential from [c](A1) */
/* ECRS = dEc/drs [c](A2) */
/* ECZET=dEc/dzeta [c](A3) */
/* FZ = dF/dzeta [c](A4) */
/* Computing 3rd power */
    ecrs = eurs*(1.-f*z4)+eprs*f*z4-alfrsm*f*(1.-z4)/fzz;
    fz = 4./3.*(pow(1.+zet,1./3.)-pow(1.-zet,1./3.))/GAM;
    eczet = 4.*pow(zet,3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu-(1.-z4)*alfm/fzz);
    comm = *ec - rs*ecrs/3. - zet*eczet;
    *vcup = comm + eczet;
    *vcdown = comm - eczet;
    if (lgga == 0) {
	return;
    }
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* PBE correlation energy */
/* G=phi(zeta), given after [a](3) */
/* DELT=bet/gamma */
/* B=A of [a](8) */
/* Computing 3rd power */
    g = 0.5*(pow(1.+zet,2./3.)+pow(1.-zet,2./3.));
    g2 = g*g;
    g3 = g2*g;
    dg = (1./3.) * (pow(1.+zet,-1./3.)-pow(1.-zet,-1./3.));
    pon = -(*ec) / (g3*gamma);
    b = delta / (exp(pon) - 1.);
    t /= g;
    b2 = b * b;
    t2 = t * t;
    t4 = t2 * t2;
    q4 = b * t2 + 1.;
    q5 = b * t2 + 1. + b2 * t4;
    *h__ = g3 * gamma * log(q4 * delta * t2 / q5 + 1.);
    if (lpot == 0) {
	return;
    }
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ENERGY DONE. NOW THE POTENTIAL. */
    tmp1 = q4/q5;
    tmp2 = b2*t4*(1.+q4)/(q5*q5);
    tmp3 = 1./(1.+delta*t2*tmp1);
    tmp4 = dg * (3.*(*h__)/g - beta*t2*g2*(2.*tmp1 - tmp2*(3.*(b+delta)*(*ec)/(beta*g3)+2.))*tmp3);

    *dvcup = *h__ - beta*g3*t2 * (7./3.*tmp1 + tmp2*((b+delta)*(*vcup-(*ec))/(beta*g3)-7./3.)) * tmp3
	+ tmp4*(1.-zet);
    *dvcdown = *h__ - beta*g3*t2 * (7./3.*tmp1 + tmp2*((b+delta)*(*vcdown-(*ec))/(beta*g3)-7./3.)) * tmp3
	- tmp4*(1.+zet);
    *dv2rho = 0.5*beta*g*phi*rs*(tmp1 - tmp2) * tmp3;
} /* corpbe_ */
