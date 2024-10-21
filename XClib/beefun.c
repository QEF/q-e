#define _XOPEN_SOURCE 500

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "beefleg.h"

#include "pbecor.h"

#ifdef _WIN32
#define srandom srand
#endif

// evaluate bee exchange energy and its derivatives de/drho and ( de/d|grad rho| ) / |grad rho|
#pragma acc routine seq
void beefx_(double *r, double *g, double *e, double *dr, double *dg, int *addlda)
{
    double s2,t,r43,r83,s,sx,dx,fx,dl,dfx;
    const int n=nmax;
    
    double L[nmax]={1.};
    double dL[nmax]={0.,1.};

    switch(beeftype) {
    case 0: //BEEF-vdW xc
    r43 = pow(*r, 4./3.);
    r83 = r43*r43;
    sx = r2e * r43;
    dx = 4./3. * sx / (*r);

    s2 = *g*pix / r83;
    s = sqrt(s2);
    t = 2.*s2/(4.+s2)-1.;

    if(beeforder==-1)
    {
	calclegdleg(t);
	if(!(*addlda)){
	    fx = ddot1(mi,L,n) - 1.;
	    }
	else{
	    fx = ddot1(mi,L,n);
	   }

	dl = ddot1(mi,dL,n);
	dfx = dl*( 4.*s / (4.+s2) - 4.*s2*s/sq(4.+s2) );
	*dr = dx*fx - 4./3.*s2/(s*(*r))*sx*dfx;
	*dg = sx*dfx*pix/(s*r83);
	*e = sx*fx;
	return;
    }
    
    if(beeforder>=0)
    {
	//(*LdLn[beeforder])(t, &fx, &dl);
	LdLnACC(t, &fx, &dl, beeforder);
	dfx = dl*( 4.*s / (4.+s2) - 4.*s2*s/sq(4.+s2) );
	*dr = dx*fx - 4./3.*s2/(s*(*r))*sx*dfx;
	*dg = sx*dfx*pix/(s*r83);
	*e = sx*fx;
    }
    else
    {
	*dr = 0.;
	*dg = 0.;
	*e = 0.;
    }
    break;
    }
}


// evaluate local part of bee correlation and its derivatives de/drho and ( de/d|grad rho| ) / |grad rho|
#pragma acc routine seq
void beeflocalcorr_(double *r, double *g, double *e, double *dr, double *dg, int *addlda)
{
    double rs, ldac, ldadr, pbec, pbedr, pbed2rho;

    if(beeforder>=0)
    {
	*e = 0.;
	*dr = 0.;
	*dg = 0.;
	return;
    }
    
    switch(beeftype) {
    case 0: //BEEF-vdW xc
    rs = invpi075tothird / pow(*r,1./3.);
    corpbe(rs, 0.5/r2k * sqrt(*g*rs) / (*r),(beeforder>-3), 1, &ldac, &ldadr, &pbec, &pbedr, &pbed2rho);

    if(beeforder==-1)
    {
	if(!(*addlda))
	{
	    *e = beefpbecfrac*pbec*(*r);
	    *dr = beefpbecfrac*pbedr;
	}
	else
	{
	    *e = (beefpbecfrac*pbec+ldac)*(*r);
	    *dr = beefpbecfrac*pbedr+ldadr;
	}
	*dg = beefpbecfrac*pbed2rho / (*r);
	return;
    }
    
    if(beeforder==-2)
    {
	*e = pbec*(*r);
	*dr = pbedr;
	*dg = pbed2rho / (*r);
    }
    else if(beeforder==-3)
    {
	*e = ldac*(*r);
	*dr = ldadr;
	*dg = 0.;
    }
    else
    {
	*e = 0.;
	*dr = 0.;
	*dg = 0.;
    }
    
    break;
    }
}

// evaluate bee exchange energy only
void beefxpot_(double *r, double *g, double *e, int *addlda)
{
    double s2,t,r43;
    const int n=nmax;
    const int i1=1;
    const int i2=1;
    
    double L[nmax]={1.};

    switch(beeftype) {
    case 0: //BEEF-vdW xc
    r43 = pow(*r, 4./3.);

    s2 = *g*pix / (r43*r43);
    t = 2.*s2/(4.+s2)-1.;

    if(beeforder==-1)
    {
	calcleg(t);
	
	if(!(*addlda))
	    *e = (ddot_(&n, mi, &i1, L, &i2) - 1.) * r2e * r43;
	else
	    *e = ddot_(&n, mi, &i1, L, &i2) * r2e * r43;
	return;
    }
    
    if(beeforder>=0)
	*e = (*Ln[beeforder])(t) * r2e * r43;
    else
	*e = 0.;

    break;
    }
}

// evaluate local part of bee correlation - energy only
void beeflocalcorrpot_(double *r, double *g, double *e, int *addlda)
{
    double rs, ldac, ldadr, pbec, pbedr, pbed2rho;
    
    if(beeforder>=0)
    {
	*e = 0.;
	return;
    }
    
    switch(beeftype) {
    case 0: //BEEF-vdW xc
    rs = invpi075tothird / pow(*r,1./3.);
    corpbe(rs, 0.5/r2k * sqrt(*g*rs) / (*r),
	(beeforder>-3), 0, &ldac, &ldadr, &pbec, &pbedr, &pbed2rho);

    if(beeforder==-1)
    {
	if(!(*addlda))
	    *e = beefpbecfrac*pbec*(*r);
	else
	    *e = (beefpbecfrac*pbec+ldac)*(*r);
	return;
    }
    
    if(beeforder==-2)
	*e = pbec*(*r);
    else if(beeforder==-3)
	*e = ldac*(*r);
    else
	*e = 0.;

    break;
    }
}



// evaluate local part of bee correlation for spin polarized system
#pragma acc routine seq
void beeflocalcorrspin_(double *r, double *z, double *g, double *e,
    double *drup, double *drdown, double *dg, int *addlda) {
    double rs, ldac, ldadrup, ldadrdown, pbec, pbedrup, pbedrdown, pbed2rho;
    
    if(beeforder>=0)
    {
	*e = 0.;
	*drup = 0.;
	*drdown = 0.;
	*dg = 0.;
	return;
    }
    
    switch(beeftype) {
    case 0: //BEEF-vdW xc
    rs = invpi075tothird / pow(*r,1./3.);
    corpbespin(rs, 0.5/r2k * sqrt(*g*rs) / (*r), *z,
	(beeforder>-3), 1, &ldac, &ldadrup, &ldadrdown, &pbec,
	&pbedrup, &pbedrdown, &pbed2rho);

    if(beeforder==-1)
    {
	if(!(*addlda))
	{
	    *e =  beefpbecfrac*pbec*(*r);
	    *drup = beefpbecfrac*pbedrup;
	    *drdown = beefpbecfrac*pbedrdown;
	}
	else
	{
	    *e =  (beefpbecfrac*pbec+ldac)*(*r);
	    *drup = beefpbecfrac*pbedrup+ldadrup;
	    *drdown = beefpbecfrac*pbedrdown+ldadrdown;
	}
	*dg = beefpbecfrac*pbed2rho / (*r);
	return;
    }
    
    if(beeforder==-2)
    {
	*e = pbec*(*r);
	*drup = pbedrup;
	*drdown = pbedrdown;
	*dg = pbed2rho / (*r);
    }
    else if(beeforder==-3)
    {
	*e = ldac*(*r);
	*drup = ldadrup;
	*drdown = ldadrdown;
	*dg = 0.;
    }
    else
    {
	*e = 0.;
	*drup = 0.;
	*drdown = 0.;
	*dg = 0.;
    }
    
    break;
    }
}

// evaluate local part of bee correlation for spin polarized system - energy only
void beeflocalcorrpotspin_(double *r, double *z, double *g, double *e, int *addlda)
{
    double rs, ldac, ldadrup, ldadrdown, pbec, pbedrup, pbedrdown, pbed2rho;
    
    if(beeforder>=0)
    {
	*e = 0.;
	return;
    }

    switch(beeftype) {
    case 0: //BEEF-vdW xc
    rs = invpi075tothird / pow(*r,1./3.);
    corpbespin(rs, 0.5/r2k * sqrt(*g*rs) / (*r), *z,
	(beeforder>-3), 0, &ldac, &ldadrup, &ldadrdown, &pbec,
	&pbedrup, &pbedrdown, &pbed2rho);

    if(beeforder==-1)
    {
	if(!(*addlda))
	    *e = beefpbecfrac*pbec*(*r);
	else
	    *e = (beefpbecfrac*pbec+ldac)*(*r);
	return;
    }
    
    if(beeforder==-2)
	*e = pbec*(*r);
    else if(beeforder==-3)
	*e = ldac*(*r);
    else
	*e = 0.;

    break;
    }
}



// mode >= 0: for perturbed parameters --- calc Legendre order mode only
// -1:        standard beefxc expansion coefficients
// -2:        PBE correlation only
// -3:        LDA correlation only
// else:      no correlation either
void beefsetmode_(int *mode)
{
    beeforder = *mode;
#pragma acc update device(beeforder)
}

// initialize pseudo random number generator
void beefrandinit_(unsigned int *seed)
{
    srandom(*seed);
}

// initialize pseudo random number generator with default seed
void beefrandinitdef_()
{
    srandom(defaultseed);
}

// calculate ensemble energies
void beefensemble_(double *beefxc, double *ensemble)
{
    double vec[nmax+2],randvec[nmax+1];
    const double alpha=1.;
    const double beta=0.;
    const int m=nmax+1;
    const int n=nmax+1;
    const int la=nmax+1;
    const int ix=1;
    const int iy=1;
    const int n2=nmax+2;

    int i=0;
    int j=0;
    for(i=0;i<nsamples;i++)
    {
	for(j=0;j<nmax+1;j++) randvec[j] = normrand();
	dgemv_("T", &m, &n, &alpha, beefmat, &la, randvec, &ix,
	    &beta, vec, &iy);
	vec[nmax+1] = -vec[nmax];
	ensemble[i] = ddot_(&n2, vec, &ix, beefxc, &iy);
    }
}


//set type of beef functional to be used
//returns true on success
//0: BEEF-vdW
int beef_set_type_(int *tbeef, int *ionode)
{
    beeftype = *tbeef;
#pragma acc update device(beeftype)
    
    if(*ionode)
    {
	puts("\n" output_spacing output_marker);
	printf(output_spacing "Initializing " PACKAGE " V" VERSION " ");
	
	switch(beeftype) {
	case 0:
	    puts("with the BEEF-vdW functional.");
	    puts(output_spacing "Citation: Wellendorff et al., PRB 85, 235149 (2012).");
	    break;

	default:
	    return 0;
	}
	
	puts(output_spacing output_marker "\n");
	fflush(stdout);
    }
    
    return 1;
}
