! GNU License
! Eyvaz Isaev
! evyaz_isaev@yahoo.com
! Eyvaz.Isaev@fysik.uu.se
!
        implicit real*8(a-h,o-z)
        integer ttr0,ttr
	character*4 atom
        common /xyz3/x(3),y(3),z(3)
     *	  /edh/e(5000)  
     *    /int1/ats(500),ads(500),a0(4,500)
     *    /int2/rog(500),roz(500)
     *    /hrog/har(100,500)
     *    /mesh1/npnt,ntet(25),nt0
     *    /mesh2/pnt(3,10000),ttr(4,81000),omg48,pnt0(3,25),ttr0(4,25)
     *    /ttrinf/ndiv,npnt0,ntet0
     *    /int/ot,ef,ts,ds
     *    /zb/zbk(3), ki
     *    /vcell/vol
     *    /zones/nzone
     *    /types/atom(25)
     *    /fordos/emax,delta_e,nstep
