!
! Copyright (C) 2006 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! GNU License
!
! Eyvaz Isaev
! 
! Theoretical Physics Department,
! Moscow State Institute of Steel and Alloys
! (Technological University)
!
! Condensed Matter Theory Group,
! Uppsala University, Sweden
!
! Eyvaz.Isaev@fysik.uu.se
! eyvaz_isaev@yahoo.com
!
! Adopted from a program published in a preprint (early 90th) of Lebedev Physical Institute (Moscow)
! Early versions of this program was used in ab initio psedopotentials code of E.I.Isaev
!  

	subroutine integration(e_init,tdos,dos,lpartial)
	include 'parameters.h'
	dimension  e0(4),pnt1(4,4), gx(3)
C
        ef=e_init
	irec=14*nzone
	jrec=14*nzone*nzone

      open(unit=17,file='eigenv',access='direct',recl=irec,
     *form='unformatted')

      open(unit=21,file='partial_DOS',access='direct',recl=jrec,
     *form='unformatted')

!  NB!!! in phonon calculations NZONE=3*NATOMS

      N6=NZONE
      NS=NZONE
      DOS=0.D0
      TDOS=0.D0

	natoms=n6/3	

c
      DO 1 I=1,NS
      rog(I)=0.D0
1     roz(i)=0.D0

      NTETMX=NTET(1)
      DO 120 KTET0=1,NT0
         LKTET=NTETMX*(KTET0-1)
         LKTET1=LKTET+1
         DO 30 I=1,4
            IPNT=TTR(I,LKTET1)
            DO 20 K=1,3
               PNT1(K,I)=PNT(K,IPNT)
 20         CONTINUE
            PNT1(4,I)=1D0
 30      CONTINUE
	ot=abs(det4(pnt1))/6.0
         NKTT=NTET(KTET0)
         DO 110 KTET=1,NKTT
            LK=LKTET+KTET
         DO 70 NZ=1,N6
            DO 60 I=1,4
               IPNT=TTR(I,LK)
	read(17,rec=ipnt)(e(j),j=1,n6)
	read(21,rec=ipnt)((har(j,k),k=1,ns),j=1,n6)
	e0(i)=e(nz)
              DO 2 J=1,NS
	a0(i,j)=har(nz,j)
2        CONTINUE
60       CONTINUE

                  CALL Tetrahedra(e0,ns,1)

                  TDOS=TDOS+TS
                  DOS=DOS+DS
              DO 3 J=1,NS
         rog(J)=rog(J)+ATS(J)
         roz(J)=roz(J)+ADS(J)
 3       CONTINUE
 70            CONTINUE
 110     CONTINUE
 120  CONTINUE
!
! We norm phonon DOS  so that  the integrated DOS is equal 3N 
! No spin polarization in phonon calcultions, so we introduce a factor 0.5
! 
      TDOS=TDOS*OMG48*0.5
      DOS=DOS*OMG48*0.5
      DO 4 I=1,NS
      rog(I)=rog(I)*OMG48*0.5
4     roz(I)=roz(I)*OMG48*0.5
       s1=roz(1)+roz(2)+roz(3)

 99     FORMAT(' Freq,Tot_DOS,DOS==',7(1X,G14.7))
        WRITE (6,99) ef,TDOS,DOS,s1,roz(1),roz(2),roz(3)

	do natom=1,natoms
        iunit=30+natom
        open(unit=iunit,file='projected_DOS.'//atom(natom),
     *	access='sequential',form='formatted')
	
        if(ef.eq.0.) then
	write(iunit,'("#",7X,"E",10X," DOS",14X,"PDOS",13X,"g_x", 
     *	13X,"g_y",13X,"g_z")')
	write(iunit,'("# ",i4,f14.6,f6.3)') nstep,emax,delta_e 
	endif
		
	pdos=0.d0
	do imode=1,3
	pdos=pdos+roz(3*(natom-1)+imode)
	gx(imode)=roz(3*(natom-1)+imode)
	enddo
	
        WRITE (iunit, '(f12.4,5(3XG14.7))') ef,dos,pdos,gx(1),
     *	gx(2),gx(3)
	enddo 
	
	close(17)
	close(21)

      RETURN
      END 
