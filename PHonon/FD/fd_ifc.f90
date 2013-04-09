!==========================================

!* File Name : ifc.f90

!* Creation Date : 25-12-2012

!* Last Modified : Tue Dec 25 13:01:29 2012

!* Created By : Marco Buongiorno Nardelli 

!==========================================

program ifc

implicit none

INTEGER, PARAMETER :: DP=KIND(1.0d0)
 
CHARACTER(4), ALLOCATABLE    :: atom_name(:)
CHARACTER(50)    :: fileout, file_out, file_in, file_force, file_force3, file_rmsd
CHARACTER(50)    :: cna,ci,cnb,cj,cnc,ck,cbx,cr,cs,ttemp,cnx

INTEGER :: na, nb, nc, i, j, k, n, ierr, idum, nrx, nrh, nat, nbx, nr, nrr, natx, nax, ncx, ns, nave, ntemp,inn
INTEGER :: nrx1, nrx2, nrx3, nrr1, nrr2, nrr3, nr1, nr2, nr3, nbb, nr1a,nr2a,nr3a,nr1b,nr2b,nr3b
REAL(KIND=DP) :: de, d_ave, d_ave_ave,sumd3, msdtot, phi_ijj, phi_ikk,phiddum

REAL(KIND=DP),    ALLOCATABLE     :: temp(:) 
REAL(KIND=DP),    ALLOCATABLE     :: rmsd(:,:) 
REAL(KIND=DP),    ALLOCATABLE     :: force0(:,:) 
REAL(KIND=DP),    ALLOCATABLE     :: force(:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid(:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid_symm(:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: force3_11(:,:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: force3_12(:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: force3(:,:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid3(:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phidD3(:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid_symm3(:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid_T(:,:,:,:)

INTEGER, ALLOCATABLE :: inat(:,:,:,:), inbt3(:,:,:), inct3(:,:,:)

REAL(KIND=DP), PARAMETER :: k_B=8.6173324d-5/13.60658
LOGICAL :: d3, verbose, offdiagonal

!
! input namelist
!
NAMELIST /input/ nat, de, nrx, nrx1, nrx2, nrx3, d3, file_force, file_force3, file_out, file_rmsd, verbose, offdiagonal

nat = 0
nrx1 = 1
nrx2 = 1
nrx2 = 1
de = 0.002 ! in Angstrom
verbose =.false.
offdiagonal=.false.
file_rmsd = ' '
file_force = ' '
file_force3 = ' '
file_out = ' '
d3 = .false.

READ(5,input,IOSTAT=ierr)
IF ( ierr /= 0 ) stop 'error in reading input'

IF (nat == 0) stop 'no atoms!'

natx=nat*nrx1*nrx2*nrx3

allocate (force0(3,natx))
allocate (force(2,3,3,natx,natx))
allocate (phid(3,3,natx,natx))
allocate (phid_symm(3,3,natx,natx))
allocate (force3_11(2,3,3,3,natx,natx,natx))
allocate (force3_12(3,3,3,natx,natx,natx))
allocate (force3(2,3,3,3,natx,natx,natx))
allocate (phid3(3,3,3,natx,natx,natx))
allocate (phidD3(3,3,3,natx,natx,natx))
allocate (phid_symm3(3,3,3,natx,natx,natx))
allocate (phid_T(3,3,natx,natx))
allocate (inat(nat,nrx1,nrx2,nrx3))
allocate (inbt3(nrx1,natx,nat))
allocate (inct3(nrx1,natx,nat))

de=de/0.529177

! read forces for null displacements (residual)

force0(:,:)=0.d0

i=0
na=0
write(cna,*) na
cna=adjustl(cna)
write(ci,*) i
ci=adjustl(ci)
file_in=TRIM(file_force)//'.'//TRIM(ci)//'.'//TRIM(cna)//'.'//TRIM(cna)
OPEN(2,FILE=TRIM(file_in),FORM='formatted')
do nb=1,natx
   read(2,*) force0(1,nb),force0(2,nb),force0(3,nb)
end do
CLOSE(2)

! read forces

force(:,:,:,:,:)=0.d0

inn=1

do nb=1,nat
   do j=1,3
      write(cnx,*) inn
      cnx=adjustl(cnx)
      write(cnb,*) nb
      cnb=adjustl(cnb)
      write(cj,*) j
      cj=adjustl(cj)
      file_in=TRIM(file_force)//'.'//TRIM(cnx)//'.'//TRIM(cj)//'.'//TRIM(cnb)
      OPEN(2,FILE=TRIM(file_in),FORM='formatted')
      do na=1,natx
         read(2,*) force(inn,1,j,na,nb),force(inn,2,j,na,nb),force(inn,3,j,na,nb)
         do i=1,3
            force(inn,i,j,na,nb)=force(inn,i,j,na,nb)-force0(i,na)
         end do
      end do
      CLOSE(2)
   end do
end do

inn=2

do nb=1,nat
   do j=1,3
      write(cnx,*) inn
      cnx=adjustl(cnx)
      write(cnb,*) nb
      cnb=adjustl(cnb)
      write(cj,*) j
      cj=adjustl(cj)
      file_in=TRIM(file_force)//'.'//TRIM(cnx)//'.'//TRIM(cj)//'.'//TRIM(cnb)
      OPEN(2,FILE=TRIM(file_in),FORM='formatted')
      do na=1,natx
         read(2,*) force(inn,1,j,na,nb),force(inn,2,j,na,nb),force(inn,3,j,na,nb)
         do i=1,3
            force(inn,i,j,na,nb)=force(inn,i,j,na,nb)-force0(i,na)
         end do
      end do
      CLOSE(2)
   end do
end do

! populate the force matrix with pbc 

nax=0
do nr1=1,nrx1
   do nr2=1,nrx2
      do nr3=1,nrx3
         do na=1,nat
            nax=nax+1
            inat(na,nr1,nr2,nr3)=nax  !the indexing has to match the order of the atoms in the supercell constructed in displacement3D.f90
         end do
      end do
   end do
end do

do nr1=1,nrx1
   do nr2=1,nrx2
      do nr3=1,nrx3
         do nb=1,nat
            nbb=inat(nb,nr1,nr2,nr3)
            do nrr1=1,nrx1
               nax=nrr1+nr1-1
               IF (nrr1+nr1-1 > nrx1) nax=nrr1+nr1-1-(nrx1+1)+1
               do nrr2=1,nrx2
                  nbx=nrr2+nr2-1
                  IF (nrr2+nr2-1 > nrx2) nbx=nrr2+nr2-1-(nrx2+1)+1
                  do nrr3=1,nrx3
                     ncx=nrr3+nr3-1
                     IF (nrr3+nr3-1 > nrx3) ncx=nrr3+nr3-1-(nrx3+1)+1
!		 print*, nrr1+nr1-1, nrr2+nr2-1, nrr3+nr3-1, nax, nbx, ncx
                     do na=1,nat
                        do j=1,3
                           do i=1,3
                              force(1,i,j,inat(na,nax,nbx,ncx),nbb)=force(1,i,j,inat(na,nrr1,nrr2,nrr3),nb)
                              force(2,i,j,inat(na,nax,nbx,ncx),nbb)=force(2,i,j,inat(na,nrr1,nrr2,nrr3),nb)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do
end do

!compute IFCs
do j=1,3
   do nb=1,natx
      do i=1,3
         do na=1,natx
            phid(i,j,na,nb)=-0.5*(force(1,i,j,na,nb)-force(2,i,j,na,nb))/de
         end do
      end do
   end do
end do

!symmetrization of the ifc's

do i=1,3
   do j=1,3
      do na=1,natx
         do nb=1,natx
           phid_symm(i,j,na,nb)=0.5*(phid(i,j,na,nb)+phid(j,i,nb,na))
         end do
      end do
   end do
end do

fileout=TRIM(file_out)//'.fc'
OPEN(3,FILE=TRIM(fileout))

WRITE (3,'(3i4)') nrx1,nrx2,nrx3
do i=1,3
   do j=1,3
      do na=1,nat
         do nb=1,nat
            WRITE (3,'(4i4)') i,j,na,nb
            do nr3=1,nrx3
               do nr2=1,nrx2
                  do nr1=1,nrx1
                     WRITE (3,'(3i4,2x,1pe18.11)') nr1,nr2,nr3,phid_symm(i,j,na,inat(nb,nr1,nr2,nr3))
                  end do
               end do
            end do
         end do
      end do
   end do
end do

CLOSE (3)

IF(verbose) THEN

! anharmonic IJJ term from one-displacement forces

phid3(:,:,:,:,:,:)=0.d0

do j=1,3
   do nb=1,natx
      do i=1,3
         do na=1,natx
            phid3(i,j,j,na,nb,nb)=-1.0*(force(1,i,j,na,nb)+force(2,i,j,na,nb))/de**2
         end do
      end do
   end do
end do

! symmetrization of the anharmonic ifc's

phid_symm3(:,:,:,:,:,:)=0.d0

do i=1,3
   do j=1,3
      do na=1,natx
         do nb=1,natx
            phid_symm3(i,j,j,na,nb,nb)=phid3(i,j,j,na,nb,nb)
            phid_symm3(j,i,j,nb,na,nb)=phid3(i,j,j,na,nb,nb)
            phid_symm3(j,j,i,nb,nb,na)=phid3(i,j,j,na,nb,nb)
         end do
      end do
   end do
end do

fileout=TRIM(file_out)//'.IJJ'//'.fc3'
OPEN(3,FILE=TRIM(fileout))

WRITE (3,'(3i4)') nrx1,nrx2,nrx3
do i=1,3 
   do j=1,3
      do k=1,3
         do na=1,nat
            do nb=1,nat
               do nc=1,nat
                  WRITE (3,'(6i4)') i,j,k,na,nb,nc
                  do nr3=1,nrx3
                     do nr2=1,nrx2
                        do nr1=1,nrx1
                        !   WRITE (3,'(3i4,2x,1pe18.11)') nr1,nr2,nr3,phid_symm3(i,j,k,na,inat(nb,nr1,nr2,nr3),inat(nc,nr1,nr2,nr3))
                           WRITE (3,'(3i4,2x,1pe18.11)') nr1,nr2,nr3,phid3(i,j,k,na,inat(nb,nr1,nr2,nr3),inat(nc,nr1,nr2,nr3))
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do
end do

CLOSE (3)

END IF

! anharmonic ifc's

IF (d3) THEN

force3(:,:,:,:,:,:,:)=0.d0
force3_11(:,:,:,:,:,:,:)=0.d0
force3_12(:,:,:,:,:,:)=0.d0

! read forces for the two displacements case

write(cr,*) 1
cr=adjustl(cr)
 
do j=1,3
   do k=1,3
      do nb=1,natx
         do nc=1,natx
            write(cnb,*) nb
            cnb=adjustl(cnb)
            write(cj,*) j
            cj=adjustl(cj)
            write(cnc,*) nc
            cnc=adjustl(cnc)
            write(ck,*) k
            ck=adjustl(ck)
            file_in=TRIM(file_force3)//'.'//TRIM(cr)//'.'//TRIM(cj)//'.'//TRIM(ck)//'.'//TRIM(cnb)//'.'//TRIM(cnc)
            OPEN(2,FILE=TRIM(file_in),FORM='formatted',IOSTAT=ierr)
            if(ierr /= 0) stop 'error in reading'
            do na=1,natx
               read(2,*) force3(1,1,j,k,na,nb,nc),force3(1,2,j,k,na,nb,nc),force3(1,3,j,k,na,nb,nc)
               do i=1,3
                  force3(1,i,j,k,na,nb,nc)=force3(1,i,j,k,na,nb,nc)-force0(i,na)
               end do
            end do
            CLOSE(2)
         end do
      end do
   end do
end do

write(cr,*) 2
cr=adjustl(cr)

do j=1,3
   do k=1,3
      do nb=1,natx
         do nc=1,natx
            write(cnb,*) nb
            cnb=adjustl(cnb)
            write(cj,*) j
            cj=adjustl(cj)
            write(cnc,*) nc
            cnc=adjustl(cnc)
            write(ck,*) k
            ck=adjustl(ck)
            file_in=TRIM(file_force3)//'.'//TRIM(cr)//'.'//TRIM(cj)//'.'//TRIM(ck)//'.'//TRIM(cnb)//'.'//TRIM(cnc)
            OPEN(2,FILE=TRIM(file_in),FORM='formatted',IOSTAT=ierr)
            if(ierr /= 0) stop 'error in reading'
            do na=1,natx
               read(2,*) force3(2,1,j,k,na,nb,nc),force3(2,2,j,k,na,nb,nc),force3(2,3,j,k,na,nb,nc)
               do i=1,3
                  force3(2,i,j,k,na,nb,nc)=force3(2,i,j,k,na,nb,nc)-force0(i,na)
               end do
            end do
            CLOSE(2)
         end do
      end do
   end do
end do

IF(offdiagonal) THEN

IF(offdiagonal) stop 'NOT IMPLEMENTED!'

END IF

! IFCs

do i=1,3
   do j=1,3
      do k=1,3
         do na=1,natx
            do nb=1,natx
               do nc=1,natx
                  phidD3(i,j,k,na,nb,nc)=(-0.5d0)*(force3(1,i,j,k,na,nb,nc)+force3(2,i,j,k,na,nb,nc))/de**2 !-0.5*(phid3(i,j,j,na,nb,nb)+phid3(i,k,k,na,nc,nc))
               end do
            end do
         end do
      end do
   end do
end do

! symmetrization of the anharmonic ifc's

do i=1,3
   do j=1,3
      do k=1,3
         do na=1,natx
            do nb=1,natx
               do nc=1,natx
                  phiddum=(1.d0/6.d0)* &
                           (phidD3(i,j,k,na,nb,nc)+phidD3(j,k,i,nb,nc,na)+phidD3(k,i,j,nc,na,nb)+ &
                            phidD3(i,k,j,na,nc,nb)+phidD3(j,i,k,nb,na,nc)+phidD3(k,j,i,nc,nb,na))
                  phid_symm3(i,j,k,na,nb,nc)=phiddum
                  phid_symm3(j,k,i,nb,nc,na)=phiddum
                  phid_symm3(k,i,j,nc,na,nb)=phiddum
                  phid_symm3(i,k,j,na,nc,nb)=phiddum
                  phid_symm3(j,i,k,nb,na,nc)=phiddum
                  phid_symm3(k,j,i,nc,nb,na)=phiddum
               end do
            end do
         end do
      end do
   end do
end do

IF (verbose) THEN

fileout=TRIM(file_out)//'.fc3'
OPEN(3,FILE=TRIM(fileout))

do i=1,3
   do j=1,3
      do k=1,3
         do na=1,natx
            do nb=1,natx 
               do nc=1,natx 
!                  WRITE (3,'(6i4)') i,j,k,na,nb,nc
                  WRITE (3,'(6i4,2x,1pe18.11)') na,i,nb,j,nc,k,phid_symm3(i,j,k,na,nb,nc)
               end do
            end do
         end do
      end do
   end do
end do

CLOSE (3)

END IF

fileout=TRIM(file_out)//'.d3.fc'
OPEN(3,FILE=TRIM(fileout))

WRITE (3,'(3i4)') nrx1,nrx2,nrx3
do i=1,3
   do j=1,3
      do k=1,3
         do na=1,nat
            do nb=1,nat
               do nc=1,nat
                  WRITE (3,'(6i4)') na,i,nb,j,nc,k
                  do nr3b=1,nrx3
                     do nr2b=1,nrx2
                        do nr1b=1,nrx1
                           do nr3a=1,nrx3
                              do nr2a=1,nrx2
                                 do nr1a=1,nrx1
                                    WRITE (3,'(6i4,2x,1pe18.11)') nr1a,nr2a,nr3a,nr1b,nr2b,nr3b, &
                                    phid_symm3(i,j,k,na,inat(nb,nr1a,nr2a,nr3a),inat(nc,nr1b,nr2b,nr3b))
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do
end do

CLOSE (3)

!WRITE (*,'(6i4,2x,1f18.11)') 1,1,1,2,1, 3,1000*phid_symm3(1,2,3,1,1,1)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,1,1,1,2, 1,1000*phid_symm3(1,1,1,1,1,2)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,1,1,2,2, 1,1000*phid_symm3(1,2,1,1,1,2)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,2,1,2,2, 1,1000*phid_symm3(2,2,1,1,1,2)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,2,1,3,2, 1,1000*phid_symm3(2,3,1,1,1,2)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,1,1,2,6, 1,1000*phid_symm3(1,2,1,1,1,6)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,1,1,1,11,1,1000*phid_symm3(1,1,1,1,1,5)/2.0
!!WRITE (*,'(6i4,2x,1f18.11)') 1,2,1,2,11,1,1000*phid_symm3(2,2,1,1,1,5)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,1,1,1,11,2,1000*phid_symm3(1,1,2,1,1,5)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,1,1,2,11,2,1000*phid_symm3(1,2,2,1,1,5)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,2,1,2,11,2,1000*phid_symm3(2,2,2,1,1,5)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,2,1,3,11,2,1000*phid_symm3(2,3,2,1,1,5)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,3,1,3,11,2,1000*phid_symm3(3,3,2,1,1,5)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,1,2,1,13,2,1000*phid_symm3(1,1,2,1,2,9)/2.0
!!WRITE (*,'(6i4,2x,1f18.11)') 1,1,2,2,13,2,1000*phid_symm3(1,2,2,1,2,9)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,1,2,3,13,2,1000*phid_symm3(1,3,2,1,2,9)/2.0
!WRITE (*,'(6i4,2x,1f18.11)') 1,2,2,3,13,3,1000*phid_symm3(2,3,3,1,2,9)/2.0

END IF

deallocate (force0)
deallocate (force)
deallocate (phid)
deallocate (phid_symm)
deallocate (force3_11)
deallocate (force3_12)
deallocate (force3)
deallocate (phid3)
deallocate (phidD3)
deallocate (phid_symm3)
deallocate (phid_T)
deallocate (inat)
deallocate (inbt3)
deallocate (inct3)

end program ifc
