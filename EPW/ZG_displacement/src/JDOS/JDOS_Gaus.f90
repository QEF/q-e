 PROGRAM JDOS_gaus
 !               Marios Zacharias & Feliciano Giustino, April 2018
 ! Department of Materials, University of Oxford, Parks Road, Oxford OX1 3PH, United Kingdom
 !
 !----------------------------------------------------------------------------------------
 !
 ! JDOS_gaus.x:
 !  Reads the eigenvalues produced by a scf or nscf calcualtion on the "ZG configuration" 
 !  and calculates the joint density of states, as described in Phys. Rev. B 94, 075125, 2016. 
 !
 ! Inputs: 
 ! 
 ! prefix            : prefix of scf/nscf calculation 
 ! 
 ! wmin, wmax, steps : define the energy range (in eV) and the grid.
 !  
 ! smearing          : smearing of the Gaussian broadening (units of eV/sqrt(2))
 ! 
 ! vb, cb            : number of valence and conduction states 
 ! 
 ! nkpts             : number of kpts
 ! 
   IMPLICIT NONE ! "implicit none" statement forces the programmer to declare all variables,
   ! Input files
        CHARACTER (LEN=256)             :: prefix !  Length of a character  entity   
      
        DOUBLE PRECISION, ALLOCATABLE   ::  trans_ene(:), vb_ene(:),cb_ene(:)
   ! Output matrices
        DOUBLE PRECISION, ALLOCATABLE   ::  wgrid_out(:) , JDOS_out(:)
        
   ! Input variables
        INTEGER                         :: steps
        INTEGER                         :: vb,cb,nkpts
        DOUBLE PRECISION                :: smearing, wmin, wmax
   ! Generic counting integers
        INTEGER                         ::  i, j, v,ctr
   ! Constants
        DOUBLE PRECISION, PARAMETER     :: PI = 3.1415927, Ha_to_eV = 13.60569253*2
   ! Variables
        DOUBLE PRECISION                :: jump

   ! Outputs
        CHARACTER (LEN=256)             :: OUTPUT_1

    ! Input namelist
      NAMELIST/LIST_NAME/ steps,wmin, wmax, prefix, vb, cb, &
                          smearing, nkpts, OUTPUT_1 !Namelist, sequential READ statements, 


  READ(*,NML=LIST_NAME)

  ALLOCATE(vb_ene(nkpts*vb),cb_ene(nkpts*cb))

  call read_eigenvalues (nkpts,vb,cb,prefix,vb_ene,cb_ene)

  vb_ene(:)=vb_ene(:)*Ha_to_eV
  cb_ene(:)=cb_ene(:)*Ha_to_eV

  ! This loop is to count the entries of matrix "trans_ene" and allocate 
  ctr=0
    DO i=1, size(cb_ene)
     DO j=1, size(vb_ene)
         IF  (cb_ene(i) - vb_ene(j) .LE. wmax) THEN
         ctr=ctr+1
         END IF
     END DO
   END DO

  ALLOCATE(trans_ene(ctr))
  ctr=1
    DO i=1, size(cb_ene)
     DO j=1, size(vb_ene)
         IF  (cb_ene(i) - vb_ene(j) .LE. wmax) THEN
         trans_ene(ctr) = (cb_ene(i) - vb_ene(j))
         ctr=ctr+1
         END IF
     END DO
   END DO


  ALLOCATE(wgrid_out(steps))
  !
  jump = (wmax - wmin)/dble(steps)
  !
  DO i = 1, steps
    wgrid_out(i) = wmin + (i-1)*jump
  END DO
  !
  ALLOCATE(JDOS_out(steps))
  !
  JDOS_out = 0.d0
     DO i=1, steps
      DO j = 1, size(trans_ene)
        JDOS_out(i) = JDOS_out(i) + EXP(-(trans_ene(j) - wgrid_out(i))**2/2/smearing**2)/smearing/SQRT(2*PI)
     END DO
    END DO
!
! Print result but divide with the total number of entries 
  OPEN(50, FILE=OUTPUT_1)
  WRITE(50,'(A,F4.2,A,F4.2,A,i14)') '# Transitions within ', wmin,'-',wmax,' eV: ',size(trans_ene)
  DO i = 1, steps
    WRITE(50,'(2F20.8)') wgrid_out(i), JDOS_out(i)
  END DO

  CLOSE(50)

  DEALLOCATE( JDOS_out, wgrid_out, trans_ene, vb_ene, cb_ene)
 

 END PROGRAM JDOS_gaus 

 SUBROUTINE read_eigenvalues (nkpts,vb,cb,prefix,vb_ene,cb_ene) 

  implicit none
  integer,           intent(in)  :: nkpts, vb,cb
  CHARACTER(LEN=256),intent(in)  :: prefix
  DOUBLE PRECISION,  intent(out) :: vb_ene(vb*nkpts), cb_ene(cb*nkpts)
  integer                        :: i, j, ctr_vb,ctr_cb
  CHARACTER(LEN=256)             :: filename, foldername

  ctr_vb = 1 
  ctr_cb = 1

 filename = TRIM( prefix ) 
 OPEN(44,FILE=filename)
  DO i=1, nkpts
     DO j=1,vb
        READ(44,*) vb_ene(ctr_vb)
        ctr_vb = ctr_vb + 1
     END DO
     DO j=1,cb
        READ(44,*) cb_ene(ctr_cb)
        ctr_cb = ctr_cb + 1
     END DO
  END DO
 CLOSE(44)

 END SUBROUTINE
