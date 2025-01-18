!#define DEBUG

PROGRAM merge
 ! Utility program to merge different U matrices in a single U
 ! Usage: merge_Umat.x <what> <path/to/what1> <path/to/what2> ... <path/to/whatN>
 ! <what>:  U, centres
 !          U       -> merge matrices or centres (seedname_u.mat)
 !          centres -> merge centres (seedname_centres.xyz)
 !
 IMPLICIT NONE
 character(256) :: arg, what
 character(256), allocatable :: filename(:)
 integer :: nargs, iarg
 logical :: l_merge_U = .FALSE.
 logical :: l_merge_centres = .FALSE.
 !
 nargs = command_argument_count()
 !write(*,*) nargs
 IF (nargs == 0 ) STOP
 allocate (filename(nargs) )
 WRITE(*,*) "Merging ",nargs-1, " blocks" !, what
 !
 iarg = 1
 CALL get_command_argument( iarg, arg )
 what=trim(arg)
 WRITE(*,*) "what = ",  trim(what)
 !
 SELECT CASE( trim( what ) )
 CASE( 'U' )
    l_merge_U  = .true.
 CASE ( 'centres')
    l_merge_centres  = .true.
 CASE DEFAULT
    WRITE (*,*) "what = ", trim( what ), " not implemented"
    WRITE (*,*) "Available options: U, centres"
    STOP
 END SELECT
 !
 IF ( l_merge_U ) CALL merge_U (nargs, filename)
 IF (l_merge_centres ) CALL merge_centres (nargs, filename)
 !
 WRITE (*,'(/, " GAME OVER", /)')
 !
END PROGRAM

SUBROUTINE merge_centres (nargs, filename)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nargs
 CHARACTER (256) :: filename(nargs)
 CHARACTER (256) :: dum,arg
 REAL*8, ALLOCATABLE :: centres (:,:)
 REAL*8, ALLOCATABLE :: atoms_pos (:,:)
 INTEGER :: iline, nlines, iarg, ncentres(nargs), icentres, i, ncentres_tot
 INTEGER :: nat
 CHARACTER (2), ALLOCATABLE :: atoms_name(:)
 !
 ncentres(:) = 0
 iarg = 1
 DO WHILE (iarg < nargs)
   iarg = iarg + 1
   CALL get_command_argument( iarg, arg )
   filename(iarg)=trim(arg)
   OPEN (27, file=filename(iarg))
   READ(27,*) nlines
   READ(27,*) dum
   DO iline = 1 , nlines
    READ (27,*) dum
    IF (TRIM(dum) =="X") ncentres(iarg)=ncentres(iarg)+1
   ENDDO
   CLOSE (27)
   WRITE(*,*) "block # ", iarg-1, "file = ", trim(filename(iarg)), " ncentres = ", ncentres(iarg)
   ncentres(1) = ncentres(1) + ncentres(iarg)
 ENDDO
 !WRITE(*,*) "number of centres = ", ncentres(2:nargs)
 WRITE(*,*) "Final number of centres = ", ncentres(1)
 !
 ALLOCATE (centres(3, ncentres(1)) )
 nat = nlines - ncentres(nargs)
 WRITE(*,*) "number of atoms = ", nat
 ALLOCATE (atoms_name (nat))
 ALLOCATE (atoms_pos (3, nat) )

 !
 icentres = 0
 centres = 0.D0
 DO iarg = 2, nargs
   WRITE(*,*) "Reading ", TRIM(filename(iarg))
   OPEN (27, file=filename(iarg))
   READ(27,*) nlines
   READ(27,*) dum
   DO iline = 1 , ncentres(iarg)
     icentres = icentres + 1
     IF (icentres .gt. ncentres(1)) THEN
        WRITE(*,*) "Something wrong, too many centres", icentres, ncentres(1)
        STOP
     ENDIF
     READ(27,*) dum, (centres (i,icentres), i =1,3)
   ENDDO
   IF (iarg == 2 ) THEN
     DO iline = 1, nat
       READ(27,*) atoms_name(iline), (atoms_pos (i,iline), i =1,3)
     ENDDO
   ENDIF
   CLOSE (27)
 ENDDO
 !WRITE(*,*) icentres, ncentres(1)
 !WRITE(*,*) (atoms_name(iline), (atoms_pos (i,iline), i =1,3) ,iline =1, nat)
 !
 WRITE(*,*) "Writing ", TRIM('wann_centres.xyz')
 open (40, file='wann_centres.xyz', form='formatted')
 write(40, '(i6)') ncentres(1)+nat
 write (40, *) "Wannier centres, written by merge.x"
 DO icentres = 1, ncentres(1)
   write (40, '("X",6x,3(f14.8,3x))') (centres(i, icentres), i=1, 3)
 ENDDO
 DO icentres = 1, nat
   WRITE(40, '(a2,5x,3(f14.8,3x))') atoms_name(icentres), (atoms_pos (i,icentres), i=1,3)
 ENDDO
 !
 DEALLOCATE (centres)
END SUBROUTINE

SUBROUTINE merge_U (nargs, filename)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nargs
 CHARACTER(256) :: filename(nargs)
 INTEGER :: iarg, nb1_, ni, nf
 character(256) :: dum,arg
 integer :: nk1, nk2, nb1, nb2
 integer :: nb
 complex*16, allocatable :: U(:,:,:)
 real*8, allocatable :: k(:,:)
 integer :: ik, ib, jb, ii, iib, jjb
 !
 nb=0
 iarg = 1
 DO WHILE ( iarg < nargs )
   iarg = iarg +1
   CALL get_command_argument( iarg, arg )
   !write (*,*) trim(arg), seedname, 'u_mat'
   filename(iarg)=trim(arg) !//'/'//trim(seedname)//trim('_u.mat')
   open (27, file=filename(iarg))
   read(27,*) dum
   read(27,*) nk1, nb1, nb1
   close(27)
   WRITE(*,*) "block # ", iarg-1, "file = ", trim(filename(iarg)), " nbnd = ", nb1
   nb = nb + nb1
 ENDDO
 !
 write(*,*) "Final # of bands   = ", nb
 write(*,*) "Final # of kpoints = ", nk1
 ALLOCATE( U(nb,nb,nk1))
 ALLOCATE (k(3,nk1))
 U=cmplx(0.D0,0.D0)
 !
 nb1_=0
 do iarg = 2, nargs
   open (27, file=filename(iarg))
   WRITE(*,*) "Reading ", TRIM(filename(iarg))
   read(27,*) dum
   read(27,*) nk1, nb1, nb1
   !write(*, '(3I12)') nk1, nb1, nb1
   ni=1+nb1_
   nf=nb1+nb1_
   write(*,*) "istart = ", ni, " iend = ", nf !, " increment =", nb1_
   do ik = 1, nk1
      read (27, *)
      !write(*,*)
      read  (27, '(f15.10,sp,f15.10,sp,f15.10)') k(:, ik)
      !write (* , '(f15.10,sp,f15.10,sp,f15.10)') k(:, ik)
      !write(*,'(f15.10,sp,f15.10,sp,f15.10)')  k(:,ik)
      read  (27, '(f15.10,sp,f15.10)') ((U(ib, jb, ik), ib=ni, nf), jb=ni, nf)
      !write (* , '(f15.10,sp,f15.10)') ((U(ib, jb, ik), ib=ni, nf), jb=ni, nf)
      !read (27, *) ((U(ib, jb, ik), ib=ni, nf), jb=ni, nf)
   enddo
   nb1_=nb1_+nb1
   close(27)
   !
 end do
 !
 write(*,*) "Final # of bands   = ", nb
 write(*,*) "Final # of kpoints = ", nk1
 !WRITE(*,'(2(f15.10,sp,f15.10))') U(1,1,1), U(nb,nb,nk1)
 WRITE(*,*) "Writing ", TRIM('wann_u.mat')
 open (40, file='wann_u.mat', form='formatted')
 write (40, *) "Written by merge.x"
 write (40, '(3I12)') nk1, nb, nb
 do ik = 1, nk1
   write(40, *)
   write(40, '(f15.10,sp,f15.10,sp,f15.10)') k(:, ik)
   write(40, '(f15.10,sp,f15.10)') ((U(ib, jb, ik), ib=1, nb), jb=1, nb)
   !write(40, *) ((U(ib, jb, ik), ib=1, nb), jb=1, nb)
   !write(*, *) ((U(ib, jb, ik), ib=1, nb), jb=1, nb)
 enddo
 close(40)
 !
 !
END subroutine
