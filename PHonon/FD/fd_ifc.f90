!==========================================

!* File Name : fd_ifc.f90

!* Creation Date : 25-12-2012

!* Last Modified : Wed Jan 29 09:50:17 2014

!* Created By : Marco Buongiorno Nardelli 

!==========================================

program fd_ifc

  use constants,  ONLY : pi, bohr_radius_angs, amu_au, amu_ry
  use io_files,   ONLY : prefix, tmp_dir, psfile, pseudo_dir
  use io_global,  ONLY : stdout, ionode, ionode_id
  USE mp_global,  ONLY : mp_startup,mp_global_end
  USE environment,ONLY : environment_start,environment_end
  USE mp,         ONLY : mp_bcast
  USE cell_base,  ONLY : tpiba2, alat,omega, at, bg, ibrav, celldm
  USE ions_base,  ONLY : amass, nat, nat, atm, zv, tau, ntyp => nsp, ityp
  USE kinds,      ONLY : dp
  USE gvecw,      ONLY : ecutwfc
  USE matrix_inversion

  USE symm_base
  USE symme
  USE rap_point_group, ONLY : code_group, nclass, nelem, elem, &
       which_irr, char_mat, name_rap, name_class, gname, ir_ram
  USE rap_point_group_so, ONLY : nrap, nelem_so, elem_so, has_e, &
       which_irr_so, char_mat_so, name_rap_so, name_class_so, d_spin, &
       name_class_so1
  USE rap_point_group_is, ONLY : nsym_is, sr_is, ftau_is, d_spin_is, &
       gname_is, sname_is, code_group_is
  USE fft_base, ONLY : dfftp

implicit none
 
character(len=9) :: code = 'FD_IFC'
CHARACTER(4), ALLOCATABLE    :: atom_name(:)
CHARACTER(50)    :: fileout, file_out, file_in, file_force, file_force3, file_rmsd
CHARACTER(50)    :: cna,ci,cnb,cj,cnc,ck,cbx,cr,cs,ttemp,cnx
CHARACTER(256)   :: outdir

INTEGER :: na, nb, nc, i, j, k, n, ierr, idum, nrx, nrh, nbx, nr, nrr, natx, nax, ncx, ns, nave, ntemp,inn,innx
INTEGER :: nrx1, nrx2, nrx3, nrr1, nrr2, nrr3, nr1, nr2, nr3, nbb, nr1a,nr2a,nr3a,nr1b,nr2b,nr3b, nat_0
REAL(KIND=DP) :: de, d_ave, d_ave_ave,sumd3, msdtot, phi_ijj, phi_ikk,phiddum

REAL(KIND=DP),    ALLOCATABLE     :: temp(:) 
REAL(KIND=DP),    ALLOCATABLE     :: rmsd(:,:) 
REAL(KIND=DP),    ALLOCATABLE     :: force0(:,:) 
REAL(KIND=DP),    ALLOCATABLE     :: force(:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: fforce(:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid(:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid_symm(:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: force3_11(:,:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: force3_12(:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: force3(:,:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid3(:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phidD3(:,:,:,:,:,:)
REAL(KIND=DP),    ALLOCATABLE     :: phid_symm3(:,:,:,:,:,:)

INTEGER, ALLOCATABLE :: inat(:,:,:,:)
INTEGER, ALLOCATABLE :: ieq(:,:),seq(:,:),neq(:), atdp(:), sydp(:)
LOGICAL :: atom_in_list
real(kind=dp)     :: r1(3),r2(3),r3(3),rr(3,3),bg_0(3,3),at_0(3,3)
REAL(KIND=DP),    ALLOCATABLE     :: taut(:,:)
INTEGER :: ipol, apol, natdp, ios


INTEGER :: nclass_ref   ! The number of classes of the point group
INTEGER :: isym
INTEGER, ALLOCATABLE  ::  sxy(:), sxz(:), syz(:)
REAL (dp) :: ft1, ft2, ft3

  REAL (dp) :: d(3,3),rd(3,3),dhex(3,3), dcub(3,3), dhexm1(3,3), tmp(3)
  REAL (dp) :: accep=1.0d-5
  LOGICAL, ALLOCATABLE :: move_sl(:,:)
  LOGICAL, EXTERNAL :: eqdisp

  REAL(KIND=DP),    ALLOCATABLE     :: atomx(:,:), ttomx(:,:)
  integer,     ALLOCATABLE    :: itypx(:)

REAL(KIND=DP), PARAMETER :: k_B=8.6173324d-5/13.60658
LOGICAL :: ld3, verbose, offdiagonal, nodispsym, lforce0,hex, noatsym
  CHARACTER(LEN=256), EXTERNAL :: trimcheck

  real(DP), PARAMETER :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
                        msin3 =-0.866025403784438597d0, mcos3 = -0.5d0

  data dcub/ 1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
  data dhex/ 1.d0, 0.d0, 0.d0, mcos3, sin3, 0.d0, 0.d0, 0.d0,  1.d0/


!
! input namelist
!
NAMELIST /input/ prefix, de, nrx, nrx1, nrx2, nrx3, ld3, file_force, file_force3, file_out, file_rmsd, verbose, &
                 offdiagonal, innx, nodispsym, lforce0,hex,noatsym, outdir

  CALL mp_startup ( )
  CALL environment_start ( code )

prefix=' '
outdir='./'
nrx1 = 1
nrx2 = 1
nrx2 = 1
innx=2
de = 0.002 ! in Angstrom
verbose =.false.
offdiagonal=.false.
file_rmsd = ' '
file_force = ' '
file_force3 = ' '
file_out = ' '
nodispsym=.true.
noatsym=.false.
ld3 = .false.
lforce0=.true.
hex=.false.

CALL input_from_file ( )
READ(5,input,IOSTAT=ios)
IF (ios /= 0) CALL errore ('FD_IFC', 'reading input namelist', ABS(ios) )
tmp_dir = trimcheck( outdir )

!reading the xml file
call read_xml_file

    if (verbose) then
    write(6,*) '**************************************************'
    write(6,*) '* Info from the preceding pw.x run:              *'
    write(6,*) '**************************************************'
    write(6,*) ''
    write(6,*) '    prefix=  ',trim(prefix)
    write(6,*) '    outdir=  ',trim(tmp_dir)
    write(6,*) '    ectuwfc= ',ecutwfc, 'Ry'
    ALLOCATE(atom_name(nat))
    ALLOCATE(taut(3,nat))

    WRITE( stdout, '(i4,2f12.6,2i4)') ibrav, alat, omega, nat, ntyp
    !
    WRITE( stdout, '(/2(3X,3(2X,"celldm(",I1,")=",F11.6),/))' ) &
      ( i, celldm(i), i = 1, 6 )
    !
    !lattice vectors in Angs
    rr = at*alat*bohr_radius_angs
    r1(:) = rr(:,1)
    r2(:) = rr(:,2)
    r3(:) = rr(:,3)
    WRITE( stdout, '(5X, &
      &     "Lattice vectors: (cart. coord. in Angs)",/, &
      &       3(15x,"a(",i1,") = (",3f11.6," )  ",/ ) )')  (apol,  &
      (rr (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
    !atoms
    do na=1,nat
      atom_name(na)=atm(ityp(na))
    enddo
    do na=1,nat
      taut(:,na) = tau(:,na)*alat*bohr_radius_angs
    enddo
    WRITE( stdout, '(/,5x,"Atomic coordiantes")')
    WRITE( stdout, '(5x,"site n.     atom                  positions (Angs)")')
    WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
      (na, atom_name(na), na, (taut(ipol,na), ipol=1,3), na=1,nat)
    
   deallocate(atom_name)
   deallocate(taut)
   end if
IF (nat == 0) stop 'no atoms!'

! find equivalent atoms and check that forces have been calculated

    write(6,*) ''
    write(6,*) '**************************************************'
    write(6,*) '* Info on equivalent atoms and displacements:    *'
    write(6,*) '**************************************************'
    write(6,*) ''

    ALLOCATE (ieq(nat,nsym))
    allocate (neq(nat))
    allocate (seq(nat,nsym))
    allocate (atdp(nat))

    ieq=0
    seq=0
    neq=0
    atdp=0

    natdp = 1
    atdp(1)=1
    
    na=1

    ! natdp             # of non equivalent atoms
    ! atdp(natdp)       index of non equivalent atoms
    ! neq(na)           # of atoms equivalent to na (including na itself (identity always isym=1))
    ! ieq(na,1:neq(na)) list of equivalent atoms (including identity) => ieq(na,2:neq(na)) are the
    !                   equivalent atoms (ieq(na,1) is displaced, the others are not)
    ! seq(na,1:neq(na)) list of symmetry operations for equivalent atoms (including identity, seq(na,1)) => 
    !                   seq(na,2:neq(na)) are the operations that map equivalent atoms in the displaced one

    call equiv_atoms (na,nsym,irt(:,na),ieq(na,:),seq(na,:),neq(na))
    outer: do na=2,nat
       do j=1,na-1
          do k=1,neq(j)
             if(na == ieq(j,k)) then
               cycle outer
             end if
          end do
       end do
       call equiv_atoms (na,nsym,irt(:,na),ieq(na,:),seq(na,:),neq(na))
       natdp=natdp+1
       atdp(natdp)=na
    end do outer

    if (noatsym) then
       natdp =nat
       do j=1,nat
          atdp(j)=j
       end do
    end if

    write(*,fmt='(a,i0)') 'Number of independent atoms: ',natdp
    write(*,advance='no',fmt='(a)') 'Read forces for displaced atom(s):'
    write(*,*) atdp(1:natdp)

    natx=nat*nrx1*nrx2*nrx3

! define equivalent atoms in the supercell
! the indexing has to match the order of the atoms in the supercell constructed in fd.f90

allocate (inat(nat,nrx1,nrx2,nrx3))

nax=0
do nr1=1,nrx1
   do nr2=1,nrx2
      do nr3=1,nrx3
         do na=1,nat
            nax=nax+1
            inat(na,nr1,nr2,nr3)=nax  
         end do
      end do
   end do
end do

allocate (force0(3,natx))

de=de/0.529177 ! de is read in A and used in bohr.

! read forces for null displacements (residual)

force0(:,:)=0.d0

if (lforce0) then
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
end if

   nat_0=nat

   allocate (force(innx,3,3,natx,nat_0))
   allocate(sxy(nat_0))
   allocate(sxz(nat_0))
   allocate(syz(nat_0))
   allocate(move_sl(3,nat_0))

   move_sl = .true.
   sxy = 0
   sxz = 0
   syz = 0

if(verbose) then
    write(6,*) ''
    write(6,*) '**************************************************'
    write(6,*) '* Info on supercell geometry:                    *'
    write(6,*) '**************************************************'
    write(6,*) ''
end if

   ALLOCATE(itypx(natx))
   ALLOCATE(atomx(3,natx))
   ALLOCATE(ttomx(3,natx))

   nax=0
   do nr1=1,nrx1
     do nr2=1,nrx2
       do nr3=1,nrx3
         do na=1,nat_0
           nax=nax+1
           do i=1,3
            atomx(i,nax)=tau(i,na)+(nr1-1)*at(i,1)+(nr2-1)*at(i,2)+(nr3-1)*at(i,3)
           enddo
           itypx(nax)=ityp(na)
         enddo
       enddo
     enddo
   enddo

if (verbose) then
    WRITE( stdout, '(/,5x,"Atomic coordiantes")')
    WRITE( stdout, '(5x,"site n.     atom                  positions (alat)")')
    WRITE( stdout, '(6x,i4,8x,a6," tau(",i4,") = (",3f12.7,"  )")') &
      (na, atm(itypx(na)), na, (atomx(ipol,na), ipol=1,3), na=1,natx)
end if

   ! define the supercell vectors

   at(:,1)=nrx1*at(:,1)
   at(:,2)=nrx2*at(:,2)
   at(:,3)=nrx3*at(:,3)
   CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
   
   if ( .not.nodispsym ) then

   IF ( ALLOCATED( irt ) ) DEALLOCATE( irt )
   ALLOCATE( irt(48,natx))
   
   do k=1,natdp
      nb=atdp(k)
 ! center on the nb atom to find the correct symmetry operations
      do na=1,natx
         ttomx(:,na) = atomx(:,na) - atomx(:,nb)
      end do

      call find_sym_ifc( natx, ttomx, itypx )

!if (verbose) then
!do na=1,natx
!write(*,'(24i3)') (irt(i,na),i=1,nsym)
!end do
!print*, '============='
!end if 

if (hex) then
   d=dhex
   call invmat (3, dhex, dhexm1)
else
   d=dcub
end if

    ! find if the cartesian displacements are independent also in the supercell (as it should be!) and identify one
    ! symmetry operation to rotate forces

    xy_sl: do isym=1,nsym
       rd(:,1) = sr(:,1,isym) * d(1,1) + &
                 sr(:,2,isym) * d(2,1) + &
                 sr(:,3,isym) * d(3,1)
       if (eqdisp(rd(:,1),d(:,2),accep) .and. &
          ft(1,isym) == 0.d0 .and. &
          ft(2,isym) == 0.d0 .and. &
          ft(3,isym) == 0.d0) then
          write(6,'("on atom",i4," x and y displacements are equivalent for symmetry operation ",i4)') nb,isym
          move_sl(2,nb) = .false.
          sxy(nb)=isym
          exit
       end if
    end do xy_sl

    xz_sl: do isym=1,nsym
       rd(:,1) = sr(:,1,isym) * d(1,1) + &
                 sr(:,2,isym) * d(2,1) + &
                 sr(:,3,isym) * d(3,1)
       if (eqdisp(rd(:,1),d(:,3),accep) .and. &
          ft(1,isym) == 0.d0 .and. &
          ft(2,isym) == 0.d0 .and. &
          ft(3,isym) == 0.d0) then
          write(6,'("on atom",i4," x and z displacements are equivalent for symmetry operation ",i4)') nb,isym
          move_sl(3,nb) = .false.
          sxz(nb)=isym
          exit
       end if
    end do xz_sl

    yz_sl: do isym=1,nsym
       rd(:,2) = sr(:,1,isym) * d(1,2) + &
                 sr(:,2,isym) * d(2,2) + &
                 sr(:,3,isym) * d(3,2)
       if (eqdisp(rd(:,2),d(:,3),accep) .and. &
          ft(1,isym) == 0.d0 .and. &
          ft(2,isym) == 0.d0 .and. &
          ft(3,isym) == 0.d0) then
          write(6,'("on atom",i4," y and z displacements are equivalent for symmetry operation ",i4)') nb,isym
          move_sl(3,nb) = .false.
          syz(nb)=isym
          exit
       end if
    end do yz_sl

  end do

end if

! read forces and complete missing polarizations

force(:,:,:,:,:)=0.d0

do inn=1,innx
   do k=1,natdp
      nb=atdp(k)
      do j=1,3
         if (move_sl(j,nb)) then 
            write(cnx,*) inn
            cnx=adjustl(cnx)
            write(cnb,*) nb
            cnb=adjustl(cnb)
            write(cj,*) j
            cj=adjustl(cj)
            file_in=TRIM(file_force)//'.'//TRIM(cnx)//'.'//TRIM(cj)//'.'//TRIM(cnb)
            print*, 'reading ',file_in
            OPEN(2,FILE=TRIM(file_in),FORM='formatted')
            do na=1,natx
               read(2,*) force(inn,1,j,na,nb),force(inn,2,j,na,nb),force(inn,3,j,na,nb)
               do i=1,3
                  force(inn,i,j,na,nb)=force(inn,i,j,na,nb)-force0(i,na)
               end do
            end do
            CLOSE(2)
         end if
      end do
   end do
end do

! rotate forces if needed

if (.not.nodispsym) then

do inn=1,innx
   do k=1,natdp
      nb=atdp(k)
 ! center on the nb atom to find the correct symmetry operations
      do na=1,natx
         ttomx(:,na) = atomx(:,na) - atomx(:,nb)
      end do

      call find_sym_ifc( natx, ttomx, itypx )

      if (.not.move_sl(2,nb) .and. sxy(nb) .ne. 0) then
         do na=1,natx
            if (.not.hex) then 
            force(inn,:,2,na,nb) = sr(:,1,sxy(nb)) * force(inn,1,1,irt(invs(sxy(nb)),na),nb) + &
                                   sr(:,2,sxy(nb)) * force(inn,2,1,irt(invs(sxy(nb)),na),nb) + &
                                   sr(:,3,sxy(nb)) * force(inn,3,1,irt(invs(sxy(nb)),na),nb)
            else 
            tmp(:) = sr(:,1,sxy(nb)) * force(inn,1,1,irt(invs(sxy(nb)),na),nb) + &
                     sr(:,2,sxy(nb)) * force(inn,2,1,irt(invs(sxy(nb)),na),nb) + &
                     sr(:,3,sxy(nb)) * force(inn,3,1,irt(invs(sxy(nb)),na),nb)
            force(inn,:,2,na,nb) = dhexm1(:,1) * tmp(1) + dhexm1(:,2) * tmp(2) + dhexm1(:,3) * tmp(3)
            end if
         end do
      end if
      if (.not.move_sl(3,nb) .and. sxz(nb) .ne. 0) then
         do na=1,natx
            force(inn,:,3,na,nb) = sr(:,1,sxz(nb)) * force(inn,1,1,irt(invs(sxz(nb)),na),nb) + &
                                   sr(:,2,sxz(nb)) * force(inn,2,1,irt(invs(sxz(nb)),na),nb) + &
                                   sr(:,3,sxz(nb)) * force(inn,3,1,irt(invs(sxz(nb)),na),nb)
         end do
      end if
      if (.not.move_sl(3,nb) .and. syz(nb) .ne. 0) then
         do na=1,natx
            force(inn,:,3,na,nb) = sr(:,1,syz(nb)) * force(inn,1,2,irt(invs(syz(nb)),na),nb) + &
                                   sr(:,2,syz(nb)) * force(inn,2,2,irt(invs(syz(nb)),na),nb) + &
                                   sr(:,3,syz(nb)) * force(inn,3,2,irt(invs(syz(nb)),na),nb)
         end do
      end if
   end do
end do

end if

if (natdp .ne. nat_0) then

! compute forces on equivalent atoms

    ! natdp             # of non equivalent atoms
    ! atdp(natdp)       index on non equivalent atoms
    ! neq(na)           # of atoms equivalent to na (including na itself (identity always isym=1))
    ! ieq(na,1:neq(na)) list of equivalent atoms (including identity) => ieq(na,2:neq(na)) are the
    !                   equivalent atoms (ieq(na,1) is displaced, the others are not)
    ! seq(na,1:neq(na)) list of symmetry operations for equivalent atoms (including identity, seq(na,1)) => 
    !                   seq(na,2:neq(na)) are the operations that map equivalent atoms in the displaced one

IF ( ALLOCATED( irt ) ) DEALLOCATE( irt )
ALLOCATE( irt(48,natx))

   call find_sym_ifc( natx, atomx, itypx )
  
   do k=1,natdp
      j=atdp(k)  ! loop over the non equivalent atoms (displaced)
      do i=2,neq(j)  ! loop over the equivalent atoms for each of the non equivalent ones
         nb=ieq(j,i)
         isym=seq(j,i)
         do inn=1,innx
            do na=1,natx
               force(inn,:,:,na,nb)=matmul(matmul(sr(:,:,isym),force(inn,:,:,irt(invs(isym),na),j)),sr(:,:,invs(isym)))
            end do
         end do
      end do
   end do
end if


! write forces so far

if (verbose) then
do nb=1,nat_0
   do j=1,3
      write(6,'(3i4)') 1,j,nb
      do na=1,natx 
         write(6,'(3f18.10)') force(1,1,j,na,nb), force(1,2,j,na,nb), force(1,3,j,na,nb)
      end do
   end do
end do
end if
! populate the force matrix with pbc 

allocate (fforce(innx,3,3,natx,natx))

do nr1=1,nrx1
   do nr2=1,nrx2
      do nr3=1,nrx3
         do nb=1,nat_0
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
                     do na=1,nat_0
                        do j=1,3
                           do i=1,3
                              do inn=1,innx
                                 fforce(inn,i,j,inat(na,nax,nbx,ncx),nbb)=force(inn,i,j,inat(na,nrr1,nrr2,nrr3),nb)
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

deallocate (force0)
deallocate (force)

allocate (phid(3,3,natx,natx))

!compute IFCs
do j=1,3
   do na=1,natx
      do i=1,3
         do nb=1,natx
            if (innx == 2) then
               phid(i,j,na,nb)=-0.5*(fforce(1,i,j,na,nb)-fforce(2,i,j,na,nb))/de
            else
               phid(i,j,na,nb)=-1.0*fforce(1,i,j,na,nb)/de
            end if
         end do
      end do
   end do
end do

deallocate (fforce)

!symmetrization of the ifc's

allocate (phid_symm(3,3,natx,natx))

do na=1,natx
   do nb=1,natx
      do i=1,3
         do j=1,3
            phid_symm(i,j,na,nb)=0.5*(phid(i,j,na,nb)+phid(j,i,nb,na))
         end do
      end do
   end do
end do

deallocate (phid)

if (innx == 1) then

   ! symmetrize

   IF ( ALLOCATED( irt ) ) DEALLOCATE( irt )
   ALLOCATE( irt(48,natx))
   call find_sym_ifc( natx, atomx, itypx )
   call symifc (natx, phid_symm, irt)

end if

fileout=TRIM(file_out)//'.fc'
OPEN(3,FILE=TRIM(fileout))

WRITE (3,'(3i4)') nrx1,nrx2,nrx3
do i=1,3
   do j=1,3
      do na=1,nat_0
         do nb=1,nat_0
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

deallocate (phid_symm)
deallocate (inat)

IF(ld3) THEN

allocate (force3_11(2,3,3,3,natx,natx,natx))
allocate (force3_12(3,3,3,natx,natx,natx))
allocate (force3(2,3,3,3,natx,natx,natx))
allocate (phid3(3,3,3,natx,natx,natx))
allocate (phidD3(3,3,3,natx,natx,natx))
allocate (phid_symm3(3,3,3,natx,natx,natx))

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

IF (ld3) THEN
   write(6,*) 'THIS IS NOT IMPLEMENTED YET'
   stop

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

stop 'NOT IMPLEMENTED!'

END IF

! IFCs

do i=1,3
   do j=1,3
      do k=1,3
         do na=1,natx
            do nb=1,natx
               do nc=1,natx
                  phidD3(i,j,k,na,nb,nc)=(-0.5d0)*(force3(1,i,j,k,na,nb,nc)+force3(2,i,j,k,na,nb,nc))/de**2
!-0.5*(phid3(i,j,j,na,nb,nb)+phid3(i,k,k,na,nc,nc))
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

IF (verbose .and. ld3) THEN

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

deallocate (force3_11)
deallocate (force3_12)
deallocate (force3)
deallocate (phid3)
deallocate (phidD3)
deallocate (phid_symm3)

END IF


   CALL environment_end( 'FD_IFC' )

end program fd_ifc

subroutine equiv_atoms(na,nsym,example,res,pos,k)
  implicit none
  integer :: nsym
  integer :: example(nsym),na         ! The input
  integer :: res(nsym),pos(nsym)  ! The output
  integer :: k                   ! The number of unique elements
  integer :: i, j

  k = 1
  res(1) = example(1)
  pos(1) = 1
  outer: do i=2,size(example)
     do j=1,k
        if (res(j) == example(i) .or. example(i) == na) then
           ! Found a match so start looking again
           cycle outer
        end if
     end do
     ! No match found so add it to the output
     k = k + 1
     res(k) = example(i)
     pos(k) = i
  end do outer

  write(*,advance='no',fmt='(a,i0,a,i0,a)') 'Atom ',na, ' has ',k,' equivalent(s): '
  write(*,*) res(1:k)
  write(*,advance='no',fmt='(a)') 'for symmetry operation(s):  '
  write(*,*) pos(1:k)
end subroutine equiv_atoms

   SUBROUTINE symifc (nat, tens, irts)
     !-----------------------------------------------------------------------
     ! Symmetrize a function f(i,j,na,nb), i,j=cartesian components, na,nb=atom index
     !

  USE kinds,      ONLY : dp
  USE symm_base
  USE symme
  USE cell_base,  ONLY : at, bg
  
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     INTEGER :: irts(48,nat)
     REAL(DP), intent(INOUT) :: tens(3,3,nat,nat)
     !
     INTEGER :: na,nb, isym, nar,nbr, i,j,k,l
     REAL(DP), ALLOCATABLE :: work (:,:,:,:)
     !
     IF (nsym == 1) RETURN
     !
     ! bring tensor to crystal axis
     !
     DO na=1,nat
        do nb=1,nat
           CALL cart_to_crys ( tens (:,:,na,nb) )
        end do
     END DO
     !
     ! symmetrize in crystal axis
     !
     ALLOCATE (work(3,3,nat,nat))
     work (:,:,:,:) = 0.0_dp
     DO na = 1, nat
        do nb = 1, nat
           DO isym = 1, nsym
              nar = irts (isym, na)
              nbr = irts (isym, nb)
              DO i = 1, 3
                 DO j = 1, 3
                    DO k = 1, 3
                       DO l = 1, 3
                          work (i,j,na,nb) = work (i,j,na,nb) + &
                             s (i,k,isym) * s (j,l,isym) * tens (k,l,nar,nbr)
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
     tens (:,:,:,:) = work (:,:,:,:) / DBLE(nsym)
     DEALLOCATE (work)
     !
     ! bring tensor back to cartesian axis
     !
     DO na=1,nat
        do nb=1,nat
           CALL crys_to_cart ( tens (:,:,na,nb) )
        end do
     END DO
     !
     !
   END SUBROUTINE symifc

!-----------------------------------------------------------------------
logical function eqdisp (x, y, accep )
  !-----------------------------------------------------------------------
  !
  !   This function test if the difference x-y is zero
  !   x, y = 3d vectors
  !
  USE kinds
  implicit none
  real(DP), intent(in) :: x (3), y (3), accep
  !
  eqdisp = abs( x(1)-y(1) ) < accep .and. &
           abs( x(2)-y(2) ) < accep .and. &
           abs( x(3)-y(3) ) < accep
  !
  return
end function eqdisp

