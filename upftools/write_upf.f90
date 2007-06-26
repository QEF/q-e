!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module upf
  !
  ! All variables to be written into the UPF file
  ! (UPF = unified pseudopotential format)
  !
  ! pp_info
  integer :: rel
  real(8) :: rcloc
  integer :: nwfs
  real(8), allocatable :: oc(:), rcut(:), rcutus(:), epseu(:)
  character(len=2), allocatable :: els(:)
  integer, allocatable:: lchi (:), nns (:)
  !
  ! pp_header
  character (len=80):: generated, date_author, comment
  character (len=2) :: psd, pseudotype
  integer :: nv = 0
  integer :: iexch, icorr, igcx, igcc
  integer :: lmax, mesh, nbeta, ntwfc
  logical :: nlcc
  real(8) :: zp, ecutrho, ecutwfc, etotps
  real(8), allocatable :: ocw(:)
  character(len=2), allocatable :: elsw(:)
  integer, allocatable:: lchiw(:)
  !
  ! pp_mesh
  real(8), allocatable :: r(:), rab(:)
  !
  ! pp_nlcc
  real(8), allocatable :: rho_atc(:)
  !
  ! pp_local
  real(8), allocatable ::  vloc0(:)
  !
  ! pp_nonlocal
  ! pp_beta
  real(8), allocatable :: betar(:,:)
  integer, allocatable:: lll(:), ikk2(:)  
  ! pp_dij
  real(8), allocatable :: dion(:,:)
  ! pp_qij
  integer ::  nqf, nqlc
  real(8), allocatable :: rinner(:), qqq(:,:), qfunc(:,:,:)
  ! pp_qfcoef
  real(8), allocatable :: qfcoef(:,:,:,:)
  !
  ! pp_pswfc
  real(8), allocatable :: chi(:,:)
  !
  ! pp_rhoatom
  real(8), allocatable :: rho_at(:)
end module upf
!
subroutine write_upf(ounps)

  use upf, only: nlcc

  integer :: ounps

  call write_pseudo_comment(ounps)  
  call write_pseudo_header(ounps)  
  call write_pseudo_mesh(ounps)
  if (nlcc)  call write_pseudo_nlcc(ounps)  
  call write_pseudo_local(ounps)  
  call write_pseudo_nl(ounps)  
  call write_pseudo_pswfc(ounps)
  call write_pseudo_rhoatom(ounps)  
  !
  print '("*** PLEASE TEST BEFORE USING!!! ***")'
  print '("review the content of the PP_INFO fields")'
  !
end subroutine write_upf

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_comment (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the comments of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  

    integer :: nb, ios  

    write (ounps, '(a9)', err = 100, iostat = ios) "<PP_INFO>"  

    write (ounps, '(a)', err = 100, iostat = ios) generated
    write (ounps, '(a)', err = 100, iostat = ios) date_author
    write (ounps, '(a)', err = 100, iostat = ios) comment
    if (rel==2) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Full-Relativistic Calculation"
    else if (rel==1) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Scalar-Relativistic Calculation"
    else if (rel==0) then 
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel, &
            & "The Pseudo was generated with a Non-Relativistic Calculation"
    endif

    if (rcloc > 0.d0) &
       write (ounps, '(1pe19.11,t24,a)', err = 100, iostat = ios) &
              rcloc, "Local Potential cutoff radius"

    if (nwfs>0) &
       write (ounps, '(a2,2a3,a6,3a19)', err = 100, iostat = ios) "nl", &
              &" pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
    do nb = 1, nwfs  
       write (ounps, '(a2,2i3,f6.2,3f19.11)') els (nb) , nns (nb) , &
            lchi (nb) , oc (nb) , rcut (nb) , rcutus (nb) , epseu(nb)

    enddo

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_INFO>"  
    return
100 write(6,'("write_pseudo_comment: error writing pseudopotential file")')
    stop

  end subroutine write_pseudo_comment

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_header (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the header of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    character (len=4) :: shortname
    character (len=20):: dft  
    integer :: nb, ios  
    !
    !
    write (ounps, '(//a11)', err = 100, iostat = ios) "<PP_HEADER>"  

    write (ounps, '(t3,i2,t24,a)', err = 100, iostat = ios) nv, &
         "Version Number"
    write (ounps, '(t3,a,t24,a)', err = 100, iostat = ios) psd , &
         "Element"
    if (pseudotype == 'NC') then  
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "NC", &
            "Norm - Conserving pseudopotential"
    else if (pseudotype == 'US') then
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "US", &
            "Ultrasoft pseudopotential"
    else
       write(6,'("write_pseudo_header: unknown PP type ",A)') pseudotype
       stop
    endif
    write (ounps, '(l5,t24,a)', err = 100, iostat = ios) nlcc , &
         "Nonlinear Core Correction"
    call dftname (iexch, icorr, igcx, igcc, dft, shortname)
    write (ounps, '(a,t24,a4,a)', err = 100, iostat = ios) &
         dft, shortname," Exchange-Correlation functional"
    write (ounps, '(f17.11,t24,a)') zp , "Z valence"  
    write (ounps, '(f17.11,t24,a)') etotps, "Total energy"  
    write (ounps, '(2f11.7,t24,a)') ecutrho, ecutwfc, &
         "Suggested cutoff for wfc and rho"  

    write (ounps, '(i5,t24,a)') lmax, "Max angular momentum component"  
    write (ounps, '(i5,t24,a)') mesh, "Number of points in mesh"
    write (ounps, '(2i5,t24,a)', err = 100, iostat = ios) ntwfc, &
         nbeta  , "Number of Wavefunctions, Number of Projectors"
    write (ounps, '(a,t24,a2,a3,a6)', err = 100, iostat = ios) &
         " Wavefunctions", "nl", "l", "occ"
    do nb = 1, ntwfc
       write (ounps, '(t24,a2,i3,f6.2)') elsw(nb), lchiw(nb), ocw(nb)
    enddo
    !---> End header writing

    write (ounps, '(a12)', err = 100, iostat = ios) "</PP_HEADER>"
    return
100 write(6,'("write_pseudo_header: error writing pseudopotential file")')
    stop

  end subroutine write_pseudo_header

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_mesh (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  
    !
    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_MESH>"  

    write (ounps, '(t3,a6)', err = 100, iostat = ios) "<PP_R>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (r(ir),  ir=1,mesh )
    write (ounps, '(t3,a7)', err = 100, iostat = ios) "</PP_R>"  
    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_RAB>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (rab(ir), ir=1,mesh )
    write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_RAB>"  

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_MESH>"  

    return

100 write(6,'("write_pseudo_mesh: error writing pseudopotential file")')
    stop

  end subroutine write_pseudo_mesh

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nlcc (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the core charge for the nonlinear core
    !     correction of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_NLCC>"  

    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                 ( rho_atc(ir), ir = 1, mesh )
    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_NLCC>"  
    return

100 write(6,'("write_pseudo_nlcc: error writing pseudopotential file")')
    stop

  end subroutine write_pseudo_nlcc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_local (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the local part of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_LOCAL>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                ( vloc0(ir), ir = 1, mesh )
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_LOCAL>"  
    return
100 write(6,'("write_pseudo_local: error writing pseudopotential file")')
    stop

  end subroutine write_pseudo_local
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nl (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the non local part of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: nb, mb, n, ir, nd, i, lp, ios  

    write (ounps, '(//a13)', err = 100, iostat = ios) "<PP_NONLOCAL>"  
    do nb = 1, nbeta  
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "<PP_BETA>"  
       write (ounps, '(2i5,t24,a)', err=100, iostat=ios) &
                                    nb, lll(nb), "Beta    L"
       write (ounps, '(i6)', err=100, iostat=ios) ikk2 (nb)  
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                    ( betar(ir,nb), ir=1,ikk2(nb) )
       write (ounps, '(t3,a10)', err = 100, iostat = ios) "</PP_BETA>"  
    enddo

    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_DIJ>"  
    nd = 0  
    do nb = 1, nbeta  
       do mb = nb, nbeta  
          if ( abs(dion(nb,mb)) .gt. 1.0d-12 )  nd = nd + 1 
       enddo
    enddo
    write (ounps, '(1p,i5,t24,a)', err=100, iostat=ios) &
                                   nd, "Number of nonzero Dij"
    do nb = 1, nbeta
       do mb = nb, nbeta  
          if ( abs(dion(nb,mb)) .gt. 1.0d-12 ) &
             write(ounps,'(1p,2i5,e19.11)', err=100, iostat=ios) &
                                   nb, mb, dion(nb,mb)
       enddo
    enddo
    write (ounps, '(t3,a9)', err=100, iostat=ios) "</PP_DIJ>"  

    if (pseudotype == 'US') then  
       write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_QIJ>"  
       write (ounps, '(i5,a)',err=100, iostat=ios) nqf,"     nqf.&
          & If not zero, Qij's inside rinner are computed using qfcoef's"
       if (nqf.gt.0) then
          write (ounps, '(t5,a11)', err=100, iostat=ios) "<PP_RINNER>"  
          write (ounps,'(i5,1pe19.11)', err=100, iostat=ios) &
                                        (i, rinner(i), i = 1, nqlc)
          write (ounps, '(t5,a12)', err=100, iostat=ios) "</PP_RINNER>"  
       end if
       do nb = 1, nbeta 
          do mb = nb, nbeta
             write (ounps, '(3i5,t24,a)', err=100, iostat=ios) &
                                          nb, mb, lll(mb) , "i  j  (l(j))"
             write (ounps, '(1pe19.11,t24,a)', err=100, iostat=ios) &
                                          qqq(nb,mb), "Q_int"
             write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                          ( qfunc (n,nb,mb), n=1,mesh )
             if (nqf.gt.0) then
                write (ounps, '(t5,a11)', err=100, iostat=ios) &
                                          "<PP_QFCOEF>"  
                write(ounps,'(1p4e19.11)', err=100, iostat=ios) &
                                 ((qfcoef(i,lp,nb,mb),i=1,nqf),lp=1,nqlc)
                write (ounps, '(t5,a12)', err=100, iostat=ios) &
                                          "</PP_QFCOEF>"
             end if
          enddo
       enddo
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_QIJ>"  

    endif
    write (ounps, '(a14)', err = 100, iostat = ios) "</PP_NONLOCAL>"  
    return

100 write(6,'("write_pseudo_nl: error writing pseudopotential file")')
    stop

  end subroutine write_pseudo_nl

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_pswfc (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the pseudo atomic functions
    !     of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: nb, ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_PSWFC>"  
    do nb = 1, ntwfc
       write (ounps,'(a2,i5,f6.2,t24,a)', err=100, iostat=ios) &
            elsw(nb), lchiw(nb), ocw(nb), "Wavefunction"
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
            ( chi(ir,nb), ir=1,mesh )
    enddo
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_PSWFC>"  
    return

100 write(6,'("write_pseudo_pswfc: error writing pseudopotential file")')
    stop

  end subroutine write_pseudo_pswfc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_rhoatom (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a12)', err = 100, iostat = ios) "<PP_RHOATOM>"  
    write (ounps, '(1p4e19.11)', err = 100, iostat = ios) &
                               ( rho_at(ir), ir=1,mesh )
    write (ounps, '(a13)', err = 100, iostat = ios) "</PP_RHOATOM>"  
    return

100 write(6,'("write_pseudo_rhoatom: error writing pseudopotential file")')
    stop

  end subroutine write_pseudo_rhoatom

  !---------------------------------------------------------------------
  subroutine dftname(iexch, icorr, igcx, igcc, longname, shortname)
  !---------------------------------------------------------------------
  implicit none
  integer iexch, icorr, igcx, igcc
  character (len=4) :: shortname
  character (len=20):: longname
  !
  ! The data used to convert iexch, icorr, igcx, igcc
  ! into a user-readable string
  !
  integer, parameter :: nxc = 1, ncc = 9, ngcx = 4, ngcc = 5
  character (len=20) :: exc, corr, gradx, gradc  
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0:ngcc)
  data exc / 'NOX ', 'SLA ' /  
  data corr / 'NOC ', 'PZ  ', 'VWN ', 'LYP ', 'PW  ', 'WIG ', 'HL  ',&
              'OBZ ', 'OBW ', 'GL  ' /
  data gradx / 'NOGX', 'B88 ', 'GGX ', 'PBE ', 'TPSS' /  
  data gradc / 'NOGC', 'P86 ', 'GGC ', 'BLYP', 'PBE ', 'TPSS' /  

  if (iexch==1.and.igcx==0.and.igcc==0) then
     shortname = corr(icorr)
  else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
     shortname = 'BLYP'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==0) then
     shortname = 'B88'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==1) then
     shortname = 'BP'
  else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
     shortname = 'PW91'
  else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
     shortname = 'PBE'
  else if (iexch==1.and.icorr==4.and.igcx==4.and.igcc==5) then
     shortname = 'TPSS'
  else
     shortname = ' '
  end if
  write(longname,'(4a5)') exc(iexch),corr(icorr),gradx(igcx),gradc(igcc)

  return
end subroutine dftname
