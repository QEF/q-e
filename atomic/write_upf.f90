!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine write_upf(ounps)

  use ld1inc, only: nlcc, rel

  integer :: ounps

  call write_pseudo_comment(ounps)  
  call write_pseudo_header(ounps)  
  call write_pseudo_mesh(ounps)
  if (nlcc)  call write_pseudo_nlcc(ounps)  
  call write_pseudo_local(ounps)  
  call write_pseudo_nl(ounps)  
  call write_pseudo_pswfc(ounps)
  call write_pseudo_rhoatom(ounps)  
  if (rel == 2) call write_pseudo_addinfo(ounps)  
  !
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
    use ld1inc
    implicit none
    integer :: ounps  

    integer :: nb, ios  
    character(len=80) :: generated, date_author, comment
    character(len=9) :: day, hour

    call date_and_tim(day,hour)
    generated='Generated using Andrea Dal Corso code (rrkj3)'
    date_author='Andrea Dal Corso   Generation date: '// day 
    comment=title

    write (ounps, '(a9)', err = 100, iostat = ios) "<PP_INFO>"  

    write (ounps, '(a)', err = 100, iostat = ios) generated
    write (ounps, '(a)', err = 100, iostat = ios) date_author
    write (ounps, '(a)', err = 100, iostat = ios) comment
    if (rel==2) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Fully-Relativistic Calculation"
    else if (rel==1) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Scalar-Relativistic Calculation"
    else  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel, &
            & "The Pseudo was generated with a Non-Relativistic Calculation"
    endif

    write (ounps, '(1pe19.11,t24,a)', err = 100, iostat = ios) &
         rcloc, "Local Potential cutoff radius"

    write (ounps, '(a2,2a3,a6,3a19)', err = 100, iostat = ios) "nl", &
         &" pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
    do nb = 1, nwfs  
          write (ounps, '(a2,2i3,f6.2,3f19.11)') els (nb) , nns (nb) , &
            lls (nb) , ocs (nb) , rcut (nb) , rcutus (nb) , enls(nb)
    enddo

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_INFO>"  
    return
100 call errore ('write_pseudo_comment', 'Writing pseudo file', abs ( &
         ios))   
  end subroutine write_pseudo_comment

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_header (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the header of the new UPF file
    !
    use ld1inc
    use funct
    use kinds, only : DP
    implicit none
    integer :: ounps  
    !
    character (len=4) :: shortname
    real(kind=dp) ecutrho, ecutwfc
    integer :: nb, ios, nv  
    !
    !
    write (ounps, '(//a11)', err = 100, iostat = ios) "<PP_HEADER>"  

    nv=0
    write (ounps, '(t3,i2,t24,a)', err = 100, iostat = ios) nv, &
         "Version Number"
    write (ounps, '(t3,a,t24,a)', err = 100, iostat = ios) title , &
         "Element"
    if (pseudotype == 1.or.pseudotype == 2) then  
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "NC", &
            "Norm - Conserving pseudopotential"
    else if (pseudotype == 3) then
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "US", &
            "Ultrasoft pseudopotential"
    else
       call errore ('write_pseudo_header',&
            'Unknown PP type: ', 1)
    endif
    write (ounps, '(l5,t24,a)', err = 100, iostat = ios) nlcc , &
         "Nonlinear Core Correction"
    call dft_name (iexch, icorr, igcx, igcc, dft, shortname)
    write (ounps, '(a,t24,a4,a)', err = 100, iostat = ios) &
         dft, shortname," Exchange-Correlation functional"
    write (ounps, '(f17.11,t24,a)') zval , "Z valence"  
    write (ounps, '(f17.11,t24,a)') etots, "Total energy"  
    ecutrho=0.0_DP
    ecutwfc=0.0_DP
    write (ounps, '(2f11.7,t24,a)') ecutrho, ecutwfc, &
         "Suggested cutoff for wfc and rho"  

    write (ounps, '(i5,t24,a)') lmax, "Max angular momentum component"  
    write (ounps, '(i5,t24,a)') mesh, "Number of points in mesh"
    write (ounps, '(2i5,t24,a)', err = 100, iostat = ios) nwfts, &
         nbeta  , "Number of Wavefunctions, Number of Projectors"
    write (ounps, '(a,t24,a2,a3,a6)', err = 100, iostat = ios) &
         " Wavefunctions", "nl", "l", "occ"
    do nb = 1, nwfts
       write (ounps, '(t24,a2,i3,f6.2)') elts(nb), llts(nb), octs(nb)
    enddo
    !---> End header writing

    write (ounps, '(a12)', err = 100, iostat = ios) "</PP_HEADER>"
    return   
100 call errore ('write_pseudo_header','Writing pseudo file', abs(ios) )

  end subroutine write_pseudo_header

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_mesh (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use ld1inc
    implicit none
    integer :: ounps  
    real(kind=dp) :: rab(ndm)
    !
    integer :: ir, ios  
    !
    rab=dx*r
    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_MESH>"  

    write (ounps, '(t3,a6)', err = 100, iostat = ios) "<PP_R>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (r(ir),  ir=1,mesh )
    write (ounps, '(t3,a7)', err = 100, iostat = ios) "</PP_R>"  
    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_RAB>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (rab(ir), ir=1,mesh )
    write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_RAB>"  

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_MESH>"  

    return

100 call errore ('write_pseudo_rhoatom','Writing pseudo file',abs(ios))

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
    use ld1inc
    use kinds, only : DP
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  
    real(kind=dp) :: fpi

    fpi=16.0_DP*atan(1.0_DP)

    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_NLCC>"  

    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                 ( rhoc(ir)/fpi/r2(ir), ir = 1, mesh )
    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_NLCC>"  
    return

100 call errore ('write_pseudo_nlcc', 'Writing pseudo file', abs (ios))

  end subroutine write_pseudo_nlcc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_local (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the local part of the new UPF file
    !
    use ld1inc
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_LOCAL>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                ( vpsloc(ir), ir = 1, mesh )
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_LOCAL>"  
    return
100 call errore ('write_pseudo_local', 'Writing pseudo file', abs(ios) )  
  end subroutine write_pseudo_local

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nl (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the non local part of the new UPF file
    !
    use ld1inc
    implicit none
    integer :: ounps  
    !
    integer :: nb, mb, n, ir, nd, i, lp, nqf, ios  

    write (ounps, '(//a13)', err = 100, iostat = ios) "<PP_NONLOCAL>"  
    do nb = 1, nbeta  
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "<PP_BETA>"  
       write (ounps, '(2i5,t24,a)', err=100, iostat=ios) &
                                       nb, lls(nb), "Beta    L"
       write (ounps, '(i6)', err=100, iostat=ios) ikk (nb)  
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                    ( betas(ir,nb), ir=1,ikk(nb) )
       write (ounps, '(t3,a10)', err = 100, iostat = ios) "</PP_BETA>"  
    enddo

    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_DIJ>"  
    nd = 0  
    do nb = 1, nbeta  
       do mb = nb, nbeta  
          if ( abs(bmat(nb,mb)) .gt. 1.0e-12_dp )  nd = nd + 1 
       enddo
    enddo
    write (ounps, '(1p,i5,t24,a)', err=100, iostat=ios) &
                                   nd, "Number of nonzero Dij"
    do nb = 1, nbeta
       do mb = nb, nbeta  
          if ( abs(bmat(nb,mb)) .gt. 1.0e-12_dp ) &
             write(ounps,'(1p,2i5,e19.11)', err=100, iostat=ios) &
                                   nb, mb, bmat(nb,mb)
       enddo
    enddo
    write (ounps, '(t3,a9)', err=100, iostat=ios) "</PP_DIJ>"  

    if (pseudotype == 3) then  
       write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_QIJ>"  
       nqf=0
       write (ounps, '(i5,a)',err=100, iostat=ios) nqf,"     nqf.&
          & If not zero, Qij's inside rinner are computed using qfcoef's"
       do nb = 1, nbeta 
          do mb = nb, nbeta
             write (ounps, '(3i5,t24,a)', err=100, iostat=ios) &
                                          nb, mb, lls(mb) , "i  j  (l(j))"
             write (ounps, '(1pe19.11,t24,a)', err=100, iostat=ios) &
                                          qq(nb,mb), "Q_int"
             write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                          ( qvan (n,nb,mb), n=1,mesh )
          enddo
       enddo
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_QIJ>"  

    endif
    write (ounps, '(a14)', err = 100, iostat = ios) "</PP_NONLOCAL>"  
    return

100 call errore ('write_pseudo_nl', 'Writing pseudo file', abs (ios) )  

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
    use ld1inc
    implicit none
    integer :: ounps  
    !
    integer :: nb, ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_PSWFC>"  
    do nb = 1, nwfts
       write (ounps,'(a2,i5,f6.2,t24,a)', err=100, iostat=ios) &
            elts(nb), llts(nb), octs(nb), "Wavefunction"
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
            ( phits(ir,nb), ir=1,mesh )
    enddo
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_PSWFC>"  
    return

100 call errore ('write_pseudo_pswfc', 'Writing pseudo file', abs(ios) )  
  end subroutine write_pseudo_pswfc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_rhoatom (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use ld1inc
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a12)', err = 100, iostat = ios) "<PP_RHOATOM>"  
    write (ounps, '(1p4e19.11)', err = 100, iostat = ios) &
                               ( rhos(ir,1), ir=1,mesh )
    write (ounps, '(a13)', err = 100, iostat = ios) "</PP_RHOATOM>"  
    return

100 call errore('write_pseudo_rhoatom','Writing pseudo file',abs(ios))
  end subroutine write_pseudo_rhoatom
!---------------------------------------------------------------------
  subroutine write_pseudo_addinfo (ounps)  
!---------------------------------------------------------------------
!
!     This routine writes the atomic charge density to the new UPF file
!
    use ld1inc
    implicit none
    integer :: ounps  
    !
    integer :: nb, ios  

    write (ounps, '(//a12)', err = 100, iostat = ios) "<PP_ADDINFO>"  
    do nb=1,nwfts
       write (ounps, '(a2,2i3,2f6.2)',err=100, iostat=ios) &
    &       elts(nb), nnts(nb), llts (nb), jjts(nb), octs (nb)
    enddo
    do nb=1,nbeta
       write (ounps, '(i5,f6.2)', err=100,iostat=ios) lls(nb), jjs(nb)
    enddo

    write(ounps, '(4f15.8)') xmin, rmax, zmesh, dx

    write (ounps, '(a13)', err = 100, iostat = ios) "</PP_ADDINFO>"  
100 call errore('write_pseudo_addinfo','Writing pseudo file',abs(ios))

end subroutine write_pseudo_addinfo

  !---------------------------------------------------------------------
  subroutine dft_name(iexch, icorr, igcx, igcc, longname, shortname)
  !---------------------------------------------------------------------
  implicit none
  integer iexch, icorr, igcx, igcc
  character (len=4) :: shortname
  character (len=20):: longname
  !
  ! The data used to convert iexch, icorr, igcx, igcc
  ! into a user-readable string
  !
  integer, parameter :: nxc = 3, ncc = 9, ngcx = 4, ngcc = 4 
  character (len=4) :: exc, corr, gradx, gradc  
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0:ngcc)
  data exc / 'NOX ', 'SLA ', 'SL1 ', 'RXC' /  
  data corr / 'NOC ', 'PZ  ', 'VWN ', 'LYP ', 'PW  ', 'WIG ', 'HL  ',&
              'OBZ ', 'OBW ', 'GL  ' /
  data gradx / 'NOGX', 'B88 ', 'GGX ', 'PBX ', 'RPB' /  
  data gradc / 'NOGC', 'P86 ', 'GGC ', 'BLYP', 'PBC ' /  

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
  else
     shortname = ' '
  end if
  write(longname,'(4a5)') exc(iexch),corr(icorr),gradx(igcx),gradc(igcc)

  return
end subroutine dft_name

