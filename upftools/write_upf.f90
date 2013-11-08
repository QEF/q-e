!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE upf
  !
  ! All variables to be written into the UPF file
  ! (UPF = unified pseudopotential format, v.1)
  !
  ! pp_info
  INTEGER :: rel
  real(8) :: rcloc
  INTEGER :: nwfs
  real(8), ALLOCATABLE :: oc(:), rcut(:), rcutus(:), epseu(:)
  CHARACTER(len=2), ALLOCATABLE :: els(:)
  INTEGER, ALLOCATABLE:: lchi (:), nns (:)
  !
  ! pp_header
  CHARACTER (len=80):: generated, date_author, comment
  CHARACTER (len=2) :: psd, pseudotype
  INTEGER :: nv = 0
  INTEGER :: iexch, icorr, igcx, igcc
  INTEGER :: lmax, mesh, nbeta, ntwfc
  LOGICAL :: nlcc
  real(8) :: zp, ecutrho, ecutwfc, etotps
  real(8), ALLOCATABLE :: ocw(:)
  CHARACTER(len=2), ALLOCATABLE :: elsw(:)
  INTEGER, ALLOCATABLE:: lchiw(:)
  !
  ! pp_mesh
  real(8), ALLOCATABLE :: r(:), rab(:)
  !
  ! pp_nlcc
  real(8), ALLOCATABLE :: rho_atc(:)
  !
  ! pp_local
  real(8), ALLOCATABLE ::  vloc0(:)
  !
  ! pp_nonlocal
  ! pp_beta
  real(8), ALLOCATABLE :: betar(:,:)
  INTEGER, ALLOCATABLE:: lll(:), ikk2(:)
  ! pp_dij
  real(8), ALLOCATABLE :: dion(:,:)
  ! pp_qij
  INTEGER ::  nqf, nqlc
  real(8), ALLOCATABLE :: rinner(:), qqq(:,:), qfunc(:,:,:)
  ! pp_qfcoef
  real(8), ALLOCATABLE :: qfcoef(:,:,:,:)
  !
  ! pp_pswfc
  real(8), ALLOCATABLE :: chi(:,:)
  !
  ! pp_rhoatom
  real(8), ALLOCATABLE :: rho_at(:)
END MODULE upf
!
SUBROUTINE write_upf_v1(ounps)

  USE upf, ONLY: nlcc

  INTEGER :: ounps

  CALL write_pseudo_comment(ounps)
  CALL write_pseudo_header(ounps)
  CALL write_pseudo_mesh(ounps)
  IF (nlcc)  CALL write_pseudo_nlcc(ounps)
  CALL write_pseudo_local(ounps)
  CALL write_pseudo_nl(ounps)
  CALL write_pseudo_pswfc(ounps)
  CALL write_pseudo_rhoatom(ounps)
  !
  PRINT '("*** PLEASE TEST BEFORE USING!!! ***")'
  PRINT '("review the content of the PP_INFO fields")'
  !
END SUBROUTINE write_upf_v1

  !
  !---------------------------------------------------------------------
  SUBROUTINE write_pseudo_comment (ounps)
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the comments of the new UPF file
    !
    USE upf
    IMPLICIT NONE
    INTEGER :: ounps

    INTEGER :: nb, ios

    WRITE (ounps, '(a9)', err = 100, iostat = ios) "<PP_INFO>"

    WRITE (ounps, '(a)', err = 100, iostat = ios) generated
    WRITE (ounps, '(a)', err = 100, iostat = ios) date_author
    WRITE (ounps, '(a)', err = 100, iostat = ios) comment
    IF (rel==2) THEN
       WRITE (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,&
            &"The Pseudo was generated with a Full-Relativistic Calculation"
    ELSEIF (rel==1) THEN
       WRITE (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,&
            &"The Pseudo was generated with a Scalar-Relativistic Calculation"
    ELSEIF (rel==0) THEN
       WRITE (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel, &
            & "The Pseudo was generated with a Non-Relativistic Calculation"
    ENDIF

    IF (rcloc > 0.d0) &
       WRITE (ounps, '(1pe19.11,t24,a)', err = 100, iostat = ios) &
              rcloc, "Local Potential cutoff radius"

    IF (nwfs>0) &
       WRITE (ounps, '(a2,2a3,a6,3a19)', err = 100, iostat = ios) "nl", &
              &" pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
    DO nb = 1, nwfs
       WRITE (ounps, '(a2,2i3,f6.2,3f19.11)') els (nb) , nns (nb) , &
            lchi (nb) , oc (nb) , rcut (nb) , rcutus (nb) , epseu(nb)

    ENDDO

    WRITE (ounps, '(a10)', err = 100, iostat = ios) "</PP_INFO>"
    RETURN
100 WRITE(6,'("write_pseudo_comment: error writing pseudopotential file")')
    STOP

  END SUBROUTINE write_pseudo_comment

  !
  !---------------------------------------------------------------------
  SUBROUTINE write_pseudo_header (ounps)
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the header of the new UPF file
    !
    USE upf
    IMPLICIT NONE
    INTEGER :: ounps
    !
    CHARACTER (len=4) :: shortname
    CHARACTER (len=20):: dft
    INTEGER :: nb, ios
    !
    !
    WRITE (ounps, '(//a11)', err = 100, iostat = ios) "<PP_HEADER>"

    WRITE (ounps, '(t3,i2,t24,a)', err = 100, iostat = ios) nv, &
         "Version Number"
    WRITE (ounps, '(t3,a,t24,a)', err = 100, iostat = ios) psd , &
         "Element"
    IF (pseudotype == 'NC') THEN
       WRITE (ounps, '(a5,t24,a)', err = 100, iostat = ios) "NC", &
            "Norm - Conserving pseudopotential"
    ELSEIF (pseudotype == 'US') THEN
       WRITE (ounps, '(a5,t24,a)', err = 100, iostat = ios) "US", &
            "Ultrasoft pseudopotential"
    ELSE
       WRITE(6,'("write_pseudo_header: unknown PP type ",A)') pseudotype
       STOP
    ENDIF
    WRITE (ounps, '(l5,t24,a)', err = 100, iostat = ios) nlcc , &
         "Nonlinear Core Correction"
    CALL dftname (iexch, icorr, igcx, igcc, dft, shortname)
    WRITE (ounps, '(a,t24,a4,a)', err = 100, iostat = ios) &
         dft, shortname," Exchange-Correlation functional"
    WRITE (ounps, '(f17.11,t24,a)') zp , "Z valence"
    WRITE (ounps, '(f17.11,t24,a)') etotps, "Total energy"
    WRITE (ounps, '(2f11.5,t24,a)') ecutwfc, ecutrho, &
         "Suggested cutoff for wfc and rho"

    WRITE (ounps, '(i5,t24,a)') lmax, "Max angular momentum component"
    WRITE (ounps, '(i5,t24,a)') mesh, "Number of points in mesh"
    WRITE (ounps, '(2i5,t24,a)', err = 100, iostat = ios) ntwfc, &
         nbeta  , "Number of Wavefunctions, Number of Projectors"
    WRITE (ounps, '(a,t24,a2,a3,a6)', err = 100, iostat = ios) &
         " Wavefunctions", "nl", "l", "occ"
    DO nb = 1, ntwfc
       WRITE (ounps, '(t24,a2,i3,f6.2)') elsw(nb), lchiw(nb), ocw(nb)
    ENDDO
    !---> End header writing

    WRITE (ounps, '(a12)', err = 100, iostat = ios) "</PP_HEADER>"
    RETURN
100 WRITE(6,'("write_pseudo_header: error writing pseudopotential file")')
    STOP

  END SUBROUTINE write_pseudo_header

  !
  !---------------------------------------------------------------------
  SUBROUTINE write_pseudo_mesh (ounps)
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    USE upf
    IMPLICIT NONE
    INTEGER :: ounps
    !
    INTEGER :: ir, ios
    !
    WRITE (ounps, '(//a9)', err = 100, iostat = ios) "<PP_MESH>"

    WRITE (ounps, '(t3,a6)', err = 100, iostat = ios) "<PP_R>"
    WRITE (ounps, '(1p4e19.11)', err=100, iostat=ios) (r(ir),  ir=1,mesh )
    WRITE (ounps, '(t3,a7)', err = 100, iostat = ios) "</PP_R>"
    WRITE (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_RAB>"
    WRITE (ounps, '(1p4e19.11)', err=100, iostat=ios) (rab(ir), ir=1,mesh )
    WRITE (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_RAB>"

    WRITE (ounps, '(a10)', err = 100, iostat = ios) "</PP_MESH>"

    RETURN

100 WRITE(6,'("write_pseudo_mesh: error writing pseudopotential file")')
    STOP

  END SUBROUTINE write_pseudo_mesh

  !
  !---------------------------------------------------------------------
  SUBROUTINE write_pseudo_nlcc (ounps)
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the core charge for the nonlinear core
    !     correction of the new UPF file
    !
    USE upf
    IMPLICIT NONE
    INTEGER :: ounps
    !
    INTEGER :: ir, ios

    WRITE (ounps, '(//a9)', err = 100, iostat = ios) "<PP_NLCC>"

    WRITE (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                 ( rho_atc(ir), ir = 1, mesh )
    WRITE (ounps, '(a10)', err = 100, iostat = ios) "</PP_NLCC>"
    RETURN

100 WRITE(6,'("write_pseudo_nlcc: error writing pseudopotential file")')
    STOP

  END SUBROUTINE write_pseudo_nlcc
  !
  !---------------------------------------------------------------------
  SUBROUTINE write_pseudo_local (ounps)
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the local part of the new UPF file
    !
    USE upf
    IMPLICIT NONE
    INTEGER :: ounps
    !
    INTEGER :: ir, ios

    WRITE (ounps, '(//a10)', err = 100, iostat = ios) "<PP_LOCAL>"
    WRITE (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                ( vloc0(ir), ir = 1, mesh )
    WRITE (ounps, '(a11)', err = 100, iostat = ios) "</PP_LOCAL>"
    RETURN
100 WRITE(6,'("write_pseudo_local: error writing pseudopotential file")')
    STOP

  END SUBROUTINE write_pseudo_local
  !
  !---------------------------------------------------------------------
  SUBROUTINE write_pseudo_nl (ounps)
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the non local part of the new UPF file
    !
    USE upf
    IMPLICIT NONE
    INTEGER :: ounps
    !
    INTEGER :: nb, mb, n, ir, nd, i, lp, ios

    WRITE (ounps, '(//a13)', err = 100, iostat = ios) "<PP_NONLOCAL>"
    DO nb = 1, nbeta
       WRITE (ounps, '(t3,a9)', err = 100, iostat = ios) "<PP_BETA>"
       WRITE (ounps, '(2i5,t24,a)', err=100, iostat=ios) &
                                    nb, lll(nb), "Beta    L"
       WRITE (ounps, '(i6)', err=100, iostat=ios) ikk2 (nb)
       WRITE (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                    ( betar(ir,nb), ir=1,ikk2(nb) )
       WRITE (ounps, '(t3,a10)', err = 100, iostat = ios) "</PP_BETA>"
    ENDDO

    WRITE (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_DIJ>"
    nd = 0
    DO nb = 1, nbeta
       DO mb = nb, nbeta
          IF ( abs(dion(nb,mb)) > 1.0d-12 )  nd = nd + 1
       ENDDO
    ENDDO
    WRITE (ounps, '(1p,i5,t24,a)', err=100, iostat=ios) &
                                   nd, "Number of nonzero Dij"
    DO nb = 1, nbeta
       DO mb = nb, nbeta
          IF ( abs(dion(nb,mb)) > 1.0d-12 ) &
             WRITE(ounps,'(1p,2i5,e19.11)', err=100, iostat=ios) &
                                   nb, mb, dion(nb,mb)
       ENDDO
    ENDDO
    WRITE (ounps, '(t3,a9)', err=100, iostat=ios) "</PP_DIJ>"

    IF (pseudotype == 'US') THEN
       WRITE (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_QIJ>"
       WRITE (ounps, '(i5,a)',err=100, iostat=ios) nqf,"     nqf.&
          & If not zero, Qij's inside rinner are computed using qfcoef's"
       IF (nqf>0) THEN
          WRITE (ounps, '(t5,a11)', err=100, iostat=ios) "<PP_RINNER>"
          WRITE (ounps,'(i5,1pe19.11)', err=100, iostat=ios) &
                                        (i, rinner(i), i = 1, nqlc)
          WRITE (ounps, '(t5,a12)', err=100, iostat=ios) "</PP_RINNER>"
       ENDIF
       DO nb = 1, nbeta
          DO mb = nb, nbeta
             WRITE (ounps, '(3i5,t24,a)', err=100, iostat=ios) &
                                          nb, mb, lll(mb) , "i  j  (l(j))"
             WRITE (ounps, '(1pe19.11,t24,a)', err=100, iostat=ios) &
                                          qqq(nb,mb), "Q_int"
             WRITE (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                          ( qfunc (n,nb,mb), n=1,mesh )
             IF (nqf>0) THEN
                WRITE (ounps, '(t5,a11)', err=100, iostat=ios) &
                                          "<PP_QFCOEF>"
                WRITE(ounps,'(1p4e19.11)', err=100, iostat=ios) &
                                 ((qfcoef(i,lp,nb,mb),i=1,nqf),lp=1,nqlc)
                WRITE (ounps, '(t5,a12)', err=100, iostat=ios) &
                                          "</PP_QFCOEF>"
             ENDIF
          ENDDO
       ENDDO
       WRITE (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_QIJ>"

    ENDIF
    WRITE (ounps, '(a14)', err = 100, iostat = ios) "</PP_NONLOCAL>"
    RETURN

100 WRITE(6,'("write_pseudo_nl: error writing pseudopotential file")')
    STOP

  END SUBROUTINE write_pseudo_nl

  !
  !---------------------------------------------------------------------
  SUBROUTINE write_pseudo_pswfc (ounps)
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the pseudo atomic functions
    !     of the new UPF file
    !
    USE upf
    IMPLICIT NONE
    INTEGER :: ounps
    !
    INTEGER :: nb, ir, ios

    WRITE (ounps, '(//a10)', err = 100, iostat = ios) "<PP_PSWFC>"
    DO nb = 1, ntwfc
       WRITE (ounps,'(a2,i5,f6.2,t24,a)', err=100, iostat=ios) &
            elsw(nb), lchiw(nb), ocw(nb), "Wavefunction"
       WRITE (ounps, '(1p4e19.11)', err=100, iostat=ios) &
            ( chi(ir,nb), ir=1,mesh )
    ENDDO
    WRITE (ounps, '(a11)', err = 100, iostat = ios) "</PP_PSWFC>"
    RETURN

100 WRITE(6,'("write_pseudo_pswfc: error writing pseudopotential file")')
    STOP

  END SUBROUTINE write_pseudo_pswfc
  !
  !---------------------------------------------------------------------
  SUBROUTINE write_pseudo_rhoatom (ounps)
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    USE upf
    IMPLICIT NONE
    INTEGER :: ounps
    !
    INTEGER :: ir, ios

    WRITE (ounps, '(//a12)', err = 100, iostat = ios) "<PP_RHOATOM>"
    WRITE (ounps, '(1p4e19.11)', err = 100, iostat = ios) &
                               ( rho_at(ir), ir=1,mesh )
    WRITE (ounps, '(a13)', err = 100, iostat = ios) "</PP_RHOATOM>"
    RETURN

100 WRITE(6,'("write_pseudo_rhoatom: error writing pseudopotential file")')
    STOP

  END SUBROUTINE write_pseudo_rhoatom

  !---------------------------------------------------------------------
  SUBROUTINE dftname(iexch, icorr, igcx, igcc, longname, shortname)
  !---------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER iexch, icorr, igcx, igcc
  CHARACTER (len=4) :: shortname
  CHARACTER (len=20):: longname
  !
  ! The data used to convert iexch, icorr, igcx, igcc
  ! into a user-readable string
  !
  integer :: nxc, ncc, ngcx, ngcc, ncnl

  parameter (nxc = 8, ncc =11, ngcx =19, ngcc = 12)

  character (len=4) :: exc, corr
  character (len=4) :: gradx, gradc
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0: ngcc)

  data exc / 'NOX', 'SLA', 'SL1', 'RXC', 'OEP', 'HF', 'PB0X', 'B3LP', 'KZK' /
  data corr / 'NOC', 'PZ', 'VWN', 'LYP', 'PW', 'WIG', 'HL', 'OBZ', &
              'OBW', 'GL' , 'B3LP', 'KZK' /

  data gradx / 'NOGX', 'B88', 'GGX', 'PBX',  'RPB', 'HCTH', 'OPTX',&
               'TPSS', 'PB0X', 'B3LP','PSX', 'WCX', 'HSE', 'RW86', 'PBE', &
               'META', 'C09X', 'SOX', 'M6LX', 'Q2DX' /

  data gradc / 'NOGC', 'P86', 'GGC', 'BLYP', 'PBC', 'HCTH', 'TPSS',&
                'B3LP', 'PSC', 'PBE', 'META', 'M6LC', 'Q2DC' /

  IF (iexch==1.and.igcx==0.and.igcc==0) THEN
     shortname = corr(icorr)
  ELSEIF (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) THEN
     shortname = 'BLYP'
  ELSEIF (iexch==1.and.icorr==1.and.igcx==1.and.igcc==0) THEN
     shortname = 'B88'
  ELSEIF (iexch==1.and.icorr==1.and.igcx==1.and.igcc==1) THEN
     shortname = 'BP'
  ELSEIF (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) THEN
     shortname = 'PW91'
  ELSEIF (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) THEN
     shortname = 'PBE'
  ELSEIF (iexch==1.and.icorr==4.and.igcx==4.and.igcc==5) THEN
     shortname = 'TPSS'
  ELSEIF (iexch==1.and.icorr==4.and.igcx==10.and.igcc==8) THEN
     shortname = 'PBESOL'
  ELSE
     shortname = ' '
  ENDIF
  WRITE(longname,'(4a5)') exc(iexch),corr(icorr),gradx(igcx),gradc(igcc)

  RETURN
END SUBROUTINE dftname
