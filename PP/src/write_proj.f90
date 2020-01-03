!
! Copyright (C) 2001-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE  write_proj_iotk (filename, lbinary, projs, lwrite_ovp, ovps )
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE io_files,         ONLY : iun => iunsat, prefix, tmp_dir, postfix
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : wk, xk, nkstot, nelec
  USE noncollin_module, ONLY : noncolin
  USE lsda_mod,         ONLY : nspin, isk
  USE ener,             ONLY : ef
  USE wvfct,            ONLY : et, nbnd
  USE iotk_module
  IMPLICIT NONE

  CHARACTER(*),  INTENT(IN) :: filename
  LOGICAL,       INTENT(IN) :: lbinary
  COMPLEX(DP),   INTENT(IN) :: projs(natomwfc,nbnd,nkstot)
  LOGICAL,       INTENT(IN) :: lwrite_ovp
  COMPLEX(DP),   INTENT(IN) :: ovps(natomwfc,natomwfc,nkstot)
  !
  CHARACTER(256)          :: tmp
  CHARACTER(iotk_attlenx) :: attr
  INTEGER :: ik, ik_eff, isp, ia, ierr, num_k_points

!
! subroutine body
!

  tmp = trim( tmp_dir ) // trim( prefix ) // postfix //trim(filename)
  !
  IF ( lbinary ) THEN
      tmp = TRIM(tmp) // ".dat"
  ELSE
      tmp = TRIM(tmp) // ".xml"
  ENDIF
  !
  CALL iotk_open_write(iun, FILE=trim(tmp), ROOT="ATOMIC_PROJECTIONS", &
                       BINARY=lbinary, IERR=ierr )
  IF ( ierr /= 0 ) RETURN
  !
  !
  num_k_points = nkstot
  IF ( nspin == 2 ) num_k_points = nkstot / 2
  !
  CALL iotk_write_begin(iun, "HEADER")
  !
  CALL iotk_write_dat(iun, "NUMBER_OF_BANDS", nbnd)
  CALL iotk_write_dat(iun, "NUMBER_OF_K-POINTS", num_k_points )
  CALL iotk_write_dat(iun, "NUMBER_OF_SPIN_COMPONENTS", nspin)
  CALL iotk_write_dat(iun, "NON-COLINEAR_CALCULATION",noncolin)
  CALL iotk_write_dat(iun, "NUMBER_OF_ATOMIC_WFC", natomwfc)
  CALL iotk_write_dat(iun, "NUMBER_OF_ELECTRONS", nelec )
  CALL iotk_write_attr(attr, "UNITS", " 2 pi / a", FIRST=.true.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_K-POINTS", ATTR=attr)
  CALL iotk_write_attr(attr, "UNITS", "Rydberg", FIRST=.true.  )
  CALL iotk_write_empty (iun,  "UNITS_FOR_ENERGY", ATTR=attr)
  CALL iotk_write_dat(iun, "FERMI_ENERGY", ef )
  !
  CALL iotk_write_end(iun, "HEADER")
  !
  !
  CALL iotk_write_dat(iun, "K-POINTS", xk(:,1:num_k_points) , COLUMNS=3 )
  CALL iotk_write_dat(iun, "WEIGHT_OF_K-POINTS", wk(1:num_k_points), COLUMNS=8 )
  !
  CALL iotk_write_begin(iun, "EIGENVALUES")
  !
  DO ik=1,num_k_points
     !
     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_dat( iun, "EIG.1", et(:,ik) )
        CALL iotk_write_dat( iun, "EIG.2", et(:,ik_eff) )
        !
     ELSE
        !
        CALL iotk_write_dat( iun, "EIG", et(:,ik) )
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "EIGENVALUES")

  !
  ! main loop atomic wfc
  !
  CALL iotk_write_begin(iun, "PROJECTIONS")
  !
  DO ik=1,num_k_points
     !
     CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
     IF ( nspin == 2 ) THEN
        !
        CALL iotk_write_begin ( iun, "SPIN.1" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.1" )
        !
        ik_eff = ik + num_k_points
        !
        CALL iotk_write_begin ( iun, "SPIN.2" )
           !
           DO ia = 1, natomwfc
               CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik_eff)  )
           ENDDO
           !
        CALL iotk_write_end ( iun, "SPIN.2" )
        !
     ELSE
        !
        DO ia = 1,natomwfc
            CALL iotk_write_dat(iun, "ATMWFC"//trim(iotk_index(ia)), projs(ia,:,ik)  )
        ENDDO
        !
     ENDIF
     !
     CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
     !
  ENDDO
  !
  CALL iotk_write_end(iun, "PROJECTIONS")

  !
  ! overlaps
  !
  IF ( lwrite_ovp ) THEN
      !
      CALL iotk_write_begin(iun, "OVERLAPS")
      !
      DO ik=1,num_k_points
          !
          CALL iotk_write_begin( iun, "K-POINT"//trim(iotk_index(ik)) )
          !
          DO isp = 1, nspin
              !
              ik_eff = ik + num_k_points * ( isp -1 )
              !
              CALL iotk_write_dat(iun, "OVERLAP"//trim(iotk_index(isp)), ovps(:,:,ik_eff)  )
              !
              !
          ENDDO
          !
          CALL iotk_write_end( iun, "K-POINT"//trim(iotk_index(ik)) )
          !
      ENDDO
      !
      CALL iotk_write_end(iun, "OVERLAPS")
      !
  ENDIF
  !
  ! closing the file
  !
  CALL iotk_close_write(iun)

END SUBROUTINE write_proj_iotk
!
!-----------------------------------------------------------------------
SUBROUTINE write_proj_file ( filproj, proj )
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE lsda_mod,  ONLY : nspin
  USE spin_orb,  ONLY : lspinorb
  USE noncollin_module, ONLY : noncolin
  USE fft_base,  ONLY : dfftp
  USE klist,     ONLY : xk, nkstot
  USE run_info,  ONLY : title
  USE cell_base, ONLY : at, ibrav, celldm
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, atm
  USE wvfct,     ONLY : nbnd
  USE basis,     ONLY : natomwfc
  USE gvect,     ONLY : gcutm 
  USE gvecs,     ONLY : dual
  USE gvecw,     ONLY : ecutwfc
  USE projections, ONLY : nlmchi
  !
  IMPLICIT NONE
  CHARACTER (len=*), INTENT(in) :: filproj
  REAL(DP), INTENT(IN) :: proj(natomwfc,nbnd,nkstot)
  !
  CHARACTER(256) :: filename
  REAL (DP), EXTERNAL :: compute_mj
  INTEGER :: is, ik, nwfc, ibnd, nksinit, nkslast, iunproj=33
  !
  IF ( TRIM(filproj) == ' ' ) RETURN
  !
  DO is=1,nspin
     IF (nspin==2) THEN
        IF (is==1) filename=trim(filproj)//'.projwfc_up'
        IF (is==2) filename=trim(filproj)//'.projwfc_down'
        nksinit=(nkstot/2)*(is-1)+1
        nkslast=(nkstot/2)*is
     ELSE
        filename=trim(filproj)//'.projwfc_up'
        nksinit=1
        nkslast=nkstot
     ENDIF
     CALL write_io_header(filename, iunproj, title, dfftp%nr1x, dfftp%nr2x, &
          dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp, ibrav, &
          celldm, at, gcutm, dual, ecutwfc, nkstot/nspin, nbnd, natomwfc)
     DO nwfc = 1, natomwfc
        IF (lspinorb) THEN
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,2i5,f5.1,1x,f5.1)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l, &
                nlmchi(nwfc)%jj, &
                compute_mj(nlmchi(nwfc)%jj,nlmchi(nwfc)%l, nlmchi(nwfc)%m)
        ELSE IF (noncolin) THEN
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,3i5,1x,f4.1)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l, &
                nlmchi(nwfc)%m, &
                0.5d0-int(nlmchi(nwfc)%ind/(2*nlmchi(nwfc)%l+2))
        ELSE
           WRITE(iunproj,'(2i5,1x,a4,1x,a2,1x,3i5)') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%els, nlmchi(nwfc)%n, nlmchi(nwfc)%l,  &
                nlmchi(nwfc)%m
        END IF
        DO ik=nksinit,nkslast
           DO ibnd=1,nbnd
              WRITE(iunproj,'(2i8,f20.10)') ik,ibnd, abs(proj(nwfc,ibnd,ik))
           ENDDO
        ENDDO
     ENDDO
     CLOSE(iunproj)
  ENDDO
  !
END SUBROUTINE write_proj_file
