!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      module vanderwaals

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

        logical :: tvdw = .false.

        PUBLIC :: vdw, tvdw

      contains


!---------------------------------
      subroutine VdW(evdw, atoms, fion, box) 

      USE constants, ONLY: au => BOHR_RADIUS_ANGS
      USE cell_base, ONLY: s_to_r, boxdimensions, pbcs
      USE mp_global, ONLY: me_image, root_image
      USE atoms_type_module, ONLY: atoms_type

! 
! nat == numero atomi 
! x,y,z == coordinate cartesiane 
! force == forze 
! evdw == energia di VdW
! nsp(1) == numero atomi specie 1 
! csp() == coeffic. di VdW 
!
      implicit none 

      TYPE (atoms_type), intent(in) :: atoms
      type(boxdimensions), intent(in) :: box
      REAL(DP), intent(out) :: evdw
      REAL(DP), intent(out) :: fion(:,:)

      REAL(DP) alp,rcc,rcut,cutoff
      parameter (alp=2.d0,rcc=6.5d0,rcut=3.0d0,cutoff=14.0d0) 

      REAL(DP) csp11, csp12, csp22
      parameter (csp11=1.0191452D0, csp12=0.2239317D0, csp22=0.04364401D0)

      REAL(DP) sij(3),rij(3),sij_image(3) 
      REAL(DP) csp1, dist, ff,dist6,fun,fact,cont
      REAL(DP) force( 3, atoms%nat )
      integer i,j,is,js,ia,ja,ix,iy,iz,iesr
      logical:: tzero,tshift

        force=0.d0
        evdw =0.d0
        iesr=1

        if(atoms%nsp.ne.2 .or. .not.tvdw) then
          return
        endif

        do i=1,atoms%nat

          if(i.le.atoms%na(1)) then
            ia = i
            is = 1
          else
            ia = i - atoms%na(1)
            is = 2 
          end if

          do j=1,atoms%nat

            if(j.le.atoms%na(1)) then
              ja = j
              js = 1
            else
              ja = j - atoms%na(1)
              js = 2 
            end if

            if (i.eq.j) then
              sij=0.d0
              tzero=.true.
            else
              tzero=.false.  
              sij = atoms%taus(:,i) - atoms%taus(:,j)
              CALL PBCS(sij(1),sij(2),sij(3),sij(1),sij(2),sij(3),1)
            end if
 
            do ix=-iesr,iesr
              sij_image(1)= sij(1)+DBLE(ix)
              do iy=-iesr,iesr
                sij_image(2)= sij(2)+DBLE(iy)
                do iz=-iesr,iesr
                  sij_image(3)= sij(3)+DBLE(iz)     
                  tshift=ix.eq.0 .and. iy.eq.0 .and. iz.eq.0
                  if(.not.(tzero.and.tshift)) then
                    call s_to_r(sij_image,rij,box)
                    dist = ( rij(1)**2 + rij(2)**2 + rij(3)**2 )**0.5d0
!
! ...               c-c vdw coefficient
!
                    CSP1 = csp11  
!
! ...               c-h vdw coefficient
!
                    if ( (i.le.atoms%na(1).and.j.gt.atoms%na(1)) .or.  &
                         (i.gt.atoms%na(1).and.j.le.atoms%na(1)) ) then
                      CSP1 = csp12 
                    end if
!
! ...               h-h vdw coefficient
!
                    if (i.gt.atoms%na(1).and.j.gt.atoms%na(1))  then
                      CSP1 = csp22 
                    end if
!
! ...               apply lower boundary cut-off
!
                    if(dist.lt.rcut) then
                      dist = rcut
                    end if

                    ff = alp * (rcc - dist) 
                    dist6 = dist**6
                    fun = - CSP1 / dist6 * cutofun_vdw(ff) / (au)**6

                    if(dist.lt.rcut) then
                      fact = 0.d0 
                    else
                      fact = (6.d0 * CSP1/dist**7 * cutofun_vdw(ff) + &
                      alp * dcutofun_vdw(ff) * CSP1/dist6) / (au)**6
                    endif

                    evdw = evdw + fun
                    force(1,i) = force(1,i)  - fact * rij(1) / dist
                    force(2,i) = force(2,i)  - fact * rij(2) / dist
                    force(3,i) = force(3,i)  - fact * rij(3) / dist
                  endif 
                enddo !iz
              enddo !iy
            enddo !ix
          enddo !j
        enddo !i
        evdw=evdw/2.d0

        IF( me_image == root_image ) THEN
          fion( :, 1:atoms%nat ) = fion( :, 1:atoms%nat ) + force( :, 1:atoms%nat )
        END IF

        return
        end subroutine vdw

!==================================================================

        function cutofun_vdw(xin) 
        implicit none 

        REAL(DP) cutofun_vdw
        REAL(DP), intent(in) :: xin
        REAL(DP) x

        if( xin .gt. 30.d0 ) then
          x = 30.d0
        else
          x = xin
        endif
        cutofun_vdw = 1.d0 / (exp(x) + 1.d0) 

        return
        end function cutofun_vdw
!================================================================== c
!==================================================================
        function dcutofun_vdw(xin) 
        implicit none 

        REAL(DP) dcutofun_vdw
        REAL(DP), intent(in) :: xin
        REAL(DP) x

        if( xin .gt. 30.d0 ) then
          x = 30.d0
        else
          x = xin
        endif
        dcutofun_vdw = - exp(x) / (exp(x) + 1.d0)**2

        return
        end function dcutofun_vdw
!==================================================================




      subroutine baricentro(bar,vectors,nvec)
      implicit none
      integer, intent(in) :: nvec
      REAL(DP), intent(out) :: bar(3)
      REAL(DP), intent(in) :: vectors(3,nvec)
      integer i,j
      do i = 1,3
        bar(i) = 0.0d0
        do j = 1,nvec
          bar(i) = bar(i) + vectors(i,j)
        end do
        bar(i) = bar(i) / DBLE(nvec)
      end do
      return
      end subroutine baricentro

      REAL(DP) function distanza(u,v)
      implicit none
      REAL(DP) u(3),v(3)
      distanza = (u(1)-v(1))**2 + (u(2)-v(2))**2 + (u(3)-v(3))**2
      distanza = sqrt(distanza)
      return
      end function distanza



!      REAL(DP) FUNCTION VDW_FORCES(C6,IESR,FION,STAU0,NA,NAX,NSP)
!
!      USE cell_base, only: R_TO_S, S_TO_R
!
!      implicit none
!
!      REAL(DP) c6
!      integer iesr
!      REAL(DP) fion(3,nax,nsp)
!      REAL(DP) stau0(3,nax,nsp)
!      integer na(nsp)
!      integer nax,nsp
!
!      REAL(DP)  EVDW
!      REAL(DP)  distanza
!      integer i,j,k,ix,iy,iz,infm,m,l,ishft,im
!      REAL(DP) XLM, YLM, ZLM, ZERO
!      REAL(DP) sxlm(3),rxlm(3),ERRE2,RLM,ADDEVDW,ADDPRE
!      REAL(DP) FXX, REPAND
!      REAL(DP) molbar(3,NAX)
!      REAL(DP) molecola(3,NAX),tau(3),rdis
!      REAL(DP) fmol(3,NAX)
!      REAL(DP) bond_len_au
!      integer iatmol(NAX,NSP),imol,nmol,natmol
!      logical TZERO
!
!
!      bond_len_au = 2.0d0
!      imol = 1
!      do i=1,na(1)
!        im = 1
!        call S_TO_R(stau0(1,i,1),molecola(1,im))
!        iatmol(i,1) = im
!        im = im + 1
!        do j = 1,na(2)
!          call S_TO_R(stau0(1,j,2),tau)
!          rdis = distanza(molecola(1,1),tau)
!          if(rdis.lt.bond_len_au) then
!            call S_TO_R(stau0(1,j,2),molecola(1,im))
!            iatmol(j,2) = im
!            im = im + 1
!          end if
!        end do
!        natmol = im - 1
!        call baricentro(tau,molecola,natmol)
!        call r_to_s(tau,molbar(1,imol))
!        imol = imol + 1 
!      end do         
!      nmol = imol - 1
!          
!
!
!      EVDW  = 0.D0
!
!      call azzera(fmol,3*nax)
!      DO L=1,nmol
!        DO M= L,nmol
!          IF(L.EQ.M) THEN
!            XLM=0.D0
!            YLM=0.D0
!            ZLM=0.D0
!            TZERO=.TRUE.
!          ELSE
!            TZERO=.FALSE.
!            XLM = molbar(1,l) - molbar(1,m)
!            YLM = molbar(2,l) - molbar(2,m)
!            ZLM = molbar(3,l) - molbar(3,m) 
!            CALL PBCS(XLM,YLM,ZLM,XLM,YLM,ZLM,1)
!          END IF
!          DO IX=-IESR,IESR
!            DO IY=-IESR,IESR
!              DO IZ=-IESR,IESR
!                ISHFT=IX*IX+IY*IY+IZ*IZ
!                IF(.NOT.(TZERO.AND.ISHFT.EQ.0)) THEN
!                  sxlm(1) = XLM + DBLE(IX)
!                  sxlm(2) = YLM + DBLE(IY)
!                  sxlm(3) = ZLM + DBLE(IZ)
!                  CALL S_TO_R(sxlm,rxlm)
!                  ERRE2 = rxlm(1)**2 + rxlm(2)**2 + rxlm(3)**2
!                  RLM   = SQRT(ERRE2)
!                  IF (TZERO) THEN
!                    ZERO=0.5D0
!                  ELSE
!                    ZERO=1.D0
!                  END IF
!                  ADDEVDW = - C6 / RLM**6
!                  EVDW = EVDW + ZERO*ADDEVDW
!                  ADDPRE = - 6.0D0 * C6 /RLM**8
!                  REPAND = ZERO*(ADDEVDW + ADDPRE)
!                  DO I=1,3
!                    FXX = REPAND*rxlm(I)
!                    FMOL(I,L) = FMOL(I,L) + FXX
!                    FMOL(I,M) = FMOL(I,M) - FXX
!                  END DO
!                END IF
!              END DO    ! IZ
!            END DO      ! IY
!          END DO        ! IX
!        END DO          ! M
!      END DO            ! L
!
!      do i=1,nsp
!        do j=1,na(i)
!          do k=1,3
!            fion(k,j,i)=fion(k,j,i)+fmol(k,iatmol(j,i))/REAL(natmol)
!          end do
!        end do
!      end do
!
!      VDW_FORCES = EVDW
!      return
!      end FUNCTION VDW_FORCES
!
!
!      subroutine VDW_STRESS(C6,IESR,STAU0,DVDW,NA,NAX,NSP)
!
!      USE cell_base, only: R_TO_S, S_TO_R
!
!      implicit none
!
!      REAL(DP) c6
!      integer iesr
!      REAL(DP) stau0(3,nax,nsp)
!      REAL(DP) dvdw(6)
!      integer na(nsp)
!      integer nax,nsp
!
!      REAL(DP)  distanza
!      integer i,j,k,ix,iy,iz,infm,m,l,ishft,im
!      REAL(DP) XLM, YLM, ZLM, ZERO
!      REAL(DP) sxlm(3),rxlm(3),ERRE2,RLM,ADDEVDW,ADDPRE
!      REAL(DP) FXX, REPAND
!      REAL(DP) molbar(3,NAX)
!      REAL(DP) molecola(3,NAX),tau(3),rdis
!      REAL(DP) bond_len_au
!      integer iatmol(NAX,NSP),imol,nmol,natmol
!      integer alpha(6),beta(6)
!      logical TZERO
!
!      ALPHA(1) = 1
!      ALPHA(2) = 2 
!      ALPHA(3) = 3
!      ALPHA(4) = 2 
!      ALPHA(5) = 3 
!      ALPHA(6) = 3
!      BETA(1)  = 1 
!      BETA(2)  = 1 
!      BETA(3)  = 1
!      BETA(4)  = 2 
!      BETA(5)  = 2
!      BETA(6)  = 3
!
!      do i=1,6
!        dvdw(i) = 0.0d0
!      end do
!
!      bond_len_au = 2.0d0
!      imol = 1
!      do i=1,na(1)
!        im = 1
!        call S_TO_R(stau0(1,i,1),molecola(1,im))
!        iatmol(i,1) = im
!        im = im + 1
!        do j = 1,na(2)
!          call S_TO_R(stau0(1,j,2),tau)
!          rdis = distanza(molecola(1,1),tau)
!          if(rdis.lt.bond_len_au) then
!            call S_TO_R(stau0(1,j,2),molecola(1,im))
!            iatmol(j,2) = im
!            im = im + 1
!          end if
!        end do
!        natmol = im - 1
!        call baricentro(tau,molecola,natmol)
!        call r_to_s(tau,molbar(1,imol))
!        imol = imol + 1 
!      end do         
!      nmol = imol - 1
!          
!
!      DO L=1,nmol
!        DO M= L,nmol
!          IF(L.EQ.M) THEN
!            XLM=0.D0
!            YLM=0.D0
!            ZLM=0.D0
!            TZERO=.TRUE.
!          ELSE
!            TZERO=.FALSE.
!            XLM = molbar(1,l) - molbar(1,m)
!            YLM = molbar(2,l) - molbar(2,m)
!            ZLM = molbar(3,l) - molbar(3,m) 
!            CALL PBCS(XLM,YLM,ZLM,XLM,YLM,ZLM,1)
!          END IF
!          DO IX=-IESR,IESR
!            DO IY=-IESR,IESR
!              DO IZ=-IESR,IESR
!                ISHFT=IX*IX+IY*IY+IZ*IZ
!                IF(.NOT.(TZERO.AND.ISHFT.EQ.0)) THEN
!                  sxlm(1) = XLM + DBLE(IX)
!                  sxlm(2) = YLM + DBLE(IY)
!                  sxlm(3) = ZLM + DBLE(IZ)
!                  CALL S_TO_R(sxlm,rxlm)
!                  ERRE2 = rxlm(1)**2 + rxlm(2)**2 + rxlm(3)**2
!                  RLM   = SQRT(ERRE2)
!                  IF (TZERO) THEN
!                    ZERO=0.5D0
!                  ELSE
!                    ZERO=1.D0
!                  END IF
!                  ADDPRE = - 6.0D0 * C6 /RLM**8
!                  REPAND = ZERO * ADDPRE
!                  DO I=1,6
!                    FXX = REPAND*rxlm(ALPHA(I))*rxlm(BETA(I))
!                    DVDW(I) = DVDW(I) - FXX
!                  END DO
!                END IF
!              END DO    ! IZ
!            END DO      ! IY
!          END DO        ! IX
!        END DO          ! M
!      END DO            ! L
!
!      return
!      end SUBROUTINE VDW_STRESS

      end module vanderwaals
