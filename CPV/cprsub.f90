!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine formf( tfirst, eself )
  !-----------------------------------------------------------------------

  !computes (a) the self-energy eself of the ionic pseudocharges;
  !         (b) the form factors of: (i) pseudopotential (vps),
  !             (ii) ionic pseudocharge (rhops)
  !         all quantities are returned in common /pseu/
  !         also calculated the derivative of vps with respect to
  !         g^2 (dvps)
  ! 
  USE kinds,           ONLY : DP
  use mp,              ONLY : mp_sum
  use control_flags,   ONLY : iprint, tpre, iprsta
  use io_global,       ONLY : stdout
  use mp_global,       ONLY : intra_image_comm
  use bhs,             ONLY : rc1, rc2, wrc2, wrc1, rcl, al, bl, lloc
  use gvecs,           ONLY : ngs
  use cell_base,       ONLY : omega, tpiba2, tpiba
  use ions_base,       ONLY : rcmax, zv, nsp, na
  use local_pseudo,    ONLY : vps, rhops, dvps, drhops
  use atom,            ONLY : r, rab, mesh, numeric
  use uspp_param,      ONLY : vloc_at, oldvan
  use pseudo_base,     ONLY : compute_rhops, formfn, formfa, compute_eself
  use pseudopotential, ONLY : tpstab, build_pstab, vps_sp, dvps_sp
  use splines,         ONLY : spline
  use reciprocal_vectors, ONLY : gstart, g
  !
  implicit none
  logical      :: tfirst
  real(DP)    :: eself
  !
  real(DP)    :: vpsum, rhopsum
  integer      :: is, ig
  REAL(DP)    :: cost1, xg

  call start_clock( 'formf' )
  !
  ! calculation of gaussian selfinteraction
  !
  eself = compute_eself( na, zv, rcmax, nsp )

  if( tfirst .or. ( iprsta >= 4 ) )then
     WRITE( stdout, 1200 ) eself
  endif
  !
  1200 format(/,3x,'formf: eself=',f10.5)
  !
  IF( tpstab ) THEN
     !
     CALL build_pstab( )
     !
  END IF
  !
  do is = 1, nsp

     IF( tpstab ) THEN
        !
        !  Use interpolation table, with cubic spline
        !
        cost1 = 1.0d0/omega
        !
        IF( gstart == 2 ) THEN
           vps (1,is) =  vps_sp(is)%y(1) * cost1
           dvps(1,is) = dvps_sp(is)%y(1) * cost1
        END IF
        !
        DO ig = gstart, ngs
           xg = SQRT( g(ig) ) * tpiba
           vps (ig,is) = spline(  vps_sp(is), xg ) * cost1
           dvps(ig,is) = spline( dvps_sp(is), xg ) * cost1
        END DO
        !
     ELSE

        if ( numeric(is) ) then

           call formfn( vps(:,is), dvps(:,is), r(:,is), rab(:,is), vloc_at(:,is), &
                        zv(is), rcmax(is), g, omega, tpiba2, mesh(is), &
                        ngs, oldvan(is), tpre )

        else

           !     bhs pseudopotentials
           !
           call formfa( vps(:,is), dvps(:,is), rc1(is), rc2(is), wrc1(is), wrc2(is), &
                        rcl(:,is,lloc(is)), al(:,is,lloc(is)), bl(:,is,lloc(is)),    &
                        zv(is), rcmax(is), g, omega, tpiba2, ngs, gstart, tpre )

        end if

     END IF
     !
     !     fourier transform of local pp and gaussian nuclear charge
     !
     call compute_rhops( rhops(:,is), drhops(:,is), zv(is), rcmax(is), g,   &
                         omega, tpiba2, ngs, tpre )

     if( tfirst .or. ( iprsta >= 4 ) )then
        vpsum = SUM( vps( 1:ngs, is ) )
        rhopsum = SUM( rhops( 1:ngs, is ) )
        call mp_sum( vpsum, intra_image_comm )
        call mp_sum( rhopsum, intra_image_comm )
        WRITE( stdout,1250) vps(1,is),rhops(1,is)
        WRITE( stdout,1300) vpsum,rhopsum
     endif
     !
  end do
  !
  call stop_clock( 'formf' )
  !
  1250 format(3x,'formf:     vps(g=0)=',f12.7,'     rhops(g=0)=',f12.7)
  1300 format(3x,'formf: sum_g vps(g)=',f12.7,' sum_g rhops(g)=',f12.7)
  !
  return
end subroutine formf
!
!-----------------------------------------------------------------------
SUBROUTINE newnlinit()
  !-----------------------------------------------------------------------
  !
  ! ... this routine calculates arrays beta, qradb, qq, qgb, rhocb
  ! ... and derivatives w.r.t. cell parameters dbeta, dqrad 
  ! ... See also comments in nlinit
  !
  use control_flags,    ONLY : tpre
  use pseudopotential,  ONLY : interpolate_beta, interpolate_qradb
  use pseudopotential,  ONLY : exact_beta, tpstab, check_tables
  USE core,             ONLY : core_charge_ftr
  !
  IMPLICIT NONE
  !
  LOGICAL :: recompute_table
  ! 
  ! ... initialization for vanderbilt species
  !
  recompute_table = tpre .AND. check_tables()
  !
  IF ( recompute_table ) &
     CALL errore( ' newnlinit', &
                  'interpolation tables recalculation, not implemented yet', 1 )
  !
  CALL interpolate_qradb( tpre )
  !
  !
  !     initialization that is common to all species
  !
  IF( tpstab ) THEN
     !
     CALL interpolate_beta( tpre )
     !
  ELSE
     !
     ! ... this is mainly for testing
     !
     CALL exact_beta( tpre )
     !
  END IF
  !
  ! ... non-linear core-correction   ( rhocb(ig,is) )
  !
  CALL core_charge_ftr( tpre )
  !
  RETURN
  !
END SUBROUTINE newnlinit
!
!-----------------------------------------------------------------------
subroutine nlfh( bec, dbec, lambda )
  !-----------------------------------------------------------------------
  !
  !     contribution to the internal stress tensor due to the constraints
  !
  USE kinds,          ONLY : DP
  use cvan,           ONLY : nvb, ish
  use uspp,           ONLY : nhsa => nkb, qq
  use uspp_param,     ONLY : nh, nhm
  use ions_base,      ONLY : na
  use electrons_base, ONLY : nbspx, nbsp, nudx, nspin, nupdwn, iupdwn
  use cell_base,      ONLY : omega, h
  use constants,      ONLY : pi, fpi
  use stre,           ONLY : stress
!
  implicit none

  real(DP), intent(in) ::  bec( nhsa, nbsp ), dbec( nhsa, nbsp, 3, 3 )
  real(DP), intent(in) ::  lambda( nudx, nudx, nspin )
!
  integer   :: i, j, ii, jj, inl, iv, jv, ia, is, iss, nss, istart
  real(DP) :: fpre(3,3)
  !
  REAL(DP), ALLOCATABLE :: tmpbec(:,:), tmpdh(:,:), temp(:,:)
  !
  ALLOCATE ( tmpbec(nhm,nudx), tmpdh(nudx,nhm), temp(nudx,nudx) )
  !
      fpre(:,:) = 0.d0
      do ii=1,3
         do jj=1,3
            do is=1,nvb
               do ia=1,na(is)
!
                  do iss = 1, nspin
                     !
                     istart = iupdwn( iss )
                     nss    = nupdwn( iss )
                     !
                     tmpbec = 0.d0
                     tmpdh  = 0.d0
!
                     do iv=1,nh(is)
                        do jv=1,nh(is)
                           inl=ish(is)+(jv-1)*na(is)+ia
                           if(abs(qq(iv,jv,is)).gt.1.e-5) then
                              do i = 1, nss
                                 tmpbec(iv,i) = tmpbec(iv,i) +             &
     &                              qq(iv,jv,is) * bec(inl, i + istart - 1)
                              end do
                           endif
                        end do
                     end do
!
                     do iv=1,nh(is)
                        inl=ish(is)+(iv-1)*na(is)+ia
                        do i = 1, nss
                           tmpdh(i,iv) = dbec( inl, i+istart-1, ii, jj )
                        end do
                     end do
!
                     if(nh(is).gt.0)then

                        temp = 0.d0
!
                        call MXMA                                          &
     &                    (tmpdh,1,nudx,tmpbec,1,nhm,temp,1,nudx,nss,nh(is),nss)
!
                        ! do j=1,nss
                        !    do i=1,nss
                        !       temp(i,j)=temp(i,j)*lambda(i,j,iss)
                        !    end do
                        ! end do
                        ! fpre(ii,jj)=fpre(ii,jj)+2.*SUM(temp(1:nss,1:nss))
                        !
                        DO j = 1, nss
                           DO i = 1, nss
                              fpre(ii,jj) = fpre(ii,jj) + 2*temp(i,j)*lambda(i,j,iss)
                           END DO
                        END DO

                     endif
                     !
                  end do
                  !
               end do
            end do
         end do
      end do
      do i=1,3
         do j=1,3
            stress(i,j)=stress(i,j)+(fpre(i,1)*h(j,1)+                  &
     &           fpre(i,2)*h(j,2)+fpre(i,3)*h(j,3))/omega
         enddo
      enddo
!
  DEALLOCATE ( tmpbec, tmpdh, temp )

  return
end subroutine nlfh


!-----------------------------------------------------------------------
subroutine nlinit
  !-----------------------------------------------------------------------
  !
  !     this routine allocates and initalizes arrays beta, qradb, qq, qgb,
  !     rhocb, and derivatives w.r.t. cell parameters dbeta, dqrad 
  !
  !       beta(ig,l,is) = 4pi/sqrt(omega) y^r(l,q^)
  !                               int_0^inf dr r^2 j_l(qr) betar(l,is,r)
  !
  !       Note that beta(g)_lm,is = (-i)^l*beta(ig,l,is) (?)
  !
  !       qradb(ig,l,k,is) = 4pi/omega int_0^r dr r^2 j_l(qr) q(r,l,k,is)
  !
  !       qq_ij=int_0^r q_ij(r)=omega*qg(g=0)
  !
  !     beta and qradb are first calculated on a fixed linear grid in |G|
  !     (betax, qradx) then calculated on the box grid by interpolation
  !     (this is done in routine newnlinit)
  !     
      use parameters,      ONLY : lmaxx
      use control_flags,   ONLY : iprint, tpre, program_name
      use io_global,       ONLY : stdout, ionode
      use gvecw,           ONLY : ngw
      use cvan,            ONLY : ish, nvb
      use core,            ONLY : rhocb, nlcc_any, allocate_core
      use constants,       ONLY : pi, fpi
      use ions_base,       ONLY : na, nsp
      use uspp,            ONLY : aainit, beta, qq, dvan, nhtol, nhtolm, indv,&
                                  nhsa => nkb, nhsavb=>nkbus
      use uspp_param,      ONLY : kkbeta, qqq, nqlc, betar, lmaxq, dion,&
                                  nbeta, nbetam, lmaxkb, lll, nhm, nh, tvanp
      use atom,            ONLY : mesh, r, rab, nlcc, numeric
      use qradb_mod,       ONLY : qradb
      use qgb_mod,         ONLY : qgb
      use gvecb,           ONLY : ngb
      use gvecp,           ONLY : ngm
      use cdvan,           ONLY : dbeta
      use dqrad_mod,       ONLY : dqrad
      use dqgb_mod,        ONLY : dqgb
      use betax,           ONLY : qradx, dqradx, refg, betagx, mmx, dbetagx
      use pseudopotential, ONLY : pseudopotential_indexes, compute_dvan, &
                                  compute_betagx, compute_qradx
      USE grid_dimensions, ONLY : nnrx

!
      implicit none
!
      integer  is, il, l, ir, iv, jv, lm, ind, ltmp, i0
      real(8), allocatable:: fint(:), jl(:),  jltmp(:), djl(:),    &
     &              dfint(:)
      real(8) xg, xrg, fac

      IF( ionode ) THEN
        WRITE( stdout, 100 )
 100    FORMAT( //, &
                3X,'Pseudopotentials initialization',/, &
                3X,'-------------------------------' )
      END IF

      !
      !   initialize indexes
      !
      CALL pseudopotential_indexes( nlcc_any )
      !
      !   initialize array ap
      !
      call aainit( lmaxkb + 1 )
      !
      CALL allocate_core( nnrx, ngm, ngb, nsp )
      !
      !
      allocate( beta( ngw, nhm, nsp ) )
      allocate( qradb( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp ) )
      allocate( qgb( ngb, nhm*(nhm+1)/2, nsp ) )
      allocate( qq( nhm, nhm, nsp ) )
      qradb(:,:,:,:) = 0.d0
      qq  (:,:,:) =0.d0
      IF (tpre) THEN
         allocate( dqrad( ngb, nbetam*(nbetam+1)/2, lmaxq, nsp, 3, 3 ) )
         allocate( dqgb( ngb, nhm*(nhm+1)/2, nsp, 3, 3 ) )
         allocate( dbeta( ngw, nhm, nsp, 3, 3 ) )
         dqrad(:,:,:,:,:,:) = 0.d0
      END IF
      !
      !     initialization for vanderbilt species
      !
      CALL compute_qradx( tpre )
      !    
      !     initialization that is common to all species
      !   
      WRITE( stdout, fmt="(//,3X,'Common initialization' )" )

      do is = 1, nsp
         WRITE( stdout, fmt="(/,3X,'Specie: ',I5)" ) is
         if ( .not. numeric(is) ) then
            fac=1.0
         else
            !     fac converts ry to hartree
            fac=0.5
         end if
         do iv = 1, nh(is)
            WRITE( stdout,901) iv, indv(iv,is), nhtol(iv,is)
         end do
 901     format(2x,i2,'  indv= ',i2,'   ang. mom= ',i2)
         !
         WRITE( stdout,*)
         WRITE( stdout,'(20x,a)') '    dion '
         do iv = 1, nbeta(is)
            WRITE( stdout,'(8f9.4)') ( fac*dion(iv,jv,is), jv = 1, nbeta(is) )
         end do
         !
      end do
      !
      !   calculation of array  betagx(ig,iv,is)
      !
      call compute_betagx( tpre )
      !
      !   calculate array  dvan(iv,jv,is)
      !
      call compute_dvan()
      !
      ! newnlinit stores qgb and qq, calculates arrays  beta  qradb  rhocb
      ! and derivatives wrt cell    dbeta dqrad
      !
      call newnlinit()

      return
end subroutine nlinit

!-------------------------------------------------------------------------
subroutine qvan2b(ngy,iv,jv,is,ylm,qg)
  !--------------------------------------------------------------------------
  !
  !     q(g,l,k) = sum_lm (-i)^l ap(lm,l,k) yr_lm(g^) qrad(g,l,l,k)
  !
  USE kinds,         ONLY : DP
  use control_flags, ONLY : iprint, tpre
  use qradb_mod,     ONLY : qradb
  use uspp,          ONLY : nlx, lpx, lpl, ap, indv, nhtolm
  use gvecb,         ONLY : ngb
  use uspp_param,    ONLY : lmaxq
! 
  implicit none
  !
  integer,      intent(in)  :: ngy, iv, jv, is
  real(DP),    intent(in)  :: ylm( ngb, lmaxq*lmaxq )
  complex(DP), intent(out) :: qg( ngb )
!
  integer      :: ivs, jvs, ijvs, ivl, jvl, i, ii, ij, l, lp, ig
  complex(DP) :: sig
  ! 
  !       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
  !       ivs = 1..4     s_1 s_2 p_1 p_2
  !       ivl = 1..4     s p_x p_z p_y
  ! 
  ivs=indv(iv,is)
  jvs=indv(jv,is)
  if (ivs >= jvs) then
     ijvs = ivs*(ivs-1)/2 + jvs
  else
     ijvs = jvs*(jvs-1)/2 + ivs
  end if
  ! ijvs is the packed index for (ivs,jvs)
  ivl=nhtolm(iv,is)
  jvl=nhtolm(jv,is)
  if (ivl > nlx .OR. jvl > nlx) &
       call errore (' qvan2b ', ' wrong dimensions', MAX(ivl,jvl))
  !
  qg(:) = (0.d0, 0.d0)
  !
  !     lpx = max number of allowed y_lm
  !     lp  = composite lm to indentify them
  !
  do i=1,lpx(ivl,jvl)
     lp=lpl(ivl,jvl,i)
     if (lp > lmaxq*lmaxq) call errore(' qvan2b ',' lp out of bounds ',lp)
     !
     !     extraction of angular momentum l from lp:  
     !     l = int ( sqrt( DBLE(l-1) + epsilon) ) + 1
     !
     if (lp == 1) then
        l=1         
     else if ((lp >= 2) .and. (lp <= 4)) then
        l=2
     else if ((lp >= 5) .and. (lp <= 9)) then
        l=3
     else if ((lp >= 10).and.(lp <= 16)) then
        l=4
     else if ((lp >= 17).and.(lp <= 25)) then
        l=5
     else if ((lp >= 26).and.(lp <= 36)) then 
        l=6
     else if ((lp >= 37).and.(lp <= 49)) then 
        l=7
     else
        call errore(' qvan2b ',' not implemented ',lp)
     endif
     !     
     !       sig= (-i)^l
     !
     sig=(0.,-1.)**(l-1)
     sig=sig*ap(lp,ivl,jvl)
     do ig=1,ngy
        qg(ig)=qg(ig)+sig*ylm(ig,lp)*qradb(ig,ijvs,l,is)
     end do
  end do

  return
end subroutine qvan2b

!-------------------------------------------------------------------------
subroutine dqvan2b(ngy,iv,jv,is,ylm,dylm,dqg)
  !--------------------------------------------------------------------------
  !
  !     dq(i,j) derivatives wrt to h(i,j) of q(g,l,k) calculated in qvan2b
  !
  USE kinds,         ONLY : DP
  use control_flags, ONLY : iprint, tpre
  use qradb_mod,     ONLY : qradb
  use uspp,          ONLY : nlx, lpx, lpl, ap, indv, nhtolm
  use gvecb,         ONLY : ngb
  use dqrad_mod,     ONLY : dqrad
  use uspp_param,    ONLY : lmaxq

  implicit none

  integer,      intent(in)  :: ngy, iv, jv, is
  REAL(DP),    INTENT(IN)  :: ylm( ngb, lmaxq*lmaxq ), dylm( ngb, lmaxq*lmaxq, 3, 3 )
  complex(DP), intent(out) :: dqg( ngb, 3, 3 )

  integer      :: ivs, jvs, ijvs, ivl, jvl, i, ii, ij, l, lp, ig
  complex(DP) :: sig
  !
  ! 
  !       iv  = 1..8     s_1 p_x1 p_z1 p_y1 s_2 p_x2 p_z2 p_y2
  !       ivs = 1..4     s_1 s_2 p_1 p_2
  !       ivl = 1..4     s p_x p_z p_y
  ! 

  ivs=indv(iv,is)
  jvs=indv(jv,is)
  if (ivs >= jvs) then
     ijvs = ivs*(ivs-1)/2 + jvs
  else
     ijvs = jvs*(jvs-1)/2 + ivs
  end if
  ! ijvs is the packed index for (ivs,jvs)
  ivl=nhtolm(iv,is)
  jvl=nhtolm(jv,is)
  if (ivl > nlx .OR. jvl > nlx) &
       call errore (' qvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl))
  !
  dqg(:,:,:) = (0.d0, 0.d0)

  !  lpx = max number of allowed y_lm
  !  lp  = composite lm to indentify them

  do i=1,lpx(ivl,jvl)
     lp=lpl(ivl,jvl,i)
     if (lp > lmaxq*lmaxq) call errore(' dqvan2b ',' lp out of bounds ',lp)

     !  extraction of angular momentum l from lp:  
     !  l = int ( sqrt( DBLE(l-1) + epsilon) ) + 1
     !
     if (lp == 1) then
        l=1         
     else if ((lp >= 2) .and. (lp <= 4)) then
        l=2
     else if ((lp >= 5) .and. (lp <= 9)) then
        l=3
     else if ((lp >= 10).and.(lp <= 16)) then
        l=4
     else if ((lp >= 17).and.(lp <= 25)) then
        l=5
     else if ((lp >= 26).and.(lp <= 36)) then 
        l=6
     else if ((lp >= 37).and.(lp <= 49)) then 
        l=7
     else
        call errore(' qvan2b ',' not implemented ',lp)
     endif
     !     
     !       sig= (-i)^l
     !
     sig=(0.,-1.)**(l-1)
     sig=sig*ap(lp,ivl,jvl)
     do ij=1,3
        do ii=1,3
           do ig=1,ngy
              dqg(ig,ii,ij) = dqg(ig,ii,ij) +  sig *                &
 &                    ( ylm(ig,lp) * dqrad(ig,ijvs,l,is,ii,ij) -    &
 &                     dylm(ig,lp,ii,ij)*qradb(ig,ijvs,l,is)   ) ! SEGNO
           end do
        end do
     end do
  end do
  !
  return
end subroutine dqvan2b

!-----------------------------------------------------------------------
subroutine dylmr2_( nylm, ngy, g, gg, ainv, dylm )
  !-----------------------------------------------------------------------
  !
  ! temporary CP interface for PW routine dylmr2
  ! dylmr2  calculates d Y_{lm} /d G_ipol
  ! dylmr2_ calculates G_ipol \sum_k h^(-1)(jpol,k) (dY_{lm} /dG_k)
  !
  USE kinds, ONLY: DP

  implicit none
  !
  integer,   intent(IN)  :: nylm, ngy
  real(DP), intent(IN)  :: g (3, ngy), gg (ngy), ainv(3,3)
  real(DP), intent(OUT) :: dylm (ngy, nylm, 3, 3)
  !
  integer :: ipol, jpol, lm, ig
  real(DP), allocatable :: dylmaux (:,:,:)
  !
  allocate ( dylmaux(ngy,nylm,3) )
  !
  dylmaux(:,:,:) = 0.d0
  !
  do ipol =1,3
     call dylmr2 (nylm, ngy, g, gg, dylmaux(1,1,ipol), ipol)
  enddo
  !
  do ipol =1,3
     do jpol =1,3
        do lm=1,nylm
           do ig = 1, ngy
              dylm (ig,lm,ipol,jpol) = (dylmaux(ig,lm,1) * ainv(jpol,1) + &
                                        dylmaux(ig,lm,2) * ainv(jpol,2) + &
                                        dylmaux(ig,lm,3) * ainv(jpol,3) ) &
                                       * g(ipol,ig)
           end do
        end do
     end do
  end do
  !
  deallocate ( dylmaux )
  !
  return
  !
end subroutine dylmr2_
