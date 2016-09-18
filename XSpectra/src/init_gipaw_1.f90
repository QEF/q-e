!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_gipaw_1
  !----------------------------------------------------------------------
  !
  ! This routine initialize the variables of the paw projector
  ! and create the projectors in radial part (paw_betar)
  !
  USE kinds ,      ONLY : dp
  USE parameters , ONLY : lqmax , lmaxx
  USE gipaw_module,ONLY : nbrx
  USE cell_base ,  ONLY : omega
  USE ions_base,   ONLY : nat, ntyp => nsp, ityp
  USE constants,   ONLY : fpi
  USE us,          ONLY : dq, nqx, tab, tab_d2y, qrad, spline_ps
  USE paw_gipaw,   ONLY : paw_recon, paw_nkb, paw_lmaxkb
  USE splinelib
  USE uspp,        ONLY : ap, aainit
  USE atom,        ONLY : rgrid, msh
  USE io_global,   ONLY : stdout
  USE mp_global,   ONLY : intra_pool_comm
  USE mp,          ONLY : mp_sum
  USE matrix_inversion
  !
  implicit none
  !
  !     here a few local variables
  !

  integer :: nt, ih, jh, nb, l, m, ir, iq, startq
  INTEGER :: lastq, na, j, n1, n2, ndm, nrs, nrc, lmaxkb
  ! various counters
  real(DP), allocatable :: aux (:), aux1 (:), besr (:)
  ! various work space
  real(DP) :: prefr, pref, qi, norm
  ! the prefactor of the q functions
  ! the prefactor of the beta functions
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  real(DP), allocatable :: s(:,:), sinv(:,:)
  ! the spherical harmonics
  real(DP) ::  vqint
  ! interpolated value

  real(DP) rc,rs,pow

  integer :: n_overlap_warnings

  real(DP), allocatable :: xdata(:)
  real(DP) :: d1, scaling_factor

  call start_clock ('init_gipaw_1')
  !
  !    Initialization of the variables
  !
  ndm = MAXVAL ( msh(1:ntyp) )
  allocate ( aux(ndm) )
  allocate ( aux1(ndm) )
  allocate ( besr(ndm) )

  paw_lmaxkb = 0
  do nt = 1, ntyp
     lmaxkb = 0
     paw_recon(nt)%paw_nh = 0
     do nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        paw_recon(nt)%paw_nh = paw_recon(nt)%paw_nh + 2 * l + 1
        lmaxkb = max ( lmaxkb, l )
     end do
     paw_lmaxkb = max ( paw_lmaxkb, lmaxkb )

     allocate ( paw_recon(nt)%paw_nhtol(paw_recon(nt)%paw_nh) )
     allocate ( paw_recon(nt)%paw_nhtom(paw_recon(nt)%paw_nh) )
     allocate ( paw_recon(nt)%paw_indv(paw_recon(nt)%paw_nh) )
     allocate ( paw_recon(nt)%paw_tab(nqx,nbrx) )
     allocate ( paw_recon(nt)%paw_nl(0:lmaxkb) )
     allocate ( paw_recon(nt)%paw_iltonh(0:lmaxkb,paw_recon(nt)%paw_nh) )
  END DO

  ! calculate the number of beta functions of the solid
  !
  paw_nkb = 0
  do na = 1, nat
     nt = ityp(na)
     paw_nkb = paw_nkb + paw_recon(nt)%paw_nh
  end do

  n_overlap_warnings = 0

  prefr = fpi / omega
  !
  !   For each pseudopotential we initialize the indices nhtol, nhtom,
  !   indv,
  !
  do nt = 1, ntyp
     paw_recon(nt)%paw_nl = 0
     paw_recon(nt)%paw_iltonh = 0
     ih = 1
     do nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        paw_recon(nt)%paw_nl(l) = paw_recon(nt)%paw_nl(l) + 1
        paw_recon(nt)%paw_iltonh(l,paw_recon(nt)%paw_nl(l)) = nb
        do m = 1, 2 * l + 1
           paw_recon(nt)%paw_nhtol(ih) = l
           paw_recon(nt)%paw_nhtom(ih) = m
           paw_recon(nt)%paw_indv(ih) = nb
           ih = ih + 1
        end do
     end do

     !
     ! Rescale the wavefunctions so that int_0^rc f|psi|^2=1
     !

     pow = 1.0_dp
     do j = 1, paw_recon(nt)%paw_nbeta
        rc = paw_recon(nt)%psphi(j)%label%rc
        !rs = 1.0_dp / 3.0_dp * rc
        rs = 2.0_dp / 3.0_dp * rc
        nrc = COUNT ( rgrid(nt)%r(1:msh(nt)) <= rc )
        nrs = COUNT ( rgrid(nt)%r(1:msh(nt)) <= rs )
        !<debug>
!        write(stdout,*) "ZZZ: ", rc, rs, nrc, nrs
        !</debug>
        IF ( nrc < 1 .OR. nrc > msh(nt) ) &
             CALL errore ( "init_gipaw_1", "impossible value for nrc", 1 )
        IF ( nrs < 1 .OR. nrs > msh(nt) ) &
             CALL errore ( "init_gipaw_1", "impossible value for nrs", 1 )
        paw_recon(nt)%psphi(j)%label%nrc = nrc
        paw_recon(nt)%aephi(j)%label%nrc = nrc
        paw_recon(nt)%psphi(j)%label%nrs = nrs
        paw_recon(nt)%aephi(j)%label%nrs = nrs

        !<apsi> Scaling removed
        !scaling_factor = paw_recon(nt)%aephi(j)%psi(nrc) &
        !     / paw_recon(nt)%psphi(j)%psi(nrc)
        !paw_recon(nt)%psphi(j)%psi(:) = paw_recon(nt)%psphi(j)%psi(:) &
        !     * scaling_factor
        !</apsi>

        call step_f ( aux, paw_recon(nt)%psphi(j)%psi**2, rgrid(nt)%r(:), &
             nrs, nrc, pow, msh(nt) )
        call simpson ( msh(nt), aux, rgrid(nt)%rab, norm )

        paw_recon(nt)%psphi(j)%psi = paw_recon(nt)%psphi(j)%psi / sqrt(norm)
        paw_recon(nt)%aephi(j)%psi = paw_recon(nt)%aephi(j)%psi / sqrt(norm)

     end do

     !
     !   calculate the overlap matrix
     !

     aux = 0.0_dp
     do l = 0, ubound(paw_recon(nt)%paw_nl,1)
        if ( paw_recon(nt)%paw_nl(l) > 0 ) then
           allocate ( s(paw_recon(nt)%paw_nl(l),paw_recon(nt)%paw_nl(l)) )
           allocate ( sinv(paw_recon(nt)%paw_nl(l),paw_recon(nt)%paw_nl(l)) )
           do ih = 1, paw_recon(nt)%paw_nl(l)
              n1 = paw_recon(nt)%paw_iltonh(l,ih)
              do jh = 1, paw_recon(nt)%paw_nl(l)
                 n2 = paw_recon(nt)%paw_iltonh(l,jh)

                 nrc = MIN ( paw_recon(nt)%psphi(n1)%label%nrc, &
                      paw_recon(nt)%psphi(n2)%label%nrc )
                 nrs = MIN ( paw_recon(nt)%psphi(n1)%label%nrs, &
                      paw_recon(nt)%psphi(n2)%label%nrs )

                 call step_f ( aux, paw_recon(nt)%psphi(n1)%psi(1:msh(nt)) &
                      * paw_recon(nt)%psphi(n2)%psi(1:msh(nt)), &
                      rgrid(nt)%r(:), nrs, nrc, pow, msh(nt) )

                 CALL simpson ( msh(nt), aux, rgrid(nt)%rab, s(ih,jh) )

                 !<apsi>
                 IF ( ih > jh ) THEN
                    IF ( ABS ( ABS ( s(ih,jh) ) - 1.0_dp ) < 1.e-5_dp ) THEN
                       WRITE ( stdout, '(5X,2A,/,5X,A,I3,A,3I2,F12.8)' ) &
                            "init_gipaw_1: ", &
                            "projectors linearly dependent:", &
                            "ntyp =", nt, ", l/n1/n2 = ", l, ih, jh, &
                            s(ih,jh)
                       FLUSH ( stdout )
                       CALL errore ( "init_gipaw_1", &
                            "two projectors are linearly dependent", +1 )
                    ELSE IF ( ABS ( ABS ( s(ih,jh) ) - 1.0_dp ) < 1.e-2_dp ) &
                         THEN
                       IF ( n_overlap_warnings == 0 ) THEN
                          WRITE ( stdout, '(A)' ) ""
                       END IF
                       n_overlap_warnings = n_overlap_warnings + 1
                       WRITE ( stdout, '(5X,2A,/,5X,A,I3,A,3I2,F12.8)' ) &
                            "init_gipaw_1: ", &
                            "projectors nearly linearly dependent:", &
                            "ntyp =", nt, ", l/n1/n2 = ", l, ih, jh, s(ih,jh)
                       FLUSH ( stdout )
                    END IF
                 END IF
                 !</apsi>

              end do
           end do

#if defined ( __DEBUG_MINE_APSI )
           !<apsi>
           IF ( iverbatim > 5 ) THEN
              do ih = 1, paw_recon(nt)%paw_nl(l)
                 do jh = ih, paw_recon(nt)%paw_nl(l)
                    write( stdout, '( A, I3, 3I2, F12.7 )' ) &
                         "PROJ: ", nt, l, ih, jh, s(ih,jh)
                 end do
              end do
           END IF
           !</apsi>
#endif

           call invmat ( paw_recon(nt)%paw_nl(l), s, sinv)

           do ih = 1, paw_recon(nt)%paw_nl(l)
              n1 = paw_recon(nt)%paw_iltonh(l,ih)

              paw_recon(nt)%paw_betar(1:msh(nt),n1) = 0.0

              do jh = 1, paw_recon(nt)%paw_nl(l)
                 n2 = paw_recon(nt)%paw_iltonh(l,jh)

                 paw_recon(nt)%paw_betar(1:msh(nt),n1) &
                      = paw_recon(nt)%paw_betar(1:msh(nt),n1) &
                      + sinv(ih,jh) * paw_recon(nt)%psphi(n2)%psi(1:msh(nt))
              end do

              nrc = paw_recon(nt)%psphi(n1)%label%nrc
              nrs = paw_recon(nt)%psphi(n1)%label%nrs

              call step_f ( aux, &
                   paw_recon(nt)%paw_betar(1:msh(nt),n1),rgrid(nt)%r(:), &
                   nrs,nrc,pow,msh(nt))
              paw_recon(nt)%paw_betar(1:ndm,n1) = aux(1:ndm)
           end do
           deallocate (sinv)
           deallocate (s)

        end if
     end do
  end do
  IF ( n_overlap_warnings > 0 ) THEN
     WRITE ( stdout, '(A)' ) ""
  END IF

!    Check the orthogonality for projectors
!
!     nt=1
!     n1=paw_iltonh(0,1,1)
!     n2=paw_iltonh(0,2,1)

!     print *,n1,n2,nt
!     aux=paw_betar(:,n1,nt)*psphi(nt,n1)%psi
!     call simpson(msh (nt), aux, rab (1, nt), norm)
!     print *,'11',norm
!     aux=paw_betar(:,n1,nt)*psphi(nt,n2)%psi
!     call simpson(msh (nt), aux, rab (1, nt), norm)
!     print *,'12',norm
!     aux=paw_betar(:,n2,nt)*psphi(nt,n2)%psi
!     call simpson(msh (nt), aux, rab (1, nt), norm)
!     print *,'11',norm


  !
  !  compute Clebsch-Gordan coefficients
  !

  call aainit ( lmaxx+1 )

  !
  !     fill the interpolation table tab
  !
  pref = fpi / sqrt ( omega )
  call divide ( intra_pool_comm, nqx, startq, lastq )
  do nt = 1, ntyp
     paw_recon(nt)%paw_tab (:,:) = 0.0_dp

     do nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        do iq = startq, lastq
           qi = ( iq - 1 ) * dq
           call sph_bes ( msh(nt), rgrid(nt)%r, qi, l, besr )
           do ir = 1, msh(nt)
              aux(ir) = paw_recon(nt)%paw_betar(ir,nb) &
                   * besr(ir) * rgrid(nt)%r(ir)
           end do
           call simpson ( msh(nt), aux, rgrid(nt)%rab, vqint )
           paw_recon(nt)%paw_tab(iq,nb) = vqint * pref
        end do
     end do

#if defined(__MPI)
     call mp_sum ( paw_recon(nt)%paw_tab(:,:), intra_pool_comm )
#endif

  end do

  ! initialize spline interpolation
  if ( spline_ps ) then
     allocate(xdata(nqx))
     do iq = 1, nqx
        xdata(iq) = (iq - 1) * dq
     end do
     do nt = 1, ntyp
        allocate ( paw_recon(nt)%paw_tab_d2y(nqx,paw_recon(nt)%paw_nbeta) )
        paw_recon(nt)%paw_tab_d2y = 0.0_dp
        do nb = 1, paw_recon(nt)%paw_nbeta
           l = paw_recon(nt)%aephi(nb)%label%l
           d1 = ( paw_recon(nt)%paw_tab(2,nb) - paw_recon(nt)%paw_tab(1,nb) ) &
                / dq
           call spline ( xdata, paw_recon(nt)%paw_tab(:,nb), 0.0_dp, d1, &
                paw_recon(nt)%paw_tab_d2y(:,nb) )
        end do
     end do
     deallocate ( xdata )
  end if

  deallocate (besr)
  deallocate (aux1)
  deallocate (aux)

  call stop_clock ('init_gipaw_1')
  return

end subroutine init_gipaw_1




subroutine step_f(f2,f,r,nrs,nrc,pow,mesh)

  use kinds , only : dp

  !
  ! This routine apply a fonction which go smoothly to zero from rs to rc
  !

  implicit none
  integer :: mesh
  real(DP), Intent(out):: f2(mesh)
  real(DP), Intent(in) :: f(mesh), r(mesh)
  real(DP), Intent(in) :: pow
  integer :: nrs, nrc

  Integer :: i
  real(DP) :: rcp, rsp

  rcp = r(nrc)
  rsp = r(nrs)

      Do i=1,mesh
       If(r(i).Le.rsp) Then
          f2(i) = f(i)
       Else
          If(r(i).Le.rcp) Then
             f2(i)=f(i)* (1.d0-3.d0*((r(i)-rsp)/(rcp-rsp))**2+ &
                  2.d0*((r(i)-rsp)/(rcp-rsp))**3)**pow
          Else
             f2(i)=0.d0
          End If
       End If

    End Do

  End subroutine step_f
