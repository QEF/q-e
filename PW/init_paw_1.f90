!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
subroutine init_paw_1
  !----------------------------------------------------------------------
  !
  ! This routine initialize the variables of the paw projector
  ! and create the projectors in radial part (paw_betar) 
  !
  USE kinds ,      ONLY : dp
  USE parameters , ONLY : lqmax , nbrx, lmaxx
  USE cell_base ,  ONLY : omega
  USE ions_base,   ONLY : nat, ntyp => nsp, ityp
  USE constants ,  ONLY : fpi
#ifdef USE_SPLINES
  USE us,          ONLY : nqxq, dq, nqx, tab, tab_d2y, qrad
  USE paw ,        ONLY : paw_nhm, paw_nh, paw_lmaxkb, paw_nkb, paw_nl, &
                          paw_iltonh, paw_tab, aephi, paw_betar, psphi, &
                          paw_indv, paw_nhtom, paw_nhtol, paw_nbeta, &
                          paw_tab_d2y
  USE splinelib
#else
  USE us,          ONLY : nqx, dq
  USE paw ,        ONLY : paw_nhm, paw_nh, paw_lmaxkb, paw_nkb, paw_nl, &
                          paw_iltonh, paw_tab, aephi, paw_betar, psphi, &
                          paw_indv, paw_nhtom, paw_nhtol, paw_nbeta 
#endif
  USE uspp,        ONLY : ap, aainit
  USE atom ,       ONLY : r, rab, msh
  !
  implicit none
  !
  !     here a few local variables
  !

  integer :: nt, ih, jh, nb, l, m, ir, iq, startq, &
       lastq, na, j, n1, n2, ndm, nrs, nrc
  ! various counters
  real(DP), allocatable :: aux (:), aux1 (:), besr (:)
  ! various work space
  real(DP) :: prefr, pref, qi, norm
  ! the prefactor of the q functions
  ! the prefactor of the beta functions
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  real(DP), allocatable :: ylmk0 (:), s(:,:), sinv(:,:)
  ! the spherical harmonics
  real(DP) ::  vqint
  ! interpolated value

  real(DP) rc,rs,pow
  
#ifdef USE_SPLINES
  integer :: paw_nbeta_max
  real(DP), allocatable :: xdata(:)
  real(DP) :: d1
#endif
  
  call start_clock ('init_paw_1')
  !
  !    Initialization of the variables
  !
  ndm = MAXVAL (msh(1:ntyp))
  paw_nhm = 0
  paw_nh = 0
  paw_lmaxkb = 0
  do nt = 1, ntyp
     do nb = 1, paw_nbeta (nt)
        paw_nh (nt) = paw_nh (nt) + 2 * aephi(nt,nb)%label%l + 1
        paw_lmaxkb = max (paw_lmaxkb,  aephi(nt,nb)%label%l)
     enddo
     if (paw_nh (nt) .gt.paw_nhm) paw_nhm = paw_nh (nt)
  enddo

  allocate (aux ( ndm))    
  allocate (aux1( ndm))    
  allocate (besr( ndm))    
  allocate (ylmk0( (paw_lmaxkb+1) ** 2 ))    
  allocate (paw_nhtol(paw_nhm, ntyp))
  allocate (paw_nhtom(paw_nhm, ntyp))
  allocate (paw_indv(paw_nhm, ntyp))
  allocate (paw_tab(nqx, nbrx, ntyp))
  allocate (paw_nl(0:paw_lmaxkb, ntyp))
  allocate (paw_iltonh(0:paw_lmaxkb,paw_nhm, ntyp))

  ! calculate the number of beta functions of the solid
  !
  paw_nkb = 0
  do na = 1, nat
     paw_nkb = paw_nkb + paw_nh (ityp(na))
  enddo

 
  prefr = fpi / omega
  !
  !   For each pseudopotential we initialize the indices nhtol, nhtom,
  !   indv, 
  !
  paw_nl=0
  paw_iltonh=0
  do nt = 1, ntyp
        ih = 1
        do nb = 1, paw_nbeta (nt)
           l = aephi(nt,nb)%label%l
           paw_nl(l,nt) = paw_nl(l,nt) + 1
           paw_iltonh(l,paw_nl(l,nt) ,nt)= nb
            do m = 1, 2 * l + 1
              paw_nhtol (ih, nt) = l
              paw_nhtom (ih, nt) = m
              paw_indv (ih, nt) = nb
              ih = ih + 1
           enddo
        enddo

  ! Rescale the wavefunctions so that int_0^rc f|psi|^2=1
  ! 
  !      rc=1.6d0

        pow=1.d0
        do j = 1, paw_nbeta (nt)
           rc=psphi(nt,j)%label%rc
           rs=1.d0/3.d0*rc
           nrc =  Count(r(1:msh(nt),nt).le.rc)
           nrs =  Count(r(1:msh(nt),nt).le.rs)
           psphi(nt,j)%label%nrc = nrc
           aephi(nt,j)%label%nrc = nrc
           call step_f(aux,psphi(nt,j)%psi**2,r(:,nt),nrs,nrc,pow,msh(nt))
           call simpson (msh (nt), aux, rab (1, nt), norm )
           
           psphi(nt,j)%psi = psphi(nt,j)%psi/ sqrt(norm)
           aephi(nt,j)%psi = aephi(nt,j)%psi / sqrt(norm)

        enddo
        
        !
        !   calculate the overlap matrix
        !

        aux=0.d0
        do l=0,paw_lmaxkb
          if (paw_nl(l,nt)>0) then
           allocate (s(paw_nl(l,nt),paw_nl(l,nt)))
           allocate (sinv(paw_nl(l,nt),paw_nl(l,nt)))
           do ih=1,paw_nl(l,nt)
              n1=paw_iltonh(l,ih,nt)
              do jh=1,paw_nl(l,nt)
                 n2=paw_iltonh(l,jh,nt)
                 call step_f(aux,psphi(nt,n1)%psi(1:msh(nt)) * &
                      psphi(nt,n2)%psi(1:msh(nt)),r(:,nt),nrs,nrc,pow,msh(nt))
                 call simpson (msh (nt), aux, rab (1, nt), s(ih,jh))
              enddo
           enddo
           call invmat (paw_nl(l,nt), s, sinv, norm) 

           do ih=1,paw_nl(l,nt)
              n1=paw_iltonh(l,ih,nt)
              do jh=1,paw_nl(l,nt)
                 n2=paw_iltonh(l,jh,nt)
                 
                 paw_betar(1:msh(nt),n1,nt)=paw_betar(1:msh(nt),n1,nt)+ &
                      sinv(ih,jh) * psphi(nt,n2)%psi(1:msh(nt))
              enddo
              call step_f(aux, &
                   paw_betar(1:msh(nt),n1,nt),r(:,nt),nrs,nrc,pow,msh(nt))
              paw_betar(1:ndm,n1,nt)=aux(1:ndm)
           enddo
           deallocate (sinv)
           deallocate (s)

          endif
        enddo
     enddo

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

  call aainit (lmaxx+1)

  !
  !     fill the interpolation table tab
  !
  pref = fpi / sqrt (omega)
  call divide (nqx, startq, lastq)
  paw_tab (:,:,:) = 0.d0
  do nt = 1, ntyp
     do nb = 1, paw_nbeta (nt)
        l = aephi(nt, nb)%label%l
        do iq = startq, lastq
           qi = (iq - 1) * dq
           call sph_bes (msh(nt), r (1, nt), qi, l, besr)
           do ir = 1, msh(nt)
              aux (ir) = paw_betar (ir, nb, nt) * besr (ir) * r (ir, nt)
           enddo
           call simpson (msh (nt), aux, rab (1, nt), vqint)
           paw_tab (iq, nb, nt) = vqint * pref
        enddo
     enddo
  enddo

#ifdef __PARA
  call reduce (nqx * nbrx * ntyp, paw_tab)
#endif

#ifdef USE_SPLINES
  !<ceres>
  ! initialize spline interpolation
  paw_nbeta_max = maxval ( paw_nbeta (:) )
  allocate ( paw_tab_d2y ( nqx, paw_nbeta_max, ntyp ) )
  
  paw_tab_d2y = 0.0
  
  allocate(xdata(lastq-startq+1))
  do iq = startq, lastq
    xdata(iq) = (iq - 1) * dq
  enddo
  do nt = 1, ntyp
     do nb = 1, paw_nbeta (nt)
        l = aephi(nt, nb)%label%l
        d1 = (paw_tab(2,nb,nt) - paw_tab(1,nb,nt)) / dq
        call spline(xdata, paw_tab(:,nb,nt), 0.d0, d1, paw_tab_d2y(:,nb,nt))
     enddo
  enddo
  deallocate(xdata)
  !</ceres>
#endif
  
  deallocate (ylmk0)
  deallocate (besr)
  deallocate (aux1)
  deallocate (aux)

  call stop_clock ('init_paw_1')
  return

end subroutine init_paw_1

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
          Endif
       Endif

    End Do

  End subroutine step_f
