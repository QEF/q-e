!----------------------------------------------------------------------
      subroutine calc_btq(ql,qr_k,idbes)
!----------------------------------------------------------------------
!
!   Calculates the Bessel-transform (or its derivative if idbes=1) 
!   of the augmented qrad charges at a given ql point.
!   Rydberg atomic units are  used.
!
      use atom, only: r, rab, dx
      USE ions_base, ONLY : ntyp => nsp
      use cell_base, only: omega
      USE parameters, only:  ndmx, nbrx
      USE kinds, only: DP
      use constants, only: fpi
      USE uspp_param, only: lmaxq, qfunc, qfcoef, nqf, rinner, lll, &
           nbeta, kkbeta, tvanp
      !
      implicit none
      !
      integer :: ik,  msh_bp, i, np, m, k, l
      integer :: n,idbes,ilmin,ilmax,iv,jv
      real(kind=DP)  :: jl(ndmx), ql, sum, jlp1(ndmx), aux(ndmx), &
             qr_k(nbrx,nbrx,lmaxq,ntyp)

! declaration readvan quantities
!      integer NBETA,KKBETA,iver,nqf,ifqopt,nqlc,lll
!      real*8 DION,BETAR,QQQ,QFUNC,qfcoef,rinner
!      COMMON/NCPRM/DION(NBRX,NBRX,NPSX),
!     C           BETAR(0:ndm,NBRX,NPSX),QQQ(NBRX,NBRX,NPSX),
!     C           QFUNC(0:ndm,NBRX,NBRX,NPSX),
!     C           NBETA(NPSX),KKBETA(NPSX),NVALES(NPSX),lll(nbrx,npsx),
!     C           iver(3,npsx),nqf(npsx),ifqopt(npsx),nqlc(npsx),
!     C           qfcoef(nqfx,lmaxq,NBRX,NBRX,npsx),rinner(lmaxq,npsx)
!      common/ncprm/dion(nbrx,nbrx,npsx),
!     +           betar(0:ndm,nbrx,npsx), qqq(nbrx,nbrx,npsx),
!     +           qfunc(0:ndm,nbrx,nbrx,npsx),
!     +           qfcoef(nqfx,lmaxq,nbrx,nbrx,npsx), rinner(lmaxq,npsx),
!     +           nbeta(npsx), kkbeta(npsx),
!     +           nqf(npsx), nqlc(npsx), ifqopt(npsx), lll(nbrx,npsx),
!     +           iver(3,npsx)


!
      do np=1,ntyp
         msh_bp=kkbeta(np)
         if (tvanp(np)) then
            do iv =1, nbeta(np)
               do jv =iv, nbeta(np)
                  ilmin = iabs(lll(iv,np)-lll(jv,np))
                  ilmax = iabs(lll(iv,np)+lll(jv,np))
!       only need to calculate for for lmin,lmin+2 ...lmax-2,lmax
                  do l = ilmin,ilmax,2
                     do i =  msh_bp,2,-1
                        if (r(i,np) .lt. rinner(l+1,np)) goto 100
                        aux(i) = qfunc(i,iv,jv,np)
                     enddo
 100                 call setqf(qfcoef(1,l+1,iv,jv,np),aux(1),r(1,np) &
                         ,nqf(np),l,i)

                        if (idbes .eq. 1) then
                           call dbess(ql,l+1,msh_bp,r(1,np), &
                               jl)
                        else
                           call bess(ql,l+1,msh_bp,r(1,np), &
                               jl)
                        endif

! jl is the Bessel function (or its derivative) calculated at ql
! now integrate qfunc*jl*r^2 = Bessel transform of qfunc

                        do i=1, msh_bp
                           jlp1(i) = jl(i)*aux(i)
                        enddo
!                        if (tlog(np)) then
                           if(tvanp(np)) then
                              call radlg1(msh_bp,jlp1,rab(1,np),sum) 
                           else
                              call radlg(msh_bp,jlp1,r(1,np),dx(np),sum)
                           endif
!                        else
!                           call radin(msh_bp,dx(np),jlp1,sum)
!                        endif
                        qr_k(iv,jv,l+1,np) = sum*fpi/omega
                        qr_k(jv,iv,l+1,np) = qr_k(iv,jv,l+1,np)
                   
!c                       WRITE( stdout,*) 'qr_k=',qr_k(iv,jv,l+1,np)

                  end do
               end do
            enddo
         endif
      enddo
!
      return
      end
