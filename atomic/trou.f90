!
! ---------------------------------------------------------------
      subroutine find_coefficients &
           (ik,psi,energy,r,dx,aenorm,vpot,c,c2,lam,ndm)
!     -----------------------------------
!
!     recherche des coefficients du polynome
      implicit none
      integer, parameter :: dp=kind(1.d0)
      integer ndm, lam, ipvt(6), ik
      real(kind=dp):: vpot(ndm), psi(ndm), r(ndm), dx, energy, &
           c(6), c2, amat(6,6), y(6), rc, aenorm
      integer i, info, n, nmax, ndiv
      real(kind=dp):: c2o, dc2, newvalue, oldvalue, precision, &
      funz, rndm
!
      do i = 1,6
         c(i) = 0.d0
      end do
      rc  = r(ik)

      call fill_matrix(amat,rc,lam)
! calcul de mat = lu
!      call dgef(amat,6,6,ipvt,info)
      call DGETRF(6,6,amat,6,ipvt,info)
! calcul de y pour ax = y
      call eval_coeff(r,psi,ik,lam,energy,dx,vpot,y)
! cherche la valeur de  c2 qui verifie l'equ 29a
! This is done by miniming with a random search funz**2
! We are looking for the smallest solution
      c2o =  0.0d0        ! starting point
      dc2 =  0.1d0        ! starting range
      n   = 0             ! number of failed attempts
      nmax=2000
      ndiv= 200
!
!     original
!      precision = 1e-15
      precision = 7.e-10
      oldvalue = funz(amat,ipvt,y,rc,ik,aenorm,c2o,c,c2, &
                      lam,r,dx,ndm)**2
 10   continue
         c2 = c2o + (0.5 - rndm())*dc2
         newvalue = funz(amat,ipvt,y,rc,ik,aenorm,c2,c,c2,lam, &
                                         r,dx,ndm)**2
         if (newvalue.lt.precision) goto 900
         if (newvalue.lt.oldvalue) then
            n=0
            c2o = c2
            oldvalue=newvalue
         else
            n=n+1
! after ndiv failed attempts reduce the size of the interval
            if (n.gt.ndiv) dc2=dc2/10.d0
! after nmax failed attempts exit
            if (n.gt.nmax) then
               write(6,110) nmax, precision
110            format ('Warning: after ',i5,' attempts ',  &
                    'the error is still',e10.4)
               c2=c2o
               newvalue = funz(amat,ipvt,y,rc,ik, &
                              aenorm,c2,c,c2,lam,r,dx,ndm)**2
               goto 900
            end if
         end if
      go to 10
   
900   continue
      return
      end

!
! ---------------------------------------------------------------
      subroutine fill_matrix(amat,rc,l)
! ---------------------------------------------------------------
!
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp) :: amat(6,6),rc
      integer :: l
!      
!     this routine fills the matrix withn the coefficients taken
!     from  p,p',p'',p''', p(iv), where p is
!     p(r)= c0 + c4 r^4 + c6 r^6 + c8 r^8 + ...
      integer pr1(6),cr1(6),pr(6),cr(6),i,j
! matrices representing the coefficients (cr) and the powers of r ( pr)
      data pr1/0,4,6,8,10,12/
      data cr1/1,1,1,1,1,1/

      do i=1,6
         pr(i) = pr1(i)
         cr(i) = cr1(i)
      end do
      do i = 1,5
!     fill matrix row by row
         do j = 1,6
            amat(i,j) = cr(j) * rc**(pr(j)*1.d0)
         end do
!     derivate polynomial expression
         do j = 1,6 
            cr(j) = cr(j) * pr(j)
            pr(j) = max(0,pr(j)-1)
         end do
      end do
      do j = 1,6
         amat(6,j) = 0.d0
      end do
      amat(6,2) = 2.d0*l +5.d0

      return
      end
!
! ----------------------------------------------------------------
      subroutine eval_coeff (r,psi,ik,l,energy,dx,vpot,y)
! ----------------------------------------------------------------
!    calcule les coefficients dependant de la fct d'onde calculee
!    avec tous les electrons. ces coefficients ervent a la resolution
!    du systeme lineaire.
!    en entree : ik,nx,vpot comme dans le programme principal
!    en sortie : une matrice colonne y contenant les coefficients
!
      implicit none
      integer,parameter :: dpp=kind(1.d0)
      integer :: ik,l,lp1
      real(kind=dpp):: r(ik+3),psi(ik+3),vpot(ik+3),energy,dx,y(6)
      real(kind=dpp) p,dp,d,vae,dvae,ddvae,rc,rc2,rc3
      real(kind=dpp)   deriv_7pts,deriv2_7pts
      external deriv_7pts,deriv2_7pts
!
!   evaluer p et p' ( dp )
      rc = r(ik)
      p  = psi(ik)
      dp = deriv_7pts(psi,ik,rc,dx)
      if (p.lt.0.0) then
         p  = -p
         dp = -dp
      endif
      d = dp /p
!   evaluer vae, dvae, ddvae
      vae   = vpot(ik)
      dvae  =  deriv_7pts(vpot,ik,rc,dx)
      ddvae = deriv2_7pts(vpot,ik,rc,dx)
!   calcul de parametres intervenant dans les calculs successifs
      lp1 = l + 1
      rc2 = rc * rc
      rc3 = rc2* rc
!
      y(1) = log ( p / rc**lp1 )
      y(2) = dp/p - (lp1 / rc)
      y(3) = vae - energy + (lp1*lp1)/rc2 - d*d
      y(4) = dvae - 2.d0*(vae - energy + l*lp1/rc2)*d - &
                   2.d0*(lp1*lp1)/rc3 &
                  + 2.d0*(d*d*d)
      y(5) = ddvae - 2.d0*(dvae - 2.d0*l*lp1/rc3)*d  &
             + 6.d0*lp1*lp1/(rc3*rc) &
             - 2.d0*(vae - energy + l*lp1/rc2 - 3.d0*d*d)* &
                 (vae - energy + l*lp1/rc2 - d*d)
      return
      end
! ----------------------------------------------------------------
       function deriv2_7pts(f,ik,rc,h)
! ----------------------------------------------------------------
!
!      evaluates the second derivative of function f, the function
!      is given numerically on logarithmic mesh r.
!      nm = dimension of mesh
!      ik = integer : position of the point in which the derivative
!      will be evaluated.
!      h is distance between x(i) and x(i+1) where
!      r(i) = exp(x(i))/znesh & r(j) = exp(x(j))/znesh
!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       integer :: a(7),n,nm,i,ik
       real(kind=dp) :: f(ik+3),rc,h,sum,sum1,deriv_7pts,deriv2_7pts
!      coefficients for the formula in abramowitz & stegun p.914
!      the formula is used for 7 points.
       data a/4,-54,540,-980,540,-54,4/   ! these are coefficients

! formula for linear mesh 
       sum = 0.d0
       do i=1,7
          sum = sum + a(i)*f(i-4+ik)
       end do
       sum = 2.d0*sum/(720.d0*h**2)
! transform to logarithmic mesh
       sum1 = deriv_7pts(f,ik,rc,h) 
       deriv2_7pts = sum/(rc*rc) - sum1 /rc

       return
       end
! ---------------------------------------------------------------     
       function deriv_7pts(f,ik,rc,h)
! ---------------------------------------------------------------
!      evaluates the first derivative of function f, the function
!      is given numerically on logarithmic mesh r.
!      nm = dimension of mesh
!      ik = integer : position of the point in which the derivative
!      will be evaluated.
!      h is distance between x(i) and x(i+1) where
!      r(i) = exp(x(i))/znesh & r(j) = exp(x(j))/znesh
!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       integer :: a(7),n,ik,i
       real(kind=dp) :: f(ik+3),rc,h,sum,deriv_7pts
!      coefficients for the formula in abramowitz & stegun p.914
       data a/-12,108,-540,0,540,-108,12/

! formula for linear mesh 
       sum = 0
       do i=1,7
          sum = sum + a(i)*f(i-4+ik)
       end do
       deriv_7pts = sum/(720.0*h)
! transform to logarithmic mesh
       deriv_7pts = deriv_7pts /rc

       return
       end

! --------------------------------------------------------
      function funz(amat,ipvt,y,rc,ik,aenorm,x,c,c2,  &
                 lam,r,dx,ndm)
! --------------------------------------------------------
!    cette fonction est la fonction de c2 qu'il faut annuler
!    pour trouver une valeur de c2 qui verifie la conservation 
!    de la charge de coeur.
!    cette fonction calcule de plus a chaque fois les ci qui 
!    verifient les equations lineaires avec un c2 donne.
!    en entree : une valeur de c2 donnee dans x
!    en sortie, la valeur de la fonction pour cette valeur de x
!    : c'est la fonction qui correspond a l'equation integrale
!
      implicit none
      integer,parameter :: dp=kind(1.d0)
      integer :: ndm,lam,ipvt(6),ik,mesh,i,istart
      real(kind=dp) :: funz, x, c(6), c2, amat(6,6), y(6), rc, aenorm, &
           r(ndm), dx, chip2, psnorm, f0, f1, f2
      external chip2


!    resolution du systeme lineaire pour cette valeur de x (=c2)
      c(1) = y(1) - x*rc**2  ! ajoute les termes en c2
      c(2) = y(2) - 2*x*rc
      c(3) = y(3) - 2*x
      c(4) = y(4)       ! pas de coeff en c2
      c(5) = y(5)
      c(6) = -x*x
!      call dges(amat,6,6,ipvt,c,0)    ! resoud le systeme
      call DGETRS('N',6,1,amat,6,ipvt,c,6,0)    ! resoud le systeme
!    calcul de la norme de la pseudo-fonction d'onde
      psnorm = 0.0
      istart= 2-mod(ik,2)
      f2 = r(istart) * chip2(c,x,lam,r(istart))
      do i = istart,ik-2,2
         f0 = f2
         f1 = r(i+1)*chip2(c,x,lam,r(i+1))
         f2 = r(i+2)*chip2(c,x,lam,r(i+2))
         psnorm = psnorm + f0 + 4*f1 + f2
      end do
      psnorm = psnorm * dx/3.0 + r(1)**(2*lam+3)/(2*lam+3)
      funz = log( psnorm / aenorm )
!
      return
      end

! --------------------------------------------------------
      function chip2(c,c2,l,r)
! --------------------------------------------------------
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp):: chip2,c(6),c2,r,r2
      integer :: l

      r2 = r**2
      chip2 = r2**(l+1) * exp(2.d0* &
       (((((((c(6)*r2+c(5))*r2+c(4))*r2+c(3))*r2+c(2))*r2+c2)*r2)+c(1)))

      return
      end

! ----------------------------------------------------------------
       function der3num(r,f,ik,nm,h)
! ----------------------------------------------------------------
!     calcule la 3eme derivee d'une fonction donnee numeriquement
!     sur un maillage logarithmique.
!     en entree : r maillage; f fonction ; ik : position de l'eval.
!                 nm : dimension de f et r ; h comme dx dans 
!                 le programme principal
!     attention ! precision pas fantastique ! 
       implicit none
       integer,parameter :: dp=kind(1.d0)
       integer :: nm, ik,i
       real(kind=dp)::  r(nm),f(nm), y(7), h, deriv_7pts, deriv2_7pts
       real(kind=dp) :: der3num

       do i=1,7
          y(i) =  deriv2_7pts(f,i+ik-4,r(i+ik-4),h)
       end do
       der3num = deriv_7pts(y,4,r(ik),h)

       return
       end
! ----------------------------------------------------------------
       function der4num(r,f,ik,nm,h)
! ----------------------------------------------------------------
!     idem que der3num, mais pour la 4eme derivee
!     attention ! precision pas terrible !!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       integer:: nm,i,ik
       real(kind=dp) :: der4num, r(nm),f(nm),y(7),h,deriv2_7pts

       do i=1,7
          y(i) = deriv2_7pts(f,i+ik-4,r(i+ik-4),h)
       end do
       der4num= deriv2_7pts(y,4,r(ik),h)

       return
       end

! ----------------------------------------------------------------
       function der3an(l,c,c2,rc)
! ----------------------------------------------------------------
!
! 3rd derivative of r^(l+1) e^p(r)
!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       integer :: l
       real(kind=dp):: c(6), c2, rc,pr,dexpr,d2expr,d3expr, der3an
       external pr,dexpr,d2expr,d3expr

       der3an= (l-1)*l*(l+1)*rc**(l-2) *exp(pr(c,c2,rc)) +     &
              3.d0 *    l*(l+1)*rc**(l-1) * dexpr(c,c2,rc)  +  &
              3.d0 *      (l+1)*rc**(l  ) *d2expr(c,c2,rc)  +  &
                             rc**(l+1) *d3expr(c,c2,rc) 

       return
       end

! ----------------------------------------------------------------
       function der4an(l,c,c2,rc)
! ----------------------------------------------------------------
!
! 4rd derivative of r^(l+1) e^p(r)
!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       integer l
       real(kind=dp):: c(6), c2, rc,pr,dexpr,d2expr,d3expr,d4expr,der4an
       external pr,dexpr,d2expr,d3expr,d4expr

       der4an = (l-2)*(l-1)*l*(l+1)*rc**(l-3) * exp(pr(c,c2,rc)) +   &
               4.d0 *    (l-1)*l*(l+1)*rc**(l-2) *  dexpr(c,c2,rc)  + &
               6.d0 *          l*(l+1)*rc**(l-1) * d2expr(c,c2,rc)  + &
               4.d0 *            (l+1)*rc**(l  ) * d3expr(c,c2,rc)  + &
                                    rc**(l+1) * d4expr(c,c2,rc)

       return
       end
! ----------------------------------------------------------------
       function dexpr(c,c2,rc)
! ----------------------------------------------------------------
!
! 1st derivative of e^p(r)
!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       real(kind=dp) ::  c(6), c2, rc, pr, dpr, dexpr
       external pr, dpr

       dexpr= exp(pr(c,c2,rc)) * dpr(c,c2,rc)

       return
       end
! ----------------------------------------------------------------
       function d2expr(c,c2,rc)
! ----------------------------------------------------------------
!
! 2nd derivative of e^p(r)
!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       real(kind=dp) ::  c(6), c2, rc, pr, dpr, d2pr, d2expr
       external pr, dpr, d2pr

       d2expr= exp(pr(c,c2,rc)) * ( dpr(c,c2,rc)**2+d2pr(c,c2,rc) )

       return
       end
! ----------------------------------------------------------------
       function d3expr(c,c2,rc)
! ----------------------------------------------------------------
!
! 3nd derivative of e^p(r)
!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       real(kind=dp)::  c(6), c2, rc, pr, dpr, d2pr, d3pr, d3expr
       external pr, dpr, d2pr, d3pr

       d3expr= exp(pr(c,c2,rc)) * (   dpr(c,c2,rc)**3  &
                                  +  3*dpr(c,c2,rc)*d2pr(c,c2,rc) &
                                  +   d3pr(c,c2,rc)               )

       return
       end
! ----------------------------------------------------------------
       function d4expr(c,c2,rc)
! ----------------------------------------------------------------
!
! 4th derivative of e^p(r)
!
       implicit none
       integer,parameter :: dp=kind(1.d0)
       real(kind=dp)::  c(6), c2, rc, pr,  dpr, d2pr, d3pr, d4pr, d4expr
       external pr, dpr, d2pr, d3pr, d4pr

       d4expr= exp(pr(c,c2,rc)) * (     dpr(c,c2,rc)**4 &
                                  +  6.d0* dpr(c,c2,rc)**2*d2pr(c,c2,rc) &
                                  +  3.d0*d2pr(c,c2,rc)**2 &
                                  +  4.d0* dpr(c,c2,rc)*d3pr(c,c2,rc) &
                                  +    d4pr(c,c2,rc)                  )

       return
       end
!
! --------------------------------------------------------
      function pr(c,c2,x)
! --------------------------------------------------------
!     cette fonction evalue le polynome dont les coefficients
!     sont dans c d'apres la forme suivante
!     p(x) =  c(2)*x^4 + c(3) x^6 + c(4) x^8 + c(5) x^10
!            +c(6) x^12 + c2*x^2 + c(1)
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp) :: c(6), c2, x, y, pr

      y = x*x
      pr = (((((c(6)*y+c(5))*y+c(4))*y+c(3))*y+c(2))*y+c2)*y &
                + c(1)

      return
      end
!
      function dpr(c,c2,x)
! --------------------------------------------------------
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp) :: c(6),c2,x,dpr

      dpr  = 2*c2*x + 4*c(2)*x**3 + 6*c(3)*x**5 + 8*c(4)*x**7 &
           + 10*c(5)*x**9 + 12*c(6)*x**11

      return
      end
! 
      function d2pr(c,c2,x)
! --------------------------------------------------------
      implicit none
      integer,parameter :: dp=kind(1.d0)
      real(kind=dp) :: c(6),c2,x,d2pr

      d2pr = 2*c2 + 12*c(2)*x**2 + 30*c(3)*x**4 + 56*c(4)*x**6 &
           + 90*c(5)*x**8 + 132*c(6)*x**10

      return
      end
!
      function d3pr(c,c2,x)
! --------------------------------------------------------
      implicit none
       integer,parameter :: dp=kind(1.d0)
      real(kind=dp) :: c(6),c2,x,d3pr

      d3pr  = 24*c(2)*x + 120*c(3)*x**3 + 336*c(4)*x**5 &
           + 720*c(5)*x**7 + 1320*c(6)*x**9

      return
      end
! 
      function d4pr(c,c2,x)
! --------------------------------------------------------
      implicit none
       integer,parameter :: dp=kind(1.d0)
      real(kind=dp) :: c(6),c2,x,d4pr

      d4pr  = 24.d0*c(2) + 360.d0*c(3)*x**2 + 1680.d0*c(4)*x**4 &
           + 5040.d0*c(5)*x**6 + 11880.d0*c(6)*x**8

      return
      end
