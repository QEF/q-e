!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine summary
  !-----------------------------------------------------------------------
  !
  !    This routine writes on output all the information obtained from
  !    the input file and from the setup routine, before starting the
  !    self-consistent calculation.
  !
  !    if iverbosity = 0 only a partial summary is done.
  !
#include "machine.h"
  USE io_global,  ONLY :  stdout
  USE kinds, ONLY: DP
  USE constants, ONLY: amconv
  USE atom
  USE cell_base
  USE basis
  USE char, ONLY: title, sname
  USE cellmd, ONLY: calc, cmass
  USE dynam, ONLY: amass
  USE gvect
  USE gsmooth
  USE lsda_mod, ONLY: lsda, starting_magnetization
  USE klist
  USE ktetra
  USE pseud, ONLY: zp, alps, alpc, cc, aps, nlc, nnl, lmax, lloc, &
       a_nlcc, b_nlcc, alpha_nlcc
  USE symme, ONLY: nsym, invsym, s, ftau
  USE control_flags
  USE us, ONLY: tvanp
  USE uspp_param, ONLY: nqf, rinner, nqlc, nbeta, iver, lll, psd
  USE spin_orb, only: lspinorb
  USE funct
  implicit none
  !
  !     declaration of the local variables
  !
  integer :: i, ipol, apol, na, isym, ik, ib, nt, l, ngmtot
  ! counter on the celldm elements
  ! counter on polarizations
  ! counter on direct or reciprocal lattice vect
  ! counter on atoms
  ! counter on symmetries
  ! counter on k points
  ! counter on beta functions
  ! counter on types
  ! counter on angular momenta
  ! total number of G-vectors (parallel executio
  real(kind=DP) :: sr (3, 3), ft1, ft2, ft3
  ! symmetry matrix in real axes
  ! fractionary translation
  real(kind=DP), allocatable :: xau (:,:)
  ! atomic coordinate referred to the crystal axes
  real(kind=DP) :: xkg (3)
  ! coordinates of the k point in crystal axes
  character :: mixing_style * 9
  character :: ps * 5
  ! name of pseudo type
  real(kind=DP) :: xp
  ! fraction contributing to a given atom type (obsolescent)
  !
  !     we start with a general description of the run
  !
  if (imix.eq.-1) mixing_style = 'potential'
  if (imix.eq. 0) mixing_style = 'plain'
  if (imix.eq. 1) mixing_style = 'TF'
  if (imix.eq. 2) mixing_style = 'local-TF'

  if (title.ne.' ') then
     WRITE( stdout,"(/,5x,'Title: ',/,5x,a75)") title
  end if
  WRITE( stdout, 100) ibrav, alat, omega, nat, ntyp, &
       ecutwfc, dual * ecutwfc, tr2, mixing_beta, nmix, &
       mixing_style

100 format (/,/,5x, &
       &     'bravais-lattice index     = ',i12,/,5x, &
       &     'lattice parameter (a_0)   = ',f12.4,'  a.u.',/,5x, &
       &     'unit-cell volume          = ',f12.4,' (a.u.)^3',/,5x, &
       &     'number of atoms/cell      = ',i12,/,5x, &
       &     'number of atomic types    = ',i12,/,5x, &
       &     'kinetic-energy cutoff     = ',f12.4,'  Ry',/,5x, &
       &     'charge density cutoff     = ',f12.4,'  Ry',/,5x, &
       &     'convergence threshold     = ',1pe12.1,/,5x, &
       &     'beta                      = ',0pf12.4,/,5x, &
       &     'number of iterations used = ',i12,2x,a,' mixing')
  WRITE( stdout, '(5x,"Exchange-correlation      = ",a, &
       &       " (",4i1,")")') trim(dft) , iexch, icorr, igcx, igcc
  if (iswitch.gt.0) then
     WRITE( stdout, '(5x,"iswitch = ",i2,"  nstep  = ",i4,/)') iswitch, nstep
  else
     WRITE( stdout, '(5x,"iswitch = ",i2/)') iswitch
  endif
  if (lspinorb) write(stdout, &
               '(5x,"Noncollinear calculation with spin-orbit",/)')

  if (qcutz.gt.0.d0) then
     WRITE( stdout, 110) ecfixed, qcutz, q2sigma
110  format   (5x,'A smooth kinetic-energy cutoff is imposed at ', &
          &             f12.4,' Ry',/5x,'height of the smooth ', &
          &             'step-function =',f21.4,' Ry',/5x, &
          &             'width of the smooth step-function  =',f21.4, &
          &             ' Ry',/)

  endif
  !
  !    and here more detailed information. Description of the unit cell
  !
  WRITE( stdout, '(2(3x,3(2x,"celldm(",i1,")=",f11.6),/))') &
       (i, celldm(i), i=1,6)
  WRITE( stdout, '(5x, &
       &     "crystal axes: (cart. coord. in units of a_0)",/, &
       &       3(15x,"a(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,  &
       (at (ipol, apol) , ipol = 1, 3) , apol = 1, 3)

  WRITE( stdout, '(5x, &
       &   "reciprocal axes: (cart. coord. in units 2 pi/a_0)",/, &
       &            3(15x,"b(",i1,") = (",3f10.6," )  ",/ ) )')  (apol,&
       &  (bg (ipol, apol) , ipol = 1, 3) , apol = 1, 3)
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        ps = '(US)'
        WRITE( stdout, '(/5x,"PSEUDO",i2," is ",a2, &
             &        1x,a5,"   zval =",f5.1,"   lmax=",i2, &
             &        "   lloc=",i2)') nt, psd (nt) , ps, zp (nt) , lmax (nt) &
             &, lloc (nt)
        WRITE( stdout, '(5x,"Version ", 3i3, " of US pseudo code")') &
             (iver (i, nt) , i = 1, 3)
        WRITE( stdout, '(5x,"Using log mesh of ", i5, " points")') mesh (nt)
        WRITE( stdout, '(5x,"The pseudopotential has ",i2, &
             &       " beta functions with: ")') nbeta (nt)
        do ib = 1, nbeta (nt)
           WRITE( stdout, '(15x," l(",i1,") = ",i3)') ib, lll (ib, nt)
        enddo
        WRITE( stdout, '(5x,"Q(r) pseudized with ", &
             &          i2," coefficients,  rinner = ",3f8.3,/ &
             &          52x,3f8.3,/ &
             &          52x,3f8.3)') nqf(nt), (rinner(i,nt), i=1,nqlc(nt) )
     else
        if (nlc (nt) .eq.1.and.nnl (nt) .eq.1) then
           ps = '(vbc)'
        elseif (nlc (nt) .eq.2.and.nnl (nt) .eq.3) then
           ps = '(bhs)'
        elseif (nlc (nt) .eq.1.and.nnl (nt) .eq.3) then
           ps = '(our)'
        else
           ps = '     '
        endif

        WRITE( stdout, '(/5x,"PSEUDO",i2," is ",a2, 1x,a5,"   zval =",f5.1,&
             &      "   lmax=",i2,"   lloc=",i2)') &
                        nt, psd(nt), ps, zp(nt), lmax(nt), lloc(nt)
        if (numeric (nt) ) then
           WRITE( stdout, '(5x,"(in numerical form: ",i5,&
                &" grid points",", xmin = ",f5.2,", dx = ",f6.4,")")')&
                & mesh (nt) , xmin (nt) , dx (nt)
        else
           WRITE( stdout, '(/14x,"i=",7x,"1",13x,"2",10x,"3")')
           WRITE( stdout, '(/5x,"core")')
           WRITE( stdout, '(5x,"alpha =",4x,3g13.5)') (alpc (i, nt) , i = 1, 2)
           WRITE( stdout, '(5x,"a(i)  =",4x,3g13.5)') (cc (i, nt) , i = 1, 2)
           do l = 0, lmax (nt)
              WRITE( stdout, '(/5x,"l = ",i2)') l
              WRITE( stdout, '(5x,"alpha =",4x,3g13.5)') (alps (i, l, nt) , &
                   i = 1, 3)
              WRITE( stdout, '(5x,"a(i)  =",4x,3g13.5)') (aps (i, l, nt) , i = 1,3)
              WRITE( stdout, '(5x,"a(i+3)=",4x,3g13.5)') (aps (i, l, nt) , i= 4, 6)
           enddo
           if ( nlcc(nt) ) WRITE( stdout, 200) a_nlcc(nt), b_nlcc(nt), alpha_nlcc(nt)
200        format(/5x,'nonlinear core correction: ', &
                &     'rho(r) = ( a + b r^2) exp(-alpha r^2)', &
                & /,5x,'a    =',4x,g11.5, &
                & /,5x,'b    =',4x,g11.5, &
                & /,5x,'alpha=',4x,g11.5)
        endif
     endif

  enddo
  WRITE( stdout, '(/5x, "atomic species   valence    mass     pseudopotential")')
  xp = 1.d0
  do nt = 1, ntyp
     if (calc.eq.' ') then
        WRITE( stdout, '(5x,a6,6x,f10.2,2x,f10.5,5x,5 (a2,"(",f5.2,")"))') &
                   atm(nt), zv(nt), amass(nt), psd(nt), xp
     else
        WRITE( stdout, '(5x,a6,6x,f10.2,2x,f10.5,5x,5 (a2,"(",f5.2,")"))') &
                   atm(nt), zv(nt), amass(nt)/amconv, psd(nt), xp
     end if
  enddo

  if (calc.eq.'cd' .or. calc.eq.'cm' ) &
     WRITE( stdout, '(/5x," cell mass =", f10.5, " UMA ")') cmass/amconv
  if (calc.eq.'nd' .or. calc.eq.'nm' ) &
     WRITE( stdout, '(/5x," cell mass =", f10.5, " UMA/(a.u.)^2 ")') cmass/amconv

  if (lsda) then
     WRITE( stdout, '(/5x,"Starting magnetic structure ", &
          &      /5x,"atomic species   magnetization")')
     do nt = 1, ntyp
        WRITE( stdout, '(5x,a6,9x,f6.3)') atm(nt), starting_magnetization(nt)
     enddo
  endif
  !
  !   description of symmetries
  !
  if (nsym.le.1) then
     WRITE( stdout, '(/5x,"No symmetry!")')
  else
     if (invsym) then
        WRITE( stdout, '(/5x,i2," Sym.Ops. (with inversion)",/)') nsym
     else
        WRITE( stdout, '(/5x,i2," Sym.Ops. (no inversion)",/)') nsym
     endif
  endif
  if (iverbosity.eq.1) then
     WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
     do isym = 1, nsym
        WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
        call s_axis_to_cart (s(1,1,isym), sr, at, bg)
        if (ftau(1,isym).ne.0.or.ftau(2,isym).ne.0.or.ftau(3,isym).ne.0) then
           ft1 = at(1,1)*ftau(1,isym)/nr1 + at(1,2)*ftau(2,isym)/nr2 + &
                 at(1,3)*ftau(3,isym)/nr3
           ft2 = at(2,1)*ftau(1,isym)/nr1 + at(2,2)*ftau(2,isym)/nr2 + &
                 at(2,3)*ftau(3,isym)/nr3
           ft3 = at(3,1)*ftau(1,isym)/nr1 + at(3,2)*ftau(2,isym)/nr2 + &
                 at(3,3)*ftau(3,isym)/nr3
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                 &        " )    f =( ",f10.7," )")') &
                 isym, (s(1,ipol,isym),ipol=1,3), dble(ftau(1,isym))/dble(nr1)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                       (s(2,ipol,isym),ipol=1,3), dble(ftau(2,isym))/dble(nr2)
           WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                       (s(3,ipol,isym),ipol=1,3), dble(ftau(3,isym))/dble(nr3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                 &        " )    f =( ",f10.7," )")') &
                 isym, (sr(1,ipol),ipol=1,3), ft1
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                       (sr(2,ipol),ipol=1,3), ft2
           WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                       (sr(3,ipol),ipol=1,3), ft3
        else
           WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                                     isym,  (s (1, ipol, isym) , ipol = 1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
           WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
           WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                                         isym,  (sr (1, ipol) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )")')  (sr (2, ipol) , ipol = 1, 3)
           WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol) , ipol = 1, 3)
        endif
     enddo

  endif
  !
  !    description of the atoms inside the unit cell
  !
  WRITE( stdout, '(/,3x,"Cartesian axes")')
  WRITE( stdout, '(/,5x,"site n.     atom                  positions (a_0 units)")')

  WRITE( stdout, '(7x,i3,8x,a6," tau(",i3,") = (",3f11.7,"  )")') &
             (na, atm(ityp(na)), na, (tau(ipol,na), ipol=1,3), na=1,nat)
  !
  !  output of starting magnetization
  !
  if (iverbosity.eq.1) then
     !
     !   allocate work space
     !
     allocate (xau(3,nat))
     !
     !     Compute the coordinates of each atom in the basis of the direct la
     !     vectors
     !
     do na = 1, nat
        do ipol = 1, 3
           xau(ipol,na) = bg(1,ipol)*tau(1,na) + bg(2,ipol)*tau(2,na) + &
                          bg(3,ipol)*tau(3,na)
        enddo
     enddo
     !
     !   description of the atoms inside the unit cell
     !   (in crystallographic coordinates)
     !
     WRITE( stdout, '(/,3x,"Crystallographic axes")')
     WRITE( stdout, '(/,5x,"site n.     atom        ", &
          &             "          positions (cryst. coord.)")')

     WRITE( stdout, '(7x,i2,8x,a6," tau(",i3,") = (",3f11.7,"  )")') &
           (na, atm(ityp(na)), na,  (xau(ipol,na), ipol=1,3), na=1,nat)
     !
     !   deallocate work space
     !
     deallocate(xau)
  endif

  if (lgauss) then
     WRITE( stdout, '(/5x,"number of k points=",i5, &
          &               "  gaussian broad. (ryd)=",f8.4,5x, &
          &               "ngauss = ",i3)') nkstot, degauss, ngauss
  else if (ltetra) then
     WRITE( stdout,'(/5x,"number of k points=",i5, &
          &        " (tetrahedron method)")') nkstot
  else
     WRITE( stdout, '(/5x,"number of k points=",i5)') nkstot

  endif
  WRITE( stdout, '(23x,"cart. coord. in units 2pi/a_0")')
  do ik = 1, nkstot
     WRITE( stdout, '(8x,"k(",i4,") = (",3f12.7,"), wk =",f12.7)') ik, &
          (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
  enddo
  if (iverbosity.eq.1) then
     WRITE( stdout, '(/23x,"cryst. coord.")')
     do ik = 1, nkstot
        do ipol = 1, 3
           xkg(ipol) = at(1,ipol)*xk(1,ik) + at(2,ipol)*xk(2,ik) + &
                       at(3,ipol)*xk(3,ik)
           ! xkg are the component in the crystal RL basis
        enddo
        WRITE( stdout, '(8x,"k(",i4,") = (",3f12.7,"), wk =",f12.7)') &
             ik, (xkg (ipol) , ipol = 1, 3) , wk (ik)
     enddo
  endif
  ngmtot = ngm
#ifdef __PARA
  call ireduce (1, ngmtot)
#endif
  WRITE( stdout, '(/5x,"G cutoff =",f10.4,"  (", &
       &       i7," G-vectors)","     FFT grid: (",i3, &
       &       ",",i3,",",i3,")")') gcutm, ngmtot, nr1, nr2, nr3
  if (doublegrid) then
     ngmtot = ngms
#ifdef __PARA
     call ireduce (1, ngmtot)
#endif
     WRITE( stdout, '(5x,"G cutoff =",f10.4,"  (", &
          &    i7," G-vectors)","  smooth grid: (",i3, &
          &    ",",i3,",",i3,")")') gcutms, ngmtot, nr1s, nr2s, nr3s
  endif

  if (isolve.eq.2) then
     WRITE( stdout, * )
     WRITE( stdout, '(5x,"threshold for starting DIIS:   ",f10.4)') diis_ethr_cg
     WRITE( stdout, '(5x,"reduced basis size: ",1i5)') diis_ndim
  endif

#ifdef FLUSH
  call flush (6)
#endif
  return
end subroutine summary

