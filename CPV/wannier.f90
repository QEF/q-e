!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!WANNIER FUNCTION DYNAMICS AND ELECTRIC FIELD
!                                - M.S
module wfparm
  !calwf : parameter to control the Wannier function calculation
  !        0:   not calculate wannier function
  !        1:   plot the electron density of wannier functions
  !       >1:   calculate the overlap matrix
  !nwf   : number of Wannier functions to be plotted
  !iplot : index of Wannier functions to be plotted
  !iwf   : index of Wannier functions
  !nw    : number of the G vectors in the general symmetry cases
  !nwrwf : file to which Wannier Function is written
  !O     : Overlap Matrix
!indexplus:
!indexminus: Indices corresponding to the translated G's

  integer :: calwf, nwf, wffort, iwf, nw, nwrwf ,jwf
  integer, allocatable :: iplot(:), wfg1(:), wfg(:,:)
  integer, allocatable :: indexplus(:,:), indexminus(:,:)
  integer, allocatable :: indexplusz(:), indexminusz(:)
  integer, allocatable :: tag(:,:), tagp(:,:)

  !weight : weights of G vectors
  real(kind=8), allocatable :: weight(:)
  real(kind=8) , allocatable :: gnx(:,:)
  integer , allocatable :: gnn(:,:)
  complex(kind=8), allocatable :: expo(:,:)
  logical :: writev

end module wfparm


module wfparm2
! nw:       the number of the G vectors
! calwf:    a parameter used to debug the code
! nit:      the number of total iteration during searching
! nsd:      the number of steepest descent iterations
! ibrav:    the structure index, the same as ibrav in CP code.
! integer :: nw, nit, nsd, ibrav
! integer :: calwf
! wfdt:     time step during searching
! maxwfdt:  the maximum time step during searching
! b1,b2,b3: the reciprocal lattice
! alat:     the lattice parameter
! a1,a2,a3: the real-space lattice
!  real(kind=8) :: wfdt, maxwfdt, b1(3), b2(3), b3(3), alat
!  real(kind=8) :: a1(3), a2(3), a3(3)
! wfg:      the G vectors involoved in general symmetry calculation
!           the units are b1, b2, b3.
!           For example:
!            the ith G vector: wfg(i,1)*b1+wfg(i,2)*b2+wfg(i,3)*b3
!  integer, allocatable :: wfg(:,:)
!
! weight:   the weight of each G vectors
!  real(kind=8), allocatable :: weight(:)
!
!       These are the Input variables for Damped Dynamics
!
! q:            fictitious mass of the Unitary Matrix
! dt:           Time Step for damped dynamics
! fric:         damping coefficient, b/w 0 and 1
! nsteps:       Max No. of MD Steps
! tolw :        Convergence Criterion
! adapt :       Automatic Adaptation of friction
  real(kind=8) :: q, dt, friction, tolw, wfdt, maxwfdt
  logical :: adapt, wfsd
  integer :: nsteps, nit, nsd

end module wfparm2

module efcalc
! switch:       whether to switch on the field adiabatically
! efield:       whether to do dynamics in the presence of an external EF
! efx0, efy0, efz0: initial magnitude of the E Field in x, y and z dir
! efx1, efy1, efz1: final magnitude of the E Field in x, y and z dir
! sw_len:       the # of steps over which field should be switched on
  logical efield, switch
  real(kind=8) :: efx1, efy1, efz1, efx, efy, efz, sw_step
  real(kind=8) :: efx0, efy0, efz0
  real(kind=8),allocatable :: xdist(:),ydist(:),zdist(:)
  integer sw_len

contains

  subroutine clear_nbeg( nbeg )

!  some more electric field stuff
!             - M.S
    IF(EFIELD) THEN
      IF(SWITCH) THEN
        WRITE(6,*) "!-------------------------------------------------------!"
        WRITE(6,*) "!                                                       !"
        WRITE(6,*) "! NBEG IS SET TO 0 FOR ADIABATIC SWITCHING OF THE FIELD !"
        WRITE(6,*) "!                                                       !"
        WRITE(6,*) "!-------------------------------------------------------!"
        nbeg=0
      END IF
    END IF
    return
  end subroutine clear_nbeg

  subroutine ef_force( fion, na, nsp, zv )
    implicit none
    real(kind=8) :: fion( :, : ), zv(:)
    integer :: na(:), nsp
    integer :: is, ia, isa
    !  Electric Feild for ions here
    !                         - M.S
     if(efield) then
        isa = 0
        do is=1,nsp
           do ia=1,na(is)
               isa = isa + 1
               fion(1,isa)=fion(1,isa)+efx*zv(is)
               fion(2,isa)=fion(2,isa)+efy*zv(is)
               fion(3,isa)=fion(3,isa)+efz*zv(is)
           end do
        end do
      end if
    !   End Electric field for ions
    !                           - M.S
    return
  end subroutine ef_force

end module efcalc


! end Wannier function and electric field module
module tune
       logical shift
       integer npts,av0,av1, xdir,ydir,zdir, start
       real(kind=8) alpha,b
end module tune
!                                     - M.S



MODULE wannier_module
  IMPLICIT NONE
  SAVE
  logical :: what1, wann_calc
  real(kind=8), allocatable :: utwf(:,:)
  real(kind=8), allocatable :: wfc(:,:)
  real(kind=8), allocatable :: rhos1(:,:), rhos2(:,:)
  complex(kind=8), allocatable :: rhogdum(:,:)
!N.B:      In the presence of an electric field every wannier state feels a different
!          potantial, which depends on the position of its center. RHOS is read in as
!          the charge density in subrouting vofrho and overwritten to be the potential.
!                                                                        -M.S
  real(kind=8) :: wfx, wfy, wfz, ionx, iony, ionz

CONTAINS

  SUBROUTINE allocate_wannier( n, nnrsx, nspin, ng )
    integer, intent(in) :: n, nnrsx, nspin, ng
    allocate(utwf(n,n))
    allocate(wfc(3,n))
    allocate(rhos1(nnrsx, nspin))
    allocate(rhos2(nnrsx, nspin))
    allocate(rhogdum(ng,nspin))
    RETURN
  END SUBROUTINE allocate_wannier
  SUBROUTINE deallocate_wannier( )
    if( allocated(utwf) ) deallocate(utwf)
    if( allocated(wfc) ) deallocate(wfc)
    if( allocated(rhos1) ) deallocate(rhos1)
    if( allocated(rhos2) ) deallocate(rhos2)
    if( allocated(rhogdum) ) deallocate(rhogdum)
    RETURN
  END SUBROUTINE deallocate_wannier

END MODULE wannier_module

MODULE electric_field_module
  IMPLICIT NONE
  SAVE
  logical :: field_tune, ft
  real(kind=8) :: efe_elec, efe_ion, prefactor, e_tuned(3)
  real(kind=8) ::  tt(3), cdmm(3), tt2(3)
  real(kind=8) :: par, alen, blen, clen, rel1(3), rel2(3)
!     ====  1 Volt / meter      = 1/(5.1412*1.e+11) a.u.            ====
END MODULE electric_field_module


MODULE wannier_subroutines

  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE read_efwan_param( nbeg )

  use wannier_module
  use electric_field_module
  use tune
  use efcalc
  use wfparm
  use wfparm2

  IMPLICIT NONE

  integer, intent(in) :: nbeg
  integer :: i

  what1=.false.
  write(6,*) wann_calc
  wann_calc=.false.
  INQUIRE (file='WANNIER', EXIST=wann_calc)
  IF(wann_calc) then
    OPEN(unit=1,file='WANNIER', status='old')
    read(1,*) efield, switch
    read(1,*) sw_len
    if(sw_len.le.0) sw_len=1
    read(1,*) efx0, efy0, efz0
    read(1,*) efx1, efy1, efz1
    read(1,*) wfsd
    read(1,*) wfdt, maxwfdt, nit, nsd
    read(1,*) q, dt, friction, nsteps
    read(1,*) tolw
    read(1,*) adapt
    read(1,*) calwf, nwf, wffort
    if(nwf.gt.0) allocate(iplot(nwf))
    do i=1,nwf
      read(1,*) iplot(i)
    end do
    read(1,*) writev
    CLOSE(1)
    if(nbeg.eq.-2.and.(efield)) then
      WRITE(6,*) "ERROR! ELECTRIC FIELD MUST BE SWITCHED ON ONLY AFTER OBTAINING THE GROUND STATE"
      WRITE(6,*) "-------------------------THE PROGRAM WILL STOP---------------------------------"
      CALL errore(' read_efwan_param ', ' electric field ', 1 )
    end if
  end if
  field_tune=.false.
  ft=.false.
  INQUIRE(file='FIELD_TUNE', EXIST=ft)
  if(ft) then
    OPEN(unit=1, file='FIELD_TUNE', status='old')
    read(1,*) field_tune
    if(field_tune) then
      efx0=0.d0
      efy0=0.d0
      efz0=0.d0
      efx1=0.d0
      efy1=0.d0
      efz1=0.d0
    end if
    read(1,*) shift, start
    read(1,*) npts, av0, av1
    read(1,*) zdir, alpha,b
    CLOSE(1)
  end if
 
END SUBROUTINE read_efwan_param

SUBROUTINE wannier_init( ibrav, alat, a1, a2, a3, b1, b2, b3 )

  use wannier_module
  use electric_field_module
  use tune
  use efcalc
  use wfparm
  use wfparm2

  IMPLICIT NONE

  integer :: ibrav
  real(kind=8) :: a1(3), a2(3), a3(3)
  real(kind=8) :: b1(3), b2(3), b3(3)
  real(kind=8) :: alat

  integer :: i

!-------------------------------------------------------------------
!More Wannier and Field Initialization
!                               - M.S
!-------------------------------------------------------------------
    if (calwf.gt.1) then
      if(calwf.eq.3) then
        write(6,*) "------------------------DYNAMICS IN THE WANNIER BASIS--------------------------"
        write(6,*) "                             DYNAMICS PARAMETERS "
        if(wfsd) then
          write(6,12132) wfdt
          write(6,12133) maxwfdt
          write(6,*) nsd,"STEPS OF STEEPEST DESCENT FOR OPTIMIZATION OF THE SPREAD"
          write(6,*) nit-nsd,"STEPS OF CONJUGATE GRADIENT FOR OPTIMIZATION OF THE SPREAD"
        else
          write(6,12125) q
          write(6,12126) dt
          write(6,12124) friction
          write(6,*) nsteps,"STEPS OF DAMPED MOLECULAR DYNAMICS FOR OPTIMIZATION OF THE SPREAD"
        end if
        write(6,*) "AVERAGE WANNIER FUNCTION SPREAD WRITTEN TO     FORT.24"
        write(6,*) "INDIVIDUAL WANNIER FUNCTION SPREAD WRITTEN TO  FORT.25"
        write(6,*) "WANNIER CENTERS WRITTEN TO                     FORT.26"
        write(6,*) "SOME PERTINENT RUN-TIME INFORMATION WRITTEN TO FORT.27"
        write(6,*) "-------------------------------------------------------------------------------"
        write(6,*)
12124   format(' DAMPING COEFFICIENT USED FOR WANNIER FUNCTION SPREAD OPTIMIZATION = ',f10.7)
12125   format(' FICTITIOUS MASS PARAMETER USED FOR SPREAD OPTIMIZATION            = ',f7.1)
12126   format(' TIME STEP USED FOR DAMPED DYNAMICS                                = ',f10.7)
!
12132   format(' SMALLEST TIMESTEP IN THE SD / CG DIRECTION FOR SPREAD OPTIMIZATION= ',f10.7)
12133   format(' LARGEST TIMESTEP IN THE SD / CG DIRECTION FOR SPREAD OPTIMIZATION = ',f10.7)
      end if
      WRITE(6,*) "IBRAV SELECTED:",ibrav

      call recips( a1, a2, a3, b1, b2, b3 )
      b1 = b1 * alat
      b2 = b2 * alat
      b3 = b3 * alat

      call wfunc_init( calwf, b1, b2, b3, ibrav)
      write (6,*) "out from wfunc_init"
      write(6,*)
      utwf=0.d0
      do i=1, SIZE( utwf, 1 )
        utwf(i, i)=1.d0
      end do
    end if
    if(efield) then
       call grid_map
       write(6,*) "GRID MAPPING DONE"
       write(6,*) "DYNAMICS IN THE PRESENCE OF AN EXTERNAL ELECTRIC FIELD"
       write(6,*)
       write(6,*) "POLARIZATION CONTRIBUTION OUTPUT TO FORT.28 IN THE FOLLOWING FORMAT"
       write(6,*)
       write(6,*) "EFX, EFY, EFZ, ELECTRIC ENTHANLPY(ELECTRONIC), ELECTRIC ENTHALPY(IONIC)"
       write(6,*)
       write(6,12121) efx0
       write(6,12122) efy0
       write(6,12123) efz0
       write(6,12128) efx1
       write(6,12129) efy1
       write(6,12130) efz1
       if(switch) then
           write(6,12127) sw_len
       end if
       write(6,*)
12121   format(' E0(x) = ',f10.7)
12122   format(' E0(y) = ',f10.7)
12123   format(' E0(z) = ',f10.7)
12128   format(' E1(x) = ',f10.7)
12129   format(' E1(y) = ',f10.7)
12130   format(' E1(z) = ',f10.7)
12131   format(' Efield Now ' ,3(f12.8,1x))
12127   format(' FIELD WILL BE TURNED ON ADIBATICALLY OVER ',i5,' STEPS')
      end if
!--------------------------------------------------------------------------
!               End of more initialization - M.S
!--------------------------------------------------------------------------


  RETURN
END SUBROUTINE wannier_init


SUBROUTINE get_wannier_center( tfirst, cm, bec, becdr, eigr, eigrb, taub, irb, ibrav, b1, b2, b3 )
  use efcalc, only: efield  
  use wfparm, only: calwf, jwf
  use wannier_module, only: what1, wfc, utwf
  IMPLICIT NONE
  logical, intent(in) :: tfirst
  complex(kind=8) :: cm(:,:,:,:)
  real(kind=8) :: bec(:,:), becdr(:,:,:)
  complex(kind=8) :: eigrb(:,:), eigr(:,:)
  integer :: irb(:,:)
  real(kind=8) :: taub(:,:)
  integer :: ibrav
  real(kind=8) :: b1(:), b2(:), b3(:)
  !--------------------------------------------------------------------------
  !Get Wannier centers for the first step if efield=true
  !--------------------------------------------------------------------------
  if(efield) then
    if(tfirst) then
      what1=.true.
      jwf=1
      call wf (calwf,cm(:,:,1,1),bec,eigr,eigrb,taub,irb,b1,b2,b3,utwf,becdr,what1,wfc,jwf,ibrav)
      write(6,*) "WFC Obtained"
      what1=.false.
    end if
  end if
  RETURN
END SUBROUTINE get_wannier_center


SUBROUTINE ef_tune( rhog, tau0 )
  use electric_field_module, only: field_tune, e_tuned
  use wannier_module, only: rhogdum
  IMPLICIT NONE
  complex(kind=8) :: rhog(:,:)
  real(kind=8) :: tau0(:,:)
!-------------------------------------------------------------------
!     Tune the Electric field               - M.S
!-------------------------------------------------------------------
  if (field_tune) then
    rhogdum=rhog
    call macroscopic_average(rhogdum,tau0,e_tuned)
  end if

  RETURN
END SUBROUTINE ef_tune


SUBROUTINE write_charge_and_exit( rhog )
  use mp, only: mp_end
  use wfparm, only: writev
  IMPLICIT NONE
  complex(kind=8) :: rhog(:,:)
!-------------------------------------------------------------------
!     Write chargedensity in g-space    - M.S
      if (writev) then
         call write_rho_g(rhog)
         call mp_end()
         stop 'write_charge_and_exit'
      end if
  RETURN
END SUBROUTINE write_charge_and_exit


SUBROUTINE wf_options( tfirst, nfi, cm, rhovan, bec, becdr, eigr, eigrb, taub, irb, &
           ibrav, b1, b2, b3, rhor, rhog, rhos, enl, ekin  )

  use efcalc, only: efield
  use wfparm, only: nwf, calwf, jwf, wffort, iplot, iwf
  use wannier_module, only: what1, wfc, utwf
  use mp, only: mp_end
  use control_flags, only: iprsta

  IMPLICIT NONE

  logical, intent(in) :: tfirst
  integer :: nfi
  complex(kind=8) :: cm(:,:,:,:)
  real(kind=8) :: bec(:,:), becdr(:,:,:)
  real(kind=8) :: rhovan(:,:,:)
  complex(kind=8) :: eigrb(:,:), eigr(:,:)
  integer :: irb(:,:)
  real(kind=8) :: taub(:,:)
  integer :: ibrav
  real(kind=8) :: b1(:), b2(:), b3(:)
  complex(kind=8) :: rhog(:,:)
  real(kind=8) :: rhor(:,:), rhos(:,:)
  real(kind=8) :: enl, ekin 


  integer :: i, j

!-------------------------------------------------------------------
! Wannier Function options            - M.S
!-------------------------------------------------------------------
    jwf=1
    if (calwf.eq.1) then
      do i=1, nwf
        iwf=iplot(i)
        j=wffort+i-1
        call rhoiofr (nfi,cm, irb, eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin,j)
        if(iprsta.gt.0) write(6,*) 'Out from rhoiofr'
        if(iprsta.gt.0) write(6,*) 
      end do
      call mp_end()
      STOP 'wf_options 1'
    end if
!---------------------------------------------------------------------
    if (calwf.eq.2) then

!     calculate the overlap matrix
!
      jwf=1
      call wf (calwf,cm(:,:,1,1),bec,eigr,eigrb,taub,irb,b1,b2,b3,utwf,becdr,what1,wfc,jwf,ibrav)

      call mp_end()
      STOP 'wf_options 2'
    end if
!---------------------------------------------------------------------
    if (calwf.eq.5) then
!
      jwf=iplot(1)
      call wf (calwf,cm(:,:,1,1),bec,eigr,eigrb,taub,irb,b1,b2,b3,utwf,becdr,what1,wfc,jwf,ibrav)

      call mp_end()
      STOP 'wf_options 5'
    end if
!----------------------------------------------------------------------
! End Wannier Function options - M.S
!
!=======================================================================

  RETURN
END SUBROUTINE wf_options


SUBROUTINE ef_potential( nfi, rhos, bec, deeq, betae, c0, cm, emadt2, emaver, verl1, verl2, c2, c3 )
  use efcalc, only: efield, efx, efy, efz, efx0, efy0, efz0, efx1, efy1, efz1, &
                    switch, sw_len, sw_step, xdist, ydist, zdist
  use electric_field_module, only: field_tune, e_tuned, par, rel1, rel2
  use wannier_module, only: rhos1, rhos2, wfc
  use smooth_grid_dimensions, only: nnrsx
  use electrons_base, only: n => nbsp, nspin, nupdwn
  use cell_base, only: ainv, a1, a2, a3
  use reciprocal_vectors, only: ng0 => gstart
  use control_flags, only: tsde
  use wave_base, only: wave_steepest, wave_verlet

  IMPLICIT NONE

  integer :: nfi
  real(kind=8) :: rhos(:,:)
  real(kind=8) :: bec(:,:)
  real(kind=8) :: deeq(:,:,:,:)
  complex(kind=8) :: betae(:,:)
  complex(kind=8) :: c0( :, : ), c2( : ), c3( : )
  complex(kind=8) :: cm( :, : )
  real(kind=8) :: emadt2(:)
  real(kind=8) :: emaver(:)
  real(kind=8) :: verl1, verl2


  integer :: i, ir
 
  !    Potential for electric field
  !                    - M.S

  if(efield) then
    if(field_tune) then
      efx=e_tuned(1)
      efy=e_tuned(2)
      efz=e_tuned(3)
      write(6,12131) efx, efy,efz
12131      format(' Efield Now ' ,3(f12.8,1x))
    else
      if(switch) then
        par=0.d0
        if(nfi.le.sw_len) then
          sw_step=1.0d0/dble(sw_len)
          par=nfi*sw_step
          if(efx1.lt.efx0) then
            efx=efx0-(efx0-efx1)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
          else
            efx=efx0+(efx1-efx0)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
          end if
          if(efy1.lt.efy0) then
            efy=efy0-(efy0-efy1)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
          else
            efy=efy0+(efy1-efy0)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
          end if
          if(efz1.lt.efz0) then
            efz=efz0-(efz0-efz1)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
          else
            efz=efz0+(efz1-efz0)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
          end if
        end if
      else
        efx=efx1
        efy=efy1
        efz=efz1
      end if
    end if
  end if
  do i=1,n,2
    if(efield) then
      rhos1=0.d0
      rhos2=0.d0
      do ir=1,nnrsx
        rel1(1)=xdist(ir)*a1(1)+ydist(ir)*a2(1)+zdist(ir)*a3(1)-wfc(1,i)
        rel1(2)=xdist(ir)*a1(2)+ydist(ir)*a2(2)+zdist(ir)*a3(2)-wfc(2,i)
        rel1(3)=xdist(ir)*a1(3)+ydist(ir)*a2(3)+zdist(ir)*a3(3)-wfc(3,i)
!  minimum image convention
        call pbc(rel1,a1,a2,a3,ainv,rel1)
        if(nspin.eq.2) then
          if(i.le.nupdwn(1)) then
            rhos1(ir,1)=rhos(ir,1)+efx*rel1(1)+efy*rel1(2)+efz*rel1(3)
          else
            rhos1(ir,2)=rhos(ir,2)+efx*rel1(1)+efy*rel1(2)+efz*rel1(3)
          end if
        else
          rhos1(ir,1)=rhos(ir,1)+efx*rel1(1)+efy*rel1(2)+efz*rel1(3)
        end if
        if(i.ne.n) then
          rel2(1)=xdist(ir)*a1(1)+ydist(ir)*a2(1)+zdist(ir)*a3(1)-wfc(1,i+1)
          rel2(2)=xdist(ir)*a1(2)+ydist(ir)*a2(2)+zdist(ir)*a3(2)-wfc(2,i+1)
          rel2(3)=xdist(ir)*a1(3)+ydist(ir)*a2(3)+zdist(ir)*a3(3)-wfc(3,i+1)
!  minimum image convention
          call pbc(rel2,a1,a2,a3,ainv,rel2)
          if(nspin.eq.2) then
            if(i+1.le.nupdwn(1)) then
              rhos2(ir,1)=rhos(ir,1)+efx*rel2(1)+efy*rel2(2)+efz*rel2(3)
            else
              rhos2(ir,2)=rhos(ir,2)+efx*rel2(1)+efy*rel2(2)+efz*rel2(3)
            end if
          else
            rhos2(ir,1)=rhos(ir,1)+efx*rel2(1)+efy*rel2(2)+efz*rel2(3)
          end if
        else
          rhos2(ir,:)=rhos1(ir,:)
        end if
      end do
      call dforce_field(bec,deeq,betae,i,c0(1,i),c0(1,i+1),c2,c3,rhos1,rhos2)
    else
      call dforce(bec,betae,i,c0(1,i),c0(1,i+1),c2,c3,rhos)
    end if
    if(tsde) then
      CALL wave_steepest( cm(:, i  ), c0(:, i  ), emadt2, c2 )
      CALL wave_steepest( cm(:, i+1), c0(:, i+1), emadt2, c3 )
    else
      CALL wave_verlet( cm(:, i  ), c0(:, i  ), verl1, verl2, emaver, c2 )
      CALL wave_verlet( cm(:, i+1), c0(:, i+1), verl1, verl2, emaver, c3 )
    endif
    if (ng0.eq.2) then
      cm(1,  i)=cmplx(real(cm(1,  i)),0.0)
      cm(1,i+1)=cmplx(real(cm(1,i+1)),0.0)
    end if
  end do

  RETURN
END SUBROUTINE ef_potential


!--------------------------------------------------------------------
!Electric Field Implementation for Electric Enthalpy
!                                              - M.S
!--------------------------------------------------------------------
SUBROUTINE ef_enthalpy( enthal, tau0 )
  use efcalc, only: efield, efx, efy, efz
  use electric_field_module, only: efe_elec, efe_ion, tt2, tt
  use wannier_module, only: wfx, wfy, wfz, ionx, iony, ionz, wfc
  use electrons_base, only: n => nbsp, f
  use cell_base, only: ainv, a1, a2, a3
  use ions_base, only: na, nsp, zv
  use io_global, only: ionode

  IMPLICIT NONE

  real(kind=8) :: enthal, tau0(:,:)
  integer :: i, is, ia, isa

  if(efield) then
    !  Electronic Contribution First
    wfx=0.d0
    wfy=0.d0
    wfz=0.d0
    efe_elec=0.d0
    do i=1,n
      tt2(1)=wfc(1,i)
      tt2(2)=wfc(2,i)
      tt2(3)=wfc(3,i)
      call pbc(tt2,a1,a2,a3,ainv,tt2)
      wfx=wfx+f(i)*tt2(1)
      wfy=wfy+f(i)*tt2(2)
      wfz=wfz+f(i)*tt2(3)
    end do
    efe_elec=efe_elec+efx*wfx+efy*wfy+efz*wfz
    !Then Ionic Contribution
    ionx=0.d0
    iony=0.d0
    ionz=0.d0
    efe_ion=0.d0
    isa = 0
    do is=1,nsp
      do ia=1,na(is)
        isa = isa + 1
        tt(1)=tau0(1,isa)
        tt(2)=tau0(2,isa)
        tt(3)=tau0(3,isa)
        call pbc(tt,a1,a2,a3,ainv,tt)
        ionx=ionx+zv(is)*tt(1)
        iony=iony+zv(is)*tt(2)
        ionz=ionz+zv(is)*tt(3)
      end do
    end do
    efe_ion=efe_ion+efx*ionx+efy*iony+efz*ionz
    if( ionode ) then
      write(28,'(f12.9,1x,f12.9,1x,f12.9,1x,f20.15,1x,f20.15)') efx, efy, efz, efe_elec,-efe_ion
    end if
  end if
  enthal=enthal+efe_elec-efe_ion

  RETURN
END SUBROUTINE ef_enthalpy


SUBROUTINE wf_closing_options( nfi, c0, cm, bec, becdr, eigr, eigrb, taub, irb, &
           ibrav, b1, b2, b3, taus, tausm, vels, velsm, acc, lambda, lambdam, xnhe0, &
           xnhem, vnhe, xnhp0, xnhpm, vnhp, nhpcl, ekincm, xnhh0, xnhhm, vnhh, velh, &
           ecut, ecutw, delt, celldm, fion, tps, mat_z, occ_f )

  use efcalc, only: efield
  use wfparm, only: nwf, calwf, jwf, wffort, iplot, iwf
  use wannier_module, only: what1, wfc, utwf
  use mp, only: mp_end
  use control_flags, only: iprsta
  use electrons_base, only: n => nbsp
  use gvecw, only: ngw
  use control_flags, only: ndw
  use cell_base, only: h, hold
  use ions_base, only: pmass
  use cvan, only: nvb
  use restart_file

  IMPLICIT NONE

  integer :: nfi
  complex(kind=8) :: c0(:,:,:,:)
  complex(kind=8) :: cm(:,:,:,:)
  real(kind=8) :: bec(:,:), becdr(:,:,:)
  complex(kind=8) :: eigrb(:,:), eigr(:,:)
  integer :: irb(:,:)
  real(kind=8) :: taub(:,:)
  integer :: ibrav
  real(kind=8) :: b1(:), b2(:), b3(:)
  real(kind=8) :: taus(:,:), tausm(:,:), vels(:,:), velsm(:,:)
  real(kind=8) :: acc(:)
  real(kind=8) :: lambda(:,:), lambdam(:,:)
  real(kind=8) :: xnhe0, xnhem, vnhe, xnhp0(:), xnhpm(:), vnhp(:), ekincm
  integer      :: nhpcl
  real(kind=8) :: velh(:,:)
  real(kind=8) :: xnhh0(:,:), xnhhm(:,:), vnhh(:,:)
  real(kind=8) :: ecut, ecutw, delt, celldm(:)
  real(kind=8) :: fion(:,:), tps
  real(kind=8) :: mat_z(:,:,:), occ_f(:)

!=============================================================
! More Wannier Function Options
!                         - M.S
!=============================================================

  if(calwf.eq.4) then
    jwf=1
    call wf(calwf,c0(:,:,1,1),bec,eigr,eigrb,taub,irb,b1,b2,b3,utwf,becdr,what1,wfc,jwf,ibrav)
    if(nvb.eq.0) then
      call wf(calwf,cm(:,:,1,1),bec,eigr,eigrb,taub,irb,b1,b2,b3,utwf,becdr,what1,wfc,jwf,ibrav)
    else
        cm(1:n,1:ngw,1,1)=c0(1:n,1:ngw,1,1)
    end if

    call writefile &
     &     ( ndw,h,hold,nfi,c0(:,:,1,1),cm(:,:,1,1),taus,tausm,vels,velsm,acc,   &
     &       lambda,lambdam,xnhe0,xnhem,vnhe,xnhp0,xnhpm,vnhp,nhpcl,ekincm,   &
     &       xnhh0,xnhhm,vnhh,velh,ecut,ecutw,delt,pmass,ibrav,celldm,fion,tps, &
     &       mat_z, occ_f )


    write(6,*) 'Wannier Functions Written to unit',ndw
    call mp_end()
    STOP 'wf_closing_options 4' 
  end if

!---------------------------------------------------------

  if(calwf.eq.3) then
!   construct overlap matrix and calculate spreads and do Localization
    jwf=1
    call wf (calwf,c0(:,:,1,1),bec,eigr,eigrb,taub,irb,b1,b2,b3,utwf,becdr,what1,wfc,jwf,ibrav)
  end if
  RETURN
END SUBROUTINE wf_closing_options


END MODULE wannier_subroutines
