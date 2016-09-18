!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE plot_xanes_dipole(a,b,xnorm,ncalcv,terminator,e1s_ry,ispectra)
  !--------------------------------------------------------------------------
  ! Calculates and plots the electric-dipole absorption K-edge cross section 
  ! as a continued fraction,
  ! from the a_i and b_i coefficients calculated for each k-point 
  !--------------------------------------------------------------------------
  USE kinds,      ONLY: DP
  USE constants,  ONLY: rytoev, fpi
  USE xspectra,   ONLY: xang_mom, xemax, xemin, xiabs, xnepoint,n_lanczos, &
                        xgamma, xonly_plot, xnitermax, xe0_ry,edge, two_edges
  !*apsi  USE uspp_param, ONLY : psd  !psd(ntypx) label for the atoms 
  USE klist,      ONLY: nkstot, & ! total number of k-points
                        nks,    & ! number of k-points per pool
                        xk,     & ! k-points coordinates
                        wk        ! k-points weight
  !USE ener,       ONLY: ef
  USE io_global,  ONLY: stdout, ionode  
  USE mp_global,  ONLY: inter_pool_comm !CG
  USE lsda_mod,   ONLY: nspin,isk
  USE mp,         ONLY: mp_sum
  USE uspp_param, ONLY: upf
  USE gamma_variable_mod, ONLY: gamma_tab, gamma_mode, &
                                gamma_value, gamma_energy, gamma_file
  USE cut_valence_green,  ONLY: cut_occ_states, cut_desmooth, &
                                memu, meml, cut_nmemu, cut_nmeml

  IMPLICIT NONE

  REAL(dp), INTENT (in) :: a(xnitermax,n_lanczos,nks),&
                           b(xnitermax,n_lanczos,nks),&
                           xnorm(n_lanczos,nks)
  REAL(dp), INTENT (in) :: e1s_ry 
  INTEGER,  INTENT (in) :: ncalcv(n_lanczos,nks), ispectra
  LOGICAL,  INTENT (in) :: terminator
  
  !... Local variables
  INTEGER  :: i,ik,n,icoord           !loops
  INTEGER  :: lmax, lanczos_i, lanczos_f
  INTEGER  :: iestart, i_lanczos
  REAL(dp) :: alpha2
  REAL(dp) :: energy,de,mod_xgamma,xemax_ryd,xemin_ryd,xgamma_ryd
  REAL(dp) :: e0 ! in Ry
  REAL(dp) :: tmp_var
  REAL(dp) :: Intensity_coord(n_lanczos,xnepoint,nspin)
  REAL(dp) :: continued_fraction 
  REAL(dp) :: paste_fermi, desmooth,t1,t2,f1,f2,df1,df2,poly(4) !CG 
  LOGICAL  :: first

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  constant and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  alpha2 = fpi/137.04

  desmooth = cut_desmooth/rytoev  ! This is in Rydberg

  e0 = xe0_ry 

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Output file for the cross section 
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF( ionode ) THEN
     if(.not.two_edges.or.(two_edges.and.ispectra.eq.1)) then
        OPEN (unit=277,file='xanes.dat',form='formatted',status='unknown')
        REWIND(277) 
        !... writes input parameters in file 277
        WRITE(277,"('# Final state angular momentum:',1x,i3)") xang_mom

        IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
           WRITE(277,"('# Broadening parameter (in eV):',1x,f8.3)") xgamma
        ELSE 
           WRITE(277,'("# Energy-dependent broadening parameter:")')
           IF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
              WRITE(277,"('# -> using gamma_file:',1x,a50)") gamma_file
           ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
              WRITE(277,"('# -> first, constant up to point (',f5.2,a1,f5.2,a)") &
                   gamma_energy(1),',',gamma_value(1),') [eV]'
              WRITE(277,"('# -> then, linear up to point (',f5.2,a1,f5.2,a)") &
                   gamma_energy(2),',',gamma_value(2),') [eV]'
              WRITE(277,"('# -> finally, constant up to xemax')")
           ENDIF
        ENDIF

        WRITE(277,"('# Absorbing atom type (xiabs):',i4)") xiabs

        IF(nspin.GT.1) THEN
           WRITE(277,"('# Energy (eV)   sigma_tot   sigma_up    sigma_down ')")
        ELSE
           WRITE(277,"('# Energy (eV)   sigma')")
        ENDIF
     endif
  ENDIF

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Converts in Rydberg most of the relevant quantities
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  xemax_ryd = xemax/rytoeV + e0
  xemin_ryd = xemin/rytoeV + e0
  xgamma_ryd= xgamma/rytoeV

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Calculates the continued fraction
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  Intensity_coord(:,:,:) = 0.d0

  de = (xemax_ryd-xemin_ryd) / REAL(xnepoint-1)

!<OB>
  if( .not. two_edges ) then
     lanczos_i = 1
     lanczos_f = n_lanczos
  else
     if( ispectra == 1 ) then
        lanczos_i = 1
        if( edge == 'L23' ) then
           lanczos_f = 2
        else if( edge == 'M45' ) then
           lanczos_f = 4
        else
           write(stdout,*) 'Output not yet programmed...'
        end if
     else
        lanczos_f = n_lanczos
        if( edge == 'L23' ) then
           lanczos_i = 3
        else if( edge == 'M45' ) then
           lanczos_i = 5
        else
           write(stdout,*) 'Output not yet programmed...'
        end if
     end if
  end if
  !<OB>
  


  IF (TRIM(gamma_mode).EQ.'constant') THEN

     IF(cut_occ_states) THEN
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        iestart=(e0-xemin_ryd)/de
        do i_lanczos = lanczos_i, lanczos_f
          DO ik=1,nks
             
             first=.true.  ! to erase the memory of paste_fermi
             !<CG>
             t1=e0-desmooth
             f1=paste_fermi(t1,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df1=paste_fermi(t1-de,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df1=(f1-df1)/de
             t2=e0+desmooth
             f2=continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),t2,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)&
                  +paste_fermi(t2,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df2=continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),t2+de,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)&
                  +paste_fermi(t2+de,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df2=(df2-f2)/de
             CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly) ! calculates interpolation polynome
  
             DO n=1,xnepoint
                energy=xemin_ryd+de*(n-1)
                IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  ! interpolation 
                   tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                   tmp_var=tmp_var*xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                   !</CG>
                ELSE
                   tmp_var=0.d0
                   IF (n>iestart) THEN
                      tmp_var=  &
                           continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),energy, &
                           xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)*  &
                           xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                   ENDIF
                   tmp_var = tmp_var + paste_fermi(energy,e0,a(1,i_lanczos,ik),&
                        b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first) &
                        *xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                ENDIF
                Intensity_coord(i_lanczos,n,isk(ik)) = Intensity_coord(i_lanczos,n,isk(ik))+tmp_var*wk(ik)
             ENDDO
          ENDDO
        end do
        DEALLOCATE(memu)
        DEALLOCATE(meml)


     ELSE
        do i_lanczos = lanczos_i, lanczos_f
          DO ik=1,nks

             DO n=1,xnepoint
                energy=xemin_ryd+de*(n-1)
                tmp_var=  &
                     continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),&
                     energy,xgamma_ryd,ncalcv(i_lanczos,ik)-1, terminator)*  &
                     xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                Intensity_coord(i_lanczos,n,isk(ik)) = Intensity_coord(i_lanczos,n,isk(ik))+tmp_var*wk(ik)
             ENDDO
          ENDDO
        end do
     ENDIF

  ELSE ! nonconstant gamma

     IF(cut_occ_states) THEN
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        iestart=(e0-xemin_ryd)/de
        do i_lanczos = lanczos_i, lanczos_f 
          DO ik=1,nks

             first=.true.  ! to erase the memory of paste_fermi
  
  
             xgamma_ryd=gamma_tab(iestart)
             t1=e0-desmooth
             f1=paste_fermi(t1,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df1=paste_fermi(t1-de,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df1=(f1-df1)/de
             t2=e0+desmooth
             f2=continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),t2,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)&
                  +paste_fermi(t2,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df2=continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),t2+de,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)&
                  +paste_fermi(t2+de,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df2=(df2-f2)/de
             CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)
  
             DO n=1,xnepoint
                energy=xemin_ryd+de*(n-1)
                xgamma_ryd=gamma_tab(n)
                IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  ! interpolation
                   tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                   tmp_var=tmp_var*xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                ELSE
                   tmp_var=0.d0
                   IF (n>iestart) THEN
                      tmp_var=  &
                           continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),&
                           energy,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)*  &
                           xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                   ENDIF
                   tmp_var = tmp_var + paste_fermi(energy,e0,a(1,i_lanczos,ik),&
                        b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first) &
                        *xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                ENDIF
                !            Intensity_tot(n)=Intensity_tot(n)+tmp_var*wk(ik)
                Intensity_coord(i_lanczos,n,isk(ik)) = Intensity_coord(i_lanczos,n,isk(ik))+tmp_var*wk(ik)
             ENDDO
          ENDDO
        end do
        DEALLOCATE(memu)
        DEALLOCATE(meml)


     ELSE
        do i_lanczos = lanczos_i, lanczos_f 
          DO ik=1,nks

             DO n=1,xnepoint
                energy=xemin_ryd+de*(n-1)
                xgamma_ryd=gamma_tab(n)
                tmp_var=  &
                     continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),&
                     energy,xgamma_ryd,ncalcv(i_lanczos,ik)-1, terminator)*  &
                     xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                Intensity_coord(i_lanczos,n,isk(ik)) = Intensity_coord(i_lanczos,n,isk(ik))+tmp_var*wk(ik)
             ENDDO
          ENDDO
       end do
     ENDIF

  ENDIF ! gamma_mode


  !... Considers the two cases of constant and non-constant broadening parameter
  !... Case 1: gamma is constant
!  IF (TRIM(gamma_mode).EQ.'constant') THEN
!
!     IF(cut_occ_states) THEN
!
!        ALLOCATE(memu(cut_nmemu,2))
!        ALLOCATE(meml(cut_nmeml,2))
!        iestart = (e0-xemin_ry)/de
!
!        DO ik = 1, nks
!           first = .true.  ! to erase the memory of paste_fermi
!           !<CG>
!           t1 = e0 - desmooth
!           f1 = paste_fermi(t1, e0, a(1,1,ik), b(1,1,ik),&
!                            xgamma_ry, ncalcv(1,ik)-1,   &
!                            terminator, first)
!           df1 = paste_fermi(t1-de, e0, a(1,1,ik), b(1,1,ik),&
!                             xgamma_ry, ncalcv(1,ik)-1,      &
!                             terminator, first)
!           df1 = (f1-df1)/de
!           t2 = e0 + desmooth
!           f2 = continued_fraction(a(1,1,ik), b(1,1,ik),         &
!                                   t2, xgamma_ry, ncalcv(1,ik)-1,&
!                                   terminator)                   &
!                + paste_fermi(t2, e0, a(1,1,ik), b(1,1,ik),      &
!                              xgamma_ry, ncalcv(1,ik)-1,         &
!                              terminator, first)
!           df2 = continued_fraction(a(1,1,ik), b(1,1,ik),             &
!                                    t2+de, xgamma_ry, ncalcv(1,ik)-1, &
!                                    terminator)                       &
!                + paste_fermi(t2+de, e0, a(1,1,ik), b(1,1,ik),        &
!                              xgamma_ry, ncalcv(1,ik)-1,              &
!                              terminator, first)
!           df2 = (df2-f2)/de
!
!           !... Calculates interpolation polynome
!           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly) 
!
!           DO n = 1, xnepoint
!              energy = xemin_ry + de*(n-1)
!              IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  
!                 ! interpolation 
!                 tmp_var = poly(1) +           &
!                           poly(2)*energy +    &
!                           poly(3)*energy**2 + &
!                           poly(4)*energy**3
!                 tmp_var = tmp_var*xnorm(1,ik)*xnorm(1,ik)
!                 !</CG>
!              ELSE
!                 tmp_var = 0.d0
!                 IF (n > iestart) &
!                   tmp_var = continued_fraction(a(1,1,ik), b(1,1,ik),       &
!                                                energy, xgamma_ry,          &
!                                                ncalcv(1,ik)-1, terminator) &
!                              *xnorm(1,ik)*xnorm(1,ik)
!                 tmp_var = tmp_var +                                     &
!                           paste_fermi(energy, e0, a(1,1,ik), b(1,1,ik), &
!                                        xgamma_ry, ncalcv(1,ik)-1,       &
!                                        terminator, first)               &
!                           *xnorm(1,ik)*xnorm(1,ik)
!              ENDIF
!              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik)) &
!                                             + tmp_var*wk(ik)
!           ENDDO
!
!        ENDDO
!
!        DEALLOCATE(memu)
!        DEALLOCATE(meml)
!
!     ELSE ! if occupied states are not cut
!        DO ik = 1, nks
!           DO n = 1, xnepoint
!              energy = xemin_ry + de*(n-1)
!              tmp_var= continued_fraction(a(1,1,ik), b(1,1,ik), &
!                                          energy, xgamma_ry,    &
!                                          ncalcv(1,ik)-1,       &
!                                          terminator)           &
!                       *xnorm(1,ik)*xnorm(1,ik)
!              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik)) &
!                                             + tmp_var*wk(ik)
!           ENDDO
!        ENDDO
!     ENDIF
!
  !... Case 2: gamma is not constant (energy-dependent)
!  ELSE 
!
!     IF(cut_occ_states) THEN
!
!        ALLOCATE(memu(cut_nmemu,2))
!        ALLOCATE(meml(cut_nmeml,2))
!        iestart=(e0-xemin_ry)/de
!        DO ik=1,nks
!           first=.true.  ! to erase the memory of paste_fermi
!           xgamma_ry = gamma_tab(iestart)
!           t1 = e0 - desmooth
!           f1 = paste_fermi(t1, e0, a(1,1,ik), b(1,1,ik), &
!                            xgamma_ry, ncalcv(1,ik)-1,    &
!                            terminator, first)
!           df1 = paste_fermi(t1-de, e0, a(1,1,ik), b(1,1,ik), &
!                             xgamma_ry, ncalcv(1,ik)-1,       &
!                             terminator, first)
!           df1 = (f1-df1)/de
!           t2 = e0 + desmooth
!           f2 = continued_fraction(a(1,1,ik), b(1,1,ik),          &
!                                   t2, xgamma_ry, ncalcv(1,ik)-1, &
!                                   terminator)                    &
!                + paste_fermi(t2, e0, a(1,1,ik), b(1,1,ik),       &
!                              xgamma_ry, ncalcv(1,ik)-1,          &
!                              terminator, first)
!           df2 = continued_fraction(a(1,1,ik), b(1,1,ik),             &
!                                    t2+de, xgamma_ry, ncalcv(1,ik)-1, &
!                                    terminator)                       &
!                 + paste_fermi(t2+de, e0, a(1,1,ik), b(1,1,ik),       &
!                               xgamma_ry, ncalcv(1,ik)-1,             &
!                               terminator, first)
!           df2 = (df2 - f2)/de
!           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)
!
!           DO n=1,xnepoint
!              energy = xemin_ry + de*(n-1)
!              xgamma_ry = gamma_tab(n)
!              IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  
!                 ! interpolation
!                 tmp_var = poly(1) + &
!                           poly(2)*energy + &
!                           poly(3)*energy**2 + &
!                           poly(4)*energy**3
!                 tmp_var = tmp_var*xnorm(1,ik)*xnorm(1,ik)
!              ELSE
!                 tmp_var=0.d0
!                 IF (n>iestart) tmp_var = &
!                                   continued_fraction(a(1,1,ik), b(1,1,ik), &
!                                                      energy, xgamma_ry,    &
!                                                      ncalcv(1,ik)-1,       &
!                                                      terminator)           &
!                                   *xnorm(1,ik)*xnorm(1,ik)
!                 tmp_var = tmp_var + &
!                           paste_fermi(energy, e0, a(1,1,ik), b(1,1,ik),&
!                                       xgamma_ry, ncalcv(1,ik)-1,       &
!                                       terminator, first)               &
!                           *xnorm(1,ik)*xnorm(1,ik)
!              ENDIF
!              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik)) &
!                                             + tmp_var*wk(ik)
!           ENDDO
!
!        ENDDO
!        DEALLOCATE(memu)
!        DEALLOCATE(meml)
!
!     ELSE ! if occupied states are not cut
!        DO ik=1,nks
!           DO n=1,xnepoint
!              energy = xemin_ry + de*(n-1)
!              xgamma_ry = gamma_tab(n)
!              tmp_var =  continued_fraction(a(1,1,ik), b(1,1,ik), &
!                                            energy, xgamma_ry,    &
!                                            ncalcv(1,ik)-1,       &
!                                            terminator)           &
!                         *xnorm(1,ik)*xnorm(1,ik)
!              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))&
!                                             + tmp_var*wk(ik)
!           ENDDO
!        ENDDO
!     ENDIF
!
!  ENDIF ! gamma_mode

  !  CALL poolreduce( nspin*xnepoint, Intensity_coord )

  !<CG>  replaces poolreduce
#if defined(__MPI)
  CALL mp_sum ( Intensity_coord, inter_pool_comm )
#endif
  !</CG>

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Writes the final cross section in file 277
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(ionode) THEN
     IF(nspin == 1) THEN
        DO n=1,xnepoint
           energy = xemin_ryd + de*(n-1)
           Intensity_coord(:,n,:) = Intensity_coord(:,n,:) * &
                                    (energy+e1s_ry) *          &
                                    alpha2
           WRITE(277,'(2f14.8)') (energy-e0)*rytoeV, sum( Intensity_coord(lanczos_i:lanczos_f,n,1) )
        ENDDO
     ELSEIF(nspin == 2) THEN
        DO n=1,xnepoint
           energy = xemin_ryd + de*(n-1)
           Intensity_coord(:,n,:) = Intensity_coord(:,n,:) * &
                                    (energy+e1s_ry)          * &
                                     alpha2 !

           WRITE(277,'(4f14.8)') (energy-e0)*rytoev, &
                sum( Intensity_coord(lanczos_i:lanczos_f,n,1) )+ sum( Intensity_coord(lanczos_i:lanczos_f,n,2) ),& 
                sum( Intensity_coord(lanczos_i:lanczos_f,n,1) ), sum( Intensity_coord(lanczos_i:lanczos_f,n,2) )
        ENDDO
     ENDIF

     CLOSE(277)
  ENDIF

END SUBROUTINE plot_xanes_dipole

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE plot_xanes_quadrupole(a,b,xnorm,ncalcv,terminator,e1s_ry)
  !--------------------------------------------------------------------------
  ! Calculates and plots the electric-quadrupole absorption cross section 
  ! as a continued fraction,
  ! from the a_i and b_i coefficients calculated for each k-point 
  !--------------------------------------------------------------------------
  USE kinds,      ONLY: DP
  USE constants,  ONLY: pi, rytoev
  USE xspectra,   ONLY: xang_mom, xemax, xemin, xiabs, xnepoint, &
                        xgamma, xonly_plot, xnitermax, xe0_ry
  !*apsi  USE uspp_param, ONLY : psd  !psd(ntypx) label for the atoms 
  USE klist,      ONLY: nkstot,& ! total number of k-points
                        nks,   & ! number of k-points per pool
                        xk,    & ! k-points coordinates
                        wk       ! k-points weight
  !USE ener,       ONLY: ef
  USE io_global,  ONLY: stdout,ionode  
  USE mp_global,  ONLY: inter_pool_comm
  USE lsda_mod,   ONLY: nspin,isk
  USE mp,         ONLY: mp_sum
  USE uspp_param, ONLY: upf
  USE gamma_variable_mod, ONLY : gamma_tab, gamma_mode, &
                                 gamma_file, gamma_value, gamma_energy
  USE cut_valence_green , ONLY : cut_occ_states, cut_desmooth, &
                                 cut_nmemu, cut_nmeml, memu, meml

  IMPLICIT NONE

  REAL(dp), INTENT (in) :: a(xnitermax,1,nks),&
                           b(xnitermax,1,nks),&
                           xnorm(1,nks)
  REAL(dp), INTENT (in) :: e1s_ry 
  INTEGER,  INTENT (in) :: ncalcv(1,nks)
  LOGICAL,  INTENT (in) :: terminator

  !... Local variables
  INTEGER  :: i,ik,n,icoord           !loops
  INTEGER  :: lmax
  INTEGER  :: iestart
  REAL(dp) :: alpha2,constantqua
  REAL(dp) :: energy,de,mod_xgamma,xemax_ry,xemin_ry,xgamma_ry
  REAL(dp) :: e0
  REAL(dp) :: tmp_var
  REAL(dp) :: Intensity_tot(xnepoint,nspin)
  REAL(dp) :: continued_fraction
  REAL(dp) :: paste_fermi,desmooth,t1,t2,f1,f2,df1,df2,poly(4)
  LOGICAL  :: first

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  constant and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  constantqua = pi/(137.04*137.04*137.04)

  alpha2=4.d0*pi/137.04

  desmooth=cut_desmooth/rytoev  ! This is in Rydberg
  
  e0 = xe0_ry 

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Output file for the cross section 
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF( ionode ) THEN

     open (unit=277,file='xanes.dat',form='formatted',status='unknown')

     REWIND(277) ! to be at the initial point of file 277
     !... writes input parameters in file 277
     WRITE(277,"('# final state angular momentum:',1x,i3)") xang_mom

     IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
        WRITE(277,"('# Broadening parameter (in eV):',1x,f8.3)") xgamma
     ELSE
        WRITE(277,'("# Energy-dependent broadening parameter:")')
        IF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
           WRITE(277,"('# -> using gamma_file:',1x,a50)") gamma_file
        ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
           WRITE(277,"('# -> first, constant up to point (',f5.2,a1,f5.2,a)") &
                 gamma_energy(1),',',gamma_value(1),') [eV]'
           WRITE(277,"('# -> then, linear up to point (',f5.2,a1,f5.2,a)") &
                 gamma_energy(2),',',gamma_value(2),') [eV]'
           WRITE(277,"('# -> finally, constant up to xemax')")
        ENDIF
     ENDIF

     WRITE(277,"('# Absorbing atom type (xiabs):',i4)") xiabs

     IF(nspin.GT.1) THEN
        WRITE(277,"('# Energy (eV)   sigma_tot   sigma_up    sigma_down ')")
     ELSE
        WRITE(277,"('# Energy (eV)   sigma')")
     ENDIF

  ENDIF

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Converts in Rydberg most of the relevant quantities
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  xemax_ry = xemax/rytoev + e0
  xemin_ry = xemin/rytoev + e0
  xgamma_ry= xgamma/rytoev

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Calculates the continued fraction
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  Intensity_tot(:,:)=0.d0

  de = (xemax_ry-xemin_ry) / REAL(xnepoint-1)

  !... Considers the two cases of constant and non-constant broadening parameter
  !... Case 1: gamma is constant
  IF (TRIM(gamma_mode).EQ.'constant') THEN

     IF(cut_occ_states) THEN

        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        DO ik=1,nks
           iestart = (e0 - xemin_ry)/de
           first = .true.
           t1 = e0 - desmooth
           f1 = paste_fermi(t1, e0, a(1,1,ik), b(1,1,ik), &
                            xgamma_ry, ncalcv(1,ik)-1,    &
                            terminator, first)
           df1 = paste_fermi(t1-de, e0, a(1,1,ik), b(1,1,ik), &
                             xgamma_ry, ncalcv(1,ik)-1,       &
                             terminator, first)
           df1 = (f1 - df1)/de
           t2 = e0 + desmooth
           f2 = continued_fraction(a(1,1,ik), b(1,1,ik), t2,  &
                                   xgamma_ry, ncalcv(1,ik)-1, &
                                   terminator)                &
                + paste_fermi(t2, e0, a(1,1,ik), b(1,1,ik),   &
                              xgamma_ry, ncalcv(1,ik)-1,      &
                              terminator, first)
           df2 = continued_fraction(a(1,1,ik), b(1,1,ik), t2+de, &
                                    xgamma_ry, ncalcv(1,ik)-1,   &
                                    terminator)                  &
                 + paste_fermi(t2+de, e0, a(1,1,ik), b(1,1,ik),  &
                               xgamma_ry,ncalcv(1,ik)-1,terminator, first)
           df2=(df2-f2)/de
           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           DO n=1,xnepoint
              energy = xemin_ry + de*(n-1)
              IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  
              !interpolation
                 tmp_var = poly(1) + &
                           poly(2)*energy + &
                           poly(3)*energy**2 + &
                           poly(4)*energy**3
                 tmp_var = tmp_var*xnorm(1,ik)*xnorm(1,ik)
              ELSE
                 tmp_var=0.d0
                 IF (n>iestart) tmp_var = &
                                   continued_fraction(a(1,1,ik), b(1,1,ik),&
                                                      energy, xgamma_ry,   &
                                                      ncalcv(1,ik)-1,      &
                                                      terminator)          &
                                   *xnorm(1,ik)*xnorm(1,ik)
                 tmp_var = tmp_var + &
                           paste_fermi(energy, e0, a(1,1,ik), b(1,1,ik), &
                                       xgamma_ry, ncalcv(1,ik)-1,        &
                                       terminator, first)                &
                           *xnorm(1,ik)*xnorm(1,ik)
              ENDIF
              Intensity_tot(n,isk(ik)) = Intensity_tot(n,isk(ik)) &
                                         + tmp_var*wk(ik)
           ENDDO
        ENDDO
        DEALLOCATE(memu)
        DEALLOCATE(meml)

     ELSE ! if occupied states are not cut
        DO ik=1,nks
           DO n=1,xnepoint
              energy = xemin_ry + de*(n-1)
              tmp_var= continued_fraction(a(1,1,ik), b(1,1,ik), &
                                          energy, xgamma_ry,    &
                                          ncalcv(1,ik)-1,       &
                                          terminator)           &
                       *xnorm(1,ik)*xnorm(1,ik)
              Intensity_tot(n,isk(ik)) = Intensity_tot(n,isk(ik)) &
                                         + tmp_var*wk(ik)
           ENDDO
        ENDDO
     ENDIF
 
  !... Case 2: gamma is not constant (energy-dependent)
  ELSE 

     IF(cut_occ_states) THEN ! if occupied states are cut
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        iestart = (e0-xemin_ry)/de
        DO ik=1,nks
           first = .true. ! to erase memory of paste_fermi
           xgamma_ry = gamma_tab(iestart)
           t1 = e0 - desmooth
           f1 = paste_fermi(t1, e0, a(1,1,ik), b(1,1,ik), &
                            xgamma_ry, ncalcv(1,ik)-1,    &
                            terminator, first)
           df1 = paste_fermi(t1-de, e0, a(1,1,ik), b(1,1,ik), &
                             xgamma_ry, ncalcv(1,ik)-1,       &
                             terminator, first)
           df1 = (f1-df1) / de
           t2 = e0 + desmooth
           f2 = continued_fraction(a(1,1,ik), b(1,1,ik), t2,  &
                                   xgamma_ry, ncalcv(1,ik)-1, &
                                   terminator)                &
                + paste_fermi(t2, e0, a(1,1,ik), b(1,1,ik),   &
                              xgamma_ry, ncalcv(1,ik)-1,      &
                              terminator, first)
           df2 = continued_fraction(a(1,1,ik), b(1,1,ik), t2+de, &
                                    xgamma_ry, ncalcv(1,ik)-1,   &
                                    terminator)                  &
                 + paste_fermi(t2+de, e0, a(1,1,ik), b(1,1,ik),  &
                               xgamma_ry, ncalcv(1,ik)-1,        &
                               terminator, first)
           df2 = (df2-f2) / de
           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           DO n=1,xnepoint
              energy = xemin_ry + de*(n-1)
              xgamma_ry = gamma_tab(n)
              IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  
                 ! interpolation
                 tmp_var = poly(1) + &
                           poly(2)*energy + &
                           poly(3)*energy**2 + &
                           poly(4)*energy**3
                 tmp_var = tmp_var*xnorm(1,ik)*xnorm(1,ik)
              ELSE
                 tmp_var = 0.d0
                 IF (n>iestart) tmp_var = &
                                  continued_fraction(a(1,1,ik), b(1,1,ik), &
                                                     energy, xgamma_ry,    &
                                                     ncalcv(1,ik)-1,       &
                                                     terminator)           &
                                  *xnorm(1,ik)*xnorm(1,ik)
                 tmp_var = tmp_var + &
                           paste_fermi(energy, e0, a(1,1,ik),b(1,1,ik),&
                                       xgamma_ry, ncalcv(1,ik)-1,      &
                                       terminator, first)              &
                           *xnorm(1,ik)*xnorm(1,ik)
              ENDIF
              Intensity_tot(n,isk(ik)) = Intensity_tot(n,isk(ik)) &
                                         + tmp_var*wk(ik)
           ENDDO
        ENDDO
        DEALLOCATE(memu)
        DEALLOCATE(meml)

     ELSE ! if occupied states are not cut
        DO ik=1,nks
           DO n=1,xnepoint
              energy = xemin_ry + de*(n-1)
              xgamma_ry = gamma_tab(n)
              tmp_var = continued_fraction(a(1,1,ik), b(1,1,ik), &
                                           energy, xgamma_ry,    &
                                           ncalcv(1,ik)-1,       &
                                           terminator)           &
                        *xnorm(1,ik)*xnorm(1,ik)
              Intensity_tot(n,isk(ik)) = Intensity_tot(n,isk(ik)) &
                                         + tmp_var*wk(ik)
           ENDDO
        ENDDO
     ENDIF

  ENDIF ! gammma_mode

  !  CALL poolreduce( nspin*xnepoint, Intensity_tot )

  !<CG>  replaces poolreduce
#if defined(__MPI)
  CALL mp_sum ( Intensity_tot, inter_pool_comm )
#endif
  !</CG>

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Writes the final cross section in file 277
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(ionode) THEN
     IF(nspin == 1) THEN
        DO n=1,xnepoint
           energy = xemin_ry + de*(n-1)
           !Intensity_tot(n,:)=Intensity_tot(n,:) &
           !                   *(energy+e_1s)*(energy+e_1s)*(energy+e_1s) &
           !                   * constantqua      !normalized
           Intensity_tot(n,:) = Intensity_tot(n,:) &
                                * (energy+e1s_ry)**3 &
                                * constantqua     
           WRITE(277,'(2f14.8)') (energy-e0)*rytoev, Intensity_tot(n,:)
        ENDDO
     ELSEIF(nspin == 2) THEN
        DO n=1,xnepoint
           energy = xemin_ry + de*(n-1)
           !Intensity_tot(n,:) = Intensity_tot(n,:) &
           !                     *(energy+e_1s)*(energy+e_1s)*(energy+e_1s) &
           !                     *constantqua     !normalized
           Intensity_tot(n,:) = Intensity_tot(n,:) &
                                * (energy+e1s_ry)**3 &
                                * constantqua     
           WRITE(277,'(4f14.8)') (energy-e0)*rytoev, &
                                  Intensity_tot(n,1)+Intensity_tot(n,2),&
                                  Intensity_tot(n,1), Intensity_tot(n,2)
        ENDDO
     ENDIF

     CLOSE(277)
  ENDIF
  

END SUBROUTINE plot_xanes_quadrupole

