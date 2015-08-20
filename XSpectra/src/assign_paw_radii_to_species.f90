subroutine assign_paw_radii_to_species(rc,r_paw)
  USE kinds, ONLY : DP
  USE io_global,       ONLY : stdout
  USE parameters,      ONLY : ntypx,lmaxx
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp
  USE paw_gipaw,       ONLY : paw_recon
  USE xspectra,        ONLY : xiabs
  implicit none
!  internal
  integer :: nt,j,il
  REAL (DP) :: rc(ntypx,0:lmaxx),r_paw(0:lmaxx)
  
  WRITE(stdout,1000) ! return+line 
  WRITE(stdout,'(5x,a)')  &
       '                         Attributing the PAW radii '
  WRITE(stdout,'(5x,a)')  &
       '                for the absorbing atom [units: Bohr radius]'
  WRITE(stdout,1001) ! line+return
  
  DO nt=1,ntyp
     IF ((.NOT.paw_recon(nt)%gipaw_data_in_upf_file)) then
        paw_recon(nt)%paw_nbeta=0
     ENDIF
     DO  j=1,paw_recon(nt)%paw_nbeta
        il=paw_recon(nt)%psphi(j)%label%l
        IF(xiabs.EQ.nt.AND.DABS(r_paw(il)).lt.1.d-6) THEN
           !*apsi  if(r(kkbeta(nt),nt).GT.1.d-3) THEN
           !*apsi     rc(nt,il)=  r(kkbeta(nt),nt)
           !*apsi  ELSE
           !*apsi     WRITE(stdout,*) 'Warning-no-kkbeta r_paw(',il,')=1.0'
           !*apsi     rc(nt,il)=1.0
           !*apsi  ENDIF
           !<CG>  to be verified
           IF (paw_recon(nt)%psphi(j)%label%rc > 1.d-3) THEN
              rc(nt, il)= paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
              WRITE(stdout,'(8x,a,i2,a,i2,a,f5.2,a)') &
                   'PAW proj', j,': r_paw(l=',il,')=',rc(nt,il),'  (1.5*r_cut)'
           ELSE
              rc(nt, il)= 1.5d0
              WRITE(stdout,'(8x,a,i2,a,i2,a,f5.2,a)') &
                   'PAW proj',j,': r_paw(l=',il,')=',rc(nt,il),'  (set to 1.5)'
           ENDIF
        ELSEIF(xiabs.EQ.nt.AND.DABS(r_paw(il)).GT.1.d-6) THEN
           rc(nt,il)=r_paw(il)
        ELSEIF(nt.NE.xiabs) THEN
           !*apsi  IF(r(kkbeta(nt),nt).GT.1.d-3) THEN
           !*apsi     rc(nt,il)=r(kkbeta(nt),nt)
           !*apsi  ELSE
           !*apsi     rc(nt,il)=1.0
           !*apsi  ENDIF
           !<CG> to be verified
           IF(paw_recon(nt)%psphi(j)%label%rc.GT.1.d-3) THEN
              rc(nt,il)=paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
           ELSE
              rc(nt,il)=1.5
           ENDIF
           !</CG> 
        ENDIF
        
     ENDDO
  ENDDO
  WRITE(stdout,*)
  WRITE(stdout,'(8x,a)')&
       'NB: The calculation will not necessary use all these r_paw values.'
  WRITE(stdout,'(8x,a)')&
       '    - For a edge in the electric-dipole approximation,'
  WRITE(stdout,'(8x,a)')&
       '      only the r_paw(l=1) values are used.'
  WRITE(stdout,'(8x,a)')&
       '    - For a K edge in the electric-quadrupole approximation,'
  WRITE(stdout,'(8x,a,/)')'      only the r_paw(l=2) values are used.'
  WRITE(stdout,'(8x,a,/)')     '    - For a L2 or L3 edge in the electric-quadrupole approximation,'
  WRITE(stdout,'(8x,a,/)')'      all projectors (s, p and d) are used.'
  
  !<CG>
  DO nt=1,ntyp
     DO j = 1,paw_recon(nt)%paw_nbeta
        paw_recon(nt)%psphi(j)%label%rc=rc(nt,paw_recon(nt)%psphi(j)%label%l)
        paw_recon(nt)%aephi(j)%label%rc=rc(nt,paw_recon(nt)%aephi(j)%label%l)
     ENDDO
  ENDDO
  !</CG>
  

 ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Formats 
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 1000 FORMAT(/,5x,&
  '-------------------------------------------------------------------------')
 1001 FORMAT(5x,&
  '-------------------------------------------------------------------------',&
  /)
end subroutine assign_paw_radii_to_species
