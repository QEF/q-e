subroutine print_bseinfo()
!prints basic info from the BSE input file
USE io_global, ONLY : ionode
USE bse_wannier, ONLY : l_truncated_coulomb, truncation_radius, &
           numw_prod,&
           dual_bse,&
           lambda,eps,&
           l_cgrad,maxit,n_eig,eps_eig, scissor,&
           l_plotexc,plotn_min,plotn_max,r_hole,l_plotaverage,&
           spectra_e_min,spectra_e_max,spectra_broad,&
           l_restart,n_eig_start, nit_lcz,l_lanczos
implicit none

if(ionode) then
   write(*,*) 'Dimension of the polarizability basis:', numw_prod
   write(*,*) 'Scissor operator (eV)=', scissor
   if(l_truncated_coulomb) then
      write(*,*) 'Using truncated Coulomb interaction'
      write(*,*) 'Truncation Radius (a.u.)=', truncation_radius
   endif
endif

return
end subroutine
