module ionic_mod
USE kinds, ONLY: DP
      use mp, only: mp_sum
      use constants, only: pi

      type ionic_init_type
         real(DP) :: I_primo, eta
         integer  :: n_max !cutoff per somme in griglia reale
         real(dp), allocatable :: I_uno_g(:, :, :), I_due_g(:)
      end type

contains

   subroutine init_ionic(init_data, eta, n_max, ngm, gstart, at, alat, omega, gg, g, tpiba2)
      use mp_pools, only: intra_pool_comm
   implicit none
   type(ionic_init_type), intent(inout) :: init_data
   integer, intent(in) :: ngm, gstart, n_max
   real(dp), intent(in) :: at(3,3), alat, eta, omega, gg(:), g(:,:), tpiba2

   real(DP) :: n(3), modul, erf_value
   real(kind=DP), external :: qe_erfc
   real(DP) ::I_primo_rec
   integer  ::  a, b, igm, n_x, n_y, n_z
           
   allocate (init_data%I_uno_g(ngm, 3, 3))
   allocate (init_data%I_due_g(ngm))
   init_data%n_max=n_max
   init_data%eta=eta

   ! I_primo == S_{2,B} (Aris note)
   init_data%I_primo = 0.d0
! these sum over n_x n_y n_z are sum_{L!=0}erfc(sqrt(eta*L)*1/L)
   do n_x = -n_max, n_max
      do n_y = -n_max, n_max
         do n_z = -n_max, n_max
            if ((n_x /= 0) .or. (n_y /= 0) .or. (n_z /= 0)) then
               n(1:3) = n_x*at(1:3, 1)*alat + n_y*at(1:3, 2)*alat + n_z*at(1:3, 3)*alat
               modul = modulus(n(:))
               erf_value = qe_erfc(sqrt(eta)*modul)
               init_data%I_primo = init_data%I_primo + erf_value/modul
            end if
         end do
      end do
   end do
!now compute the other addendum of S_{2,B}
   init_data%I_primo = init_data%I_primo - 2.d0*sqrt(eta/pi)
!
   I_primo_rec = 0.d0
   do igm = gstart, ngm
      I_primo_rec = I_primo_rec + 2.d0*(4.d0*pi)/omega*alpha_0_lr(eta, gg(igm)*tpiba2, 1)
   end do
   call mp_sum(I_primo_rec, intra_pool_comm)
   I_primo_rec = I_primo_rec - pi/(eta*omega)
!
   init_data%I_primo = init_data%I_primo + I_primo_rec
!
!calcolo di I_due_g
!
   do igm = gstart, ngm
!
      init_data%I_due_g(igm) = (4.d0*pi)/omega*alpha_0_lr(eta, gg(igm)*tpiba2, 1)
!
   end do
!
   if (gstart == 2) then
      init_data%I_due_g(1) = -pi/(eta*omega)
   end if
!
!analogamente a prima inizializzo I_uno_g
   do a = 1, 3
      do b = 1, 3
         if (a >= b) then
            do igm = gstart, ngm
               init_data%I_uno_g(igm, a, b) = (4.d0*pi)/omega*exp(-tpiba2*gg(igm)/(4.d0*eta))/(tpiba2*gg(igm))*&
 & (2 + gg(igm)*tpiba2/(2.d0*eta))*g(a, igm)*g(b, igm)/gg(igm)
            end do
!
            if (gstart == 2) then
               init_data%I_uno_g(1, a, b) = 0.d0
            end if
            if (a /= b) &
                init_data%I_uno_g(:, b, a) = init_data%I_uno_g(:, a, b)
         end if
      end do
   end do


   end subroutine

   real(DP) function alpha_0_lr(eta, G_sq, flag)
      real(DP), intent(in)  ::eta, G_sq
      integer, intent(in) ::flag
      if (flag == 1) then
         alpha_0_lr = exp(-G_sq/(4.d0*eta))/G_sq
      else
         alpha_0_lr = 1.d0/G_sq
      end if
   end function alpha_0_lr

   real(kind=DP) function h_(x)
      real(kind=DP), intent(in) :: x
      real(kind=DP), external :: qe_erfc
      h_ = qe_erfc(x)/x
   end function h_

   real(kind=DP) function hp_(x)
      real(kind=DP), intent(in) :: x
      real(kind=DP), external :: qe_erfc
      hp_ = -(2.d0/sqrt(pi))*(1.d0/x)*exp(-x*x) - 1.d0/(x*x)*qe_erfc(x)
   end function hp_

   real(kind=DP) function modulus(vector)
      real(kind=DP), intent(in) ::vector(3)
      modulus = sqrt(vector(1)**2 + vector(2)**2 + vector(3)**2)
   end function modulus

   subroutine pbc(vin, vout, alat, at, bg)
      implicit none
      real(DP), intent(in) :: vin(3), alat, at(3,3), bg(3,3)
      real(DP), intent(out) :: vout(3)
      !local
      real(DP):: alatdir ! box side in a given direction
      integer :: i, n
      !vout used also as help variable
      vout(:) = matmul(vin(:)/alat, bg(:, :)) ! vout = vin in S space
      vout(:) = vout(:) - ANINT(vout(:))
      vout(:) = MATMUL(at(:, :), vout(:))*alat ! back from S space to R space
   end subroutine pbc

   subroutine I_uno_value(init_data, y, x, a, b, flag, tpiba, at, alat, gstart, g, gg, npw)
      use mp_pools, only: intra_pool_comm
      implicit none
   type(ionic_init_type), intent(in) :: init_data
      integer, intent(in) ::flag, npw, gstart
      real(DP), intent(out) ::y
      real(DP), intent(in)  ::x(3), tpiba, g(:,:), gg(:), at(3,3),alat
      integer, intent(in)   ::a, b
      integer ::igm
      real(DP) :: comp_iso

      call start_clock('rec')

      y = 0.d0
      do igm = gstart, npw
         if ((gg(igm)*tpiba*tpiba)/(4.d0*init_data%eta) > 20.d0) exit
         y = y + 2.d0*(init_data%I_uno_g(igm, a, b)*cos(DOT_PRODUCT(g(1:3, igm), x(1:3))*tpiba))
      end do
      if (gstart == 2) y = y + init_data%I_uno_g(1, a, b)
      call mp_sum(y, intra_pool_comm)
      call stop_clock('rec')
      if (flag == 1) then
         call add_local_uno(y, init_data%n_max, init_data%eta, x, a, b, at, alat)
      end if
!change for opt
!  call I_due_value(comp_iso,x,1)
!  if (a==b) then
!     y=y-comp_iso
!  end if
      !if (a == b) then
      !   call I_due_value(comp_iso, x, 1)
      !   y = y - comp_iso
      !end if

   end subroutine I_uno_value

   subroutine add_local_uno(value, n_max, eta, pos, a, b, at, alat)
      use mp_world, only: nproc, mpime
      use mp, only: mp_sum
      use mp_pools, only: intra_pool_comm
      implicit none
      integer, intent(in) ::a, b, n_max
      real(DP), intent(in) ::pos(3), at(3,3), alat, eta
      real(DP), intent(inout) ::value
      real(DP) ::u(3), u_mod, n(3)
      integer  ::n_x, n_y, n_z
      integer   :: l_blk, nbegin, nend
      real(DP)  :: value1

      call start_clock('real')
      l_blk = (2*n_max + 1)/nproc
      if (l_blk*nproc < (2*n_max + 1)) l_blk = l_blk + 1
      nbegin = mpime*l_blk - n_max
      nend = nbegin + l_blk - 1
      if (nend > n_max) nend = n_max
      value1 = 0.d0
      do n_x = nbegin, nend
         do n_y = -n_max, n_max
            do n_z = -n_max, n_max
               n(1:3) = n_x*at(1:3, 1)*alat + n_y*at(1:3, 2)*alat + n_z*at(1:3, 3)*alat
               u(1:3) = pos(1:3) - n(1:3)
               u_mod = modulus(u)
               value1 = value1 + eta*hp_(sqrt(eta)*u_mod)*u(a)*u(b)/u_mod
               if (a == b) then
                  value1 = value1 + sqrt(eta)*h_(sqrt(eta)*u_mod)
               end if
            end do
         end do
      end do
      call mp_sum(value1, intra_pool_comm)
      value = value + value1
      call stop_clock('real')
!!!!!!!!!!!!!!!!!!!
      !call start_clock( 'real' )
      !do n_x=-n_max,n_max
      !   do n_y=-n_max,n_max
      !      do n_z=-n_max,n_max
      !         n(1:3) = n_x * at(1:3,1)*alat + n_y * at(1:3,2)*alat + n_z * at(1:3,3)*alat
      !         u(1:3) = pos(1:3)-n(1:3)
      !         u_mod=modulus(u)
      !         value=value+eta*hp_(sqrt(eta)*u_mod)*u(a)*u(b)/u_mod
      !         if (a==b) then
      !            value=value+sqrt(eta)*h_(sqrt(eta)*u_mod)
      !         end if
      !      end do
      !   end do
      !end do
      !call stop_clock( 'real' )
   end subroutine add_local_uno

   subroutine I_due_value(init_data, y, x, flag, tpiba, at, alat, gstart, g, gg, npw)
      use mp, only: mp_sum
      use mp_pools, only: intra_pool_comm
      implicit none
   type(ionic_init_type), intent(in) :: init_data
      integer, intent(in) ::flag, npw, gstart
      real(DP), intent(out) ::y
      real(DP), intent(in)  ::x(3), at(3,3), alat,g(:,:), gg(:), tpiba
      integer ::igm!,ng_max
      real(DP) ::scalar

      call start_clock('rec')
      y = 0.d0
      do igm = gstart, npw
         if ((gg(igm)*tpiba*tpiba)/(4.d0*init_data%eta) > 20.d0) exit
         y = y + 2.d0*(init_data%I_due_g(igm)*cos(DOT_PRODUCT(g(1:3, igm), x(1:3))*tpiba))
      end do
      if (gstart == 2) y = y + init_data%I_due_g(1)
      call mp_sum(y, intra_pool_comm)
      if (flag == 1) then
         call add_local_due(y, init_data%n_max, init_data%eta, x, at, alat)
      end if
      call stop_clock('rec')
   end subroutine i_due_value

   subroutine add_local_due(value, n_max, eta, pos, at, alat)
      !This routines computes sum_{L}erfc(sqrt(eta)abs(d-L))/abs(d-L)
      ! the sum is in real space over a small number of shells
      use mp, only: mp_sum
      use mp_world, ONLY: nproc, mpime
      use mp_pools, ONLY: intra_pool_comm
      implicit none
      real(DP), intent(in) ::pos(3), at(3,3), alat, eta
      real(DP), intent(inout) ::value
      integer, intent(in) :: n_max
      real(DP) ::modul, n(3), erf_value
      integer  ::n_x, n_y, n_z
      real(kind=DP), external :: qe_erfc
      integer   :: l_blk, nbegin, nend
      real(DP)  :: value1

      call start_clock('real')
      l_blk = (2*n_max + 1)/nproc
      if (l_blk*nproc < (2*n_max + 1)) l_blk = l_blk + 1
      nbegin = mpime*l_blk - n_max
      nend = nbegin + l_blk - 1
      if (nend > n_max) nend = n_max
      value1 = 0.d0
      do n_x = nbegin, nend
         do n_y = -n_max, n_max
            do n_z = -n_max, n_max
               n(1:3) = n_x*at(1:3, 1)*alat + n_y*at(1:3, 2)*alat + n_z*at(1:3, 3)*alat
               modul = modulus(pos(1:3) - n(1:3))
               erf_value = qe_erfc(sqrt(eta)*modul)
               value1 = value1 + erf_value/modul
            end do
         end do
      end do
      call mp_sum(value1, intra_pool_comm)
      value = value + value1
      call stop_clock('real')
!!!!!!!!!!
!  call start_clock( 'real' )
!  do n_x=-n_max,n_max         !primo ciclo   su n
!     do n_y=-n_max,n_max      !secondo ciclo su n
!        do n_z=-n_max,n_max   !terzo ciclo   su n
!           n(1:3) = n_x * at(1:3,1)*alat + n_y * at(1:3,2)*alat + n_z * at(1:3,3)*alat
!           modul=modulus(pos(1:3)-n(1:3))
!           erf_value=qe_erfc(sqrt(eta)*modul)
!           value=value+ erf_value/modul
!        end do
!     end do
!  end do
!  call stop_clock( 'real' )

   end subroutine add_local_due

   subroutine current_ionic(init_data, current, current_a, current_b, current_c, current_d, current_e, add_current_b, &
                         nat, tau, vel, zv, ityp, alat, at, bg, tpiba, gstart, g, gg, npw, amass)
   use kinds, only: DP
   use constants, only: AMU_RY, e2
   use io_global, only: ionode

   implicit none

   type(ionic_init_type), intent(in) :: init_data
   real(dp), intent(out) :: current(3), current_a(3), current_b(3), current_c(3),&
                            current_d(3), current_e(3)
   logical, intent(in) :: add_current_b
   INTEGER, intent(in) :: npw, nat, ityp(:), gstart
   REAL(DP), intent(in) :: zv(:), tau(:,:), vel(:,:), &
                           tpiba, alat, at(:,:), &
                           g(:,:), gg(:), bg(:,:), amass(:)

   real(dp) :: u(3), u_pbc(3), val
   integer :: iatom, jatom, a, b

   current = 0.d0
   current_a = 0.d0
   current_b = 0.d0
   current_c = 0.d0
   current_d = 0.d0
   current_e = 0.d0

   call start_clock('calcolo_i')
   do iatom = 1, nat
      current_a(:) = current_a(:) + vel(:, iatom)*(1./2.*amu_ry*amass(ityp(iatom))*(vel(1, iatom)**2 +&
                       vel(2, iatom)**2 + vel(3, iatom)**2))
      current_b(:) = current_b(:) + 2./3.*e2*zv(ityp(iatom))**2*vel(:, iatom)*init_data%I_primo
   end do

   do iatom = 1, nat
      do jatom = 1, nat
         if (iatom > jatom) then
            u(1:3) = (tau(:, iatom) - tau(:, jatom))*alat
            call pbc(u(1:3), u_pbc(1:3), alat, at, bg)
            call I_due_value(init_data, val, u_pbc, 1, tpiba, at, alat, gstart, g, gg, npw)
            !current_c(:) = current_c(:) + 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))* &
            current_c(:) = current_c(:) + 1.d0*e2*zv(ityp(iatom))*zv(ityp(jatom))* &
                             (vel(:, iatom) + vel(:, jatom))*val

            do a = 1, 3
               do b = 1, 3
                  if (a > b) then
                     call I_uno_value(init_data, val, u_pbc, a, b, 1, tpiba, at, alat, gstart, g, gg, npw)
                     current_e(a) = current_e(a) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(b, jatom) + vel(b, iatom))*val
                     current_e(b) = current_e(b) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(a, jatom) + vel(a, iatom))*val
                  end if
                  if (a == b) then
                     call I_uno_value(init_data, val, u_pbc, a, b, 1, tpiba, at, alat, gstart, g, gg, npw)
                     current_d(a) = current_d(a) - 1./2.*e2*zv(ityp(iatom))*zv(ityp(jatom))*&
                                      &(vel(b, jatom) + vel(b, iatom))*val
                  end if
               end do
            end do
         end if
      end do
   end do
   current_a = current_a * alat**3
   current_b = current_b * alat
   current_c = current_c * alat
   current_d = current_d * alat
   current_e = current_e * alat
   current   = current_a + current_c + current_d + &
                 current_e
   if (add_current_b) &
      current = current + current_b


   call stop_clock('calcolo_i')
   call print_clock('calcolo_i')
   if (ionode) print *, 'IONIC CURRENT CALCULATED'

   end subroutine

end module
