!
! Copyright (C) 2013-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
program generate_rVV10_kernel_table_
  !--------------------------------------------------------------------------------------------
  use mp_global,            ONLY : mp_startup, mp_global_end
  use mp_world,             ONLY : world_comm
  
  implicit none  

  call mp_startup ()
  CALL generate_rVV10_kernel_table ( world_comm )
  call mp_global_end( )

END program GENERATE_RVV10_KERNEL_TABLE_

!----------------------------------------------------------------------------
SUBROUTINE generate_rVV10_kernel_table ( my_comm )
  !--------------------------------------------------------------------------------------------
  !! These are the user set-able parameters.  
  use mp,                   ONLY : mp_get, mp_barrier, mp_size, mp_rank
  use kinds,                ONLY : dp
  use io_global,            ONLY : ionode
  use constants,            ONLY : pi  
  !
  implicit none
  !
  integer, intent (in) :: my_comm
  !
  integer, parameter :: Nr_points = 1024         !! The number of radial points (also the number of k points) used
  real(dp), parameter :: r_max =100.0D0         !! The value of the maximum radius to use for the real-space kernel functions 
  
  !!-------------------------------------------------------------------------------------------------

  CHARACTER(LEN=30) :: double_format = "(1p4e23.14)"

  !!        integer, parameter :: Nqs = 20
  !!     
  !!        real(dp), dimension(Nqs):: q_mesh = (/ 1.0D-4, 1.0D-3, 2.180901130847297D-3, 3.730376109554437D-3, &
  !!             5.763461392415386D-3, 8.431097736459265D-3, 1.193133637576012D-2, 1.652404278396654D-2, &
  !!             2.255018966296721D-2, 3.045717151192698D-2, 4.083202018639992D-2, 5.444498744423845D-2, &
  !!             7.230673014309263D-2, 9.574334364529304D-2, 0.126494814076381D0, 0.166844198751147D0, &
  !!             0.219787125408285D0, 0.289254194252484D0, 0.380402794424442D0, 0.5D0 /)
  !!   

  integer, parameter :: Nqs = 20
     
  real(dp), dimension(Nqs):: q_mesh = (/ 1.0D-4, 3.0D-4, 5.893850845618885D-4, 1.008103720396345D-3, &
             1.613958359589310D-3, 2.490584839564653D-3, 3.758997979748929D-3, 5.594297198907115D-3, &
             8.249838297569416D-3, 1.209220822453922D-2, 1.765183095571029D-2, 2.569619042667097D-2, &
             3.733577865542191D-2, 5.417739477463518D-2, 7.854595729872216D-2, 0.113805449932145D0, &
             0.164823306218807D0, 0.238642339497217D0, 0.345452975434964D0, 0.5D0 /)
   
   
  !!        integer, parameter :: Nqs = 20
  !!     
  !!        real(dp), dimension(Nqs):: q_mesh = (/ 1.0D-4, 5.0D-3, 1.069893648719707D-2, 1.732707466783098D-2, &
  !!             2.503591926370824D-2, 3.400167757831576D-2, 4.442928720984789D-2, 5.655710047584504D-2, &
  !!             7.066233262513155D-2, 8.706739837124464D-2, 0.106147281586631D0, 0.128338106612676D0, &
  !!             0.154147107106948D0, 0.184164220293662D0, 0.219075571636537D0, 0.259679158164142D0, &
  !!             0.306903088934189D0, 0.361826799573841D0, 0.425705725813956D0, 0.5D0 /)
  !!   
  !!   
  !!   
  !!        integer, parameter :: Nqs = 20
  !!     
  !!        real(dp), dimension(Nqs):: q_mesh = (/ 1.0D-4, 1.0D-3, 2.236206697581317D-3, 3.934214474408992D-3, &
  !!             6.266535125808476D-3, 9.470124470438768D-3, 1.387045625280778D-2, 1.991458916496835D-2, &
  !!             2.821658648395121D-2, 3.961990280509427D-2, 5.528307615046694D-2, 7.679743147816347D-2, &
  !!             0.106348753867322D0, 0.146939356822725D0, 0.202693107080873D0, 0.279274395396871D0, &
  !!             0.384463619314253D0, 0.528947645003195D0, 0.727405556392285D0, 1.0D0 /)
  !!   
  !!   
  !!        integer, parameter :: Nqs = 20
  !!     
  !!        real(dp), dimension(Nqs):: q_mesh = (/ 1.0D-4, 5.0D-4, 1.078773837542156D-3, 1.916221725100339D-3, &
  !!             3.127954044159392D-3, 4.881251455098552D-3, 7.418158132303510D-3, 1.108889616493409D-2, &
  !!             1.640021400932794D-2, 2.408534353734285D-2, 3.520522330968714D-2, 5.129496203180489D-2, &
  !!             7.457576159493165D-2, 0.108261555855433D0, 0.157002696892448D0, 0.227527940002890D0, &
  !!             0.329573353999449D0, 0.477226393655365D0, 0.690870684621412D0, 1.0D0 /)
  !!   
  
  !!                                DO NOT CHANGE ANYTHING BELOW THIS LINE
  !! #########################################################################################################


  integer  :: q1_i, q2_i, r_i, count                       !! Indexing variables
  real(dp) :: dr, d1, d2                                   !! Intermediate values
  
  real(dp) :: gamma = 4.0D0*pi/9.0D0                       !! Multiplicative factor for exponent in the functions called
  !                                                        !! "h" in DION

  real(dp), parameter :: small = 1.0D-15                   !! Number at which to employ special algorithms to avoid numerical
  !                                                        !! problems.  This is probably not needed but I like to be careful.

  !! The following sets up a parallel run.
  !! ------------------------------------------------------------------------------------------------------------------------------------------

  integer :: my_start_q, my_end_q, Ntotal                  !! starting and ending q value for each processor, also the total number of 
  !                                                        !! calculations to do (  (Nqs^2 + Nqs)/2  )

  real(dp), allocatable :: phi(:,:), d2phi_dk2(:,:)        !! Arrays to store the kernel functions and their second derivatives.  They are
  !                                                        !! stored as phi(radial_point, index)

  integer, allocatable :: indices(:,:), proc_indices(:,:)  !! indices holds the values of q1 and q2 as partitioned out to the processors.  It is an 
  !                                                        !! Ntotal x 2 array stored as indices(index of point number, q1:q2).  
  !                                                        !! Proc_indices holds the section of the indices array that is assigned to each processor.  
  !                                                        !! This is a Nprocs x 2 array, stored as proc_indices(processor number, starting_index:ending_index)

  integer :: Nper, Nextra, start_q, end_q                  !! Baseline number of jobs per processor, number of processors that get an extra job in case the 
  !                                                        !! number of jobs doesn't split evenly over the number of processors, starting index into the
  !                                                        !! indices array, ending index into the indices array.

  integer :: index, proc_i, kernel_file, my_Nqs
  integer :: nproc, Nprocs, mpime                          !! Variables holding information about the parallel run. Available and total number of processors, the rank of
  !                                                        !! this particular processor

  ! The total number of phi_alpha_beta functions that have to be calculated
  Ntotal = (Nqs**2 + Nqs)/2

  allocate( indices(Ntotal, 2) )

  count = 1

  ! This part fills in the indices array.  It just loops through the q1 and q2 values and stores them.  Sections
  ! of this array will be assigned to each of the processors later.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  do q1_i = 1, Nqs
     do q2_i = 1, q1_i

        indices(count, 1) = q1_i
        indices(count, 2) = q2_i

        count = count + 1

     end do
  end do

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! Figure out the baseline number of functions to be calculated by each processor and how many processors get 1 extra job.
  nproc = mp_size (my_comm)
  mpime = mp_rank (my_comm)
  Nper = Ntotal/nproc
  Nextra = mod(Ntotal, nproc)


  allocate(proc_indices(nproc,2) )

  start_q = 0
  end_q = 0

  ! Loop over all the processors and figure out which section of the indices array each processor should do.  All processors
  ! figure this out for every processor so there is no need to communicate results.
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  do proc_i = 1, nproc

     start_q = end_q + 1
     end_q = start_q + (Nper - 1)

     if (proc_i <= Nextra) end_q = end_q + 1

     if (proc_i == (mpime+1)) then

        my_start_q = start_q
        my_end_q = end_q

     end if

     proc_indices(proc_i, 1) = start_q
     proc_indices(proc_i, 2) = end_q

  end do

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Store how many jobs are assigned to me
  my_Nqs = my_end_q - my_start_q + 1
  
  
  !! ------------------------------------------------------------------------------------------------------------------------------------------


  allocate( phi(0:Nr_points, my_Nqs), d2phi_dk2(0:Nr_points, my_Nqs) )
  
  phi = 0.0D0
  d2phi_dk2 = 0.0D0
  
  dr = (r_max)/(Nr_points)



  !! Now, we loop over all the pairs q1,q2 that are assigned to us and perform our calculations
  !! -----------------------------------------------------------------------------------------------------

  do index = 1, my_Nqs

     ! First, get the value of phi(q1*r, q2*r) for each r and the particular values of q1 and q2 we are using
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     do r_i = 1, Nr_points

        d1 = q_mesh(indices(index+my_start_q-1, 1)) * (dr * r_i)**2      !! Different definition of d1 and d2 for vv10 !!!!  
        d2 = q_mesh(indices(index+my_start_q-1, 2)) * (dr * r_i)**2      !! Different definition of d1 and d2 for vv10 !!!! 

        phi(r_i, index) = phi_value(d1, d2)
     
     end do

     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Now, perform a radial FFT to turn our phi_alpha_beta(r) into phi_alpha_beta(k) needed for SOLER 
     ! equation 11
     call radial_fft( phi(:,index) )

     ! Determine the spline interpolation coefficients for the Fourier transformed kernel function
     call set_up_splines( phi(:, index), d2phi_dk2(:, index) )                                      
     
  end do

 
  !! -----------------------------------------------------------------------------------------------------
  !! Finally, we write out the results, after letting everybody catch up
  !! -----------------------------------------------------------------------------------------------------

  call mp_barrier(my_comm)

  call write_kernel_table_file(phi, d2phi_dk2)

  !! -----------------------------------------------------------------------------------------------------

  deallocate( phi, d2phi_dk2, indices, proc_indices )


CONTAINS



  !! ###########################################################################################################
  !!                                           |                  |
  !!                                           |  SET UP SPLINES  |
  !!                                           |__________________|

  !! This subroutine accepts a function (phi) and finds at each point the second derivative (D2) for use with
  !! spline interpolation.  This function assumes we are using the expansion described in SOLER 3 and 4.  That
  !! is, the derivatives are those needed to interpolate Kronecker delta functions at each of the q values
  !! Other than some special modification to speed up the algorithm in our particular case, this algorithm is
  !! taken directly from NUMERICAL_RECIPES pages 96-97.

  subroutine set_up_splines(phi, D2)

    real(dp), intent(in) :: phi(0:Nr_points)          !! The k-space kernel function for a particular q1 and q2

    real(dp), intent(inout) :: D2(0:Nr_points)        !! The second derivatives to be used in the interpolation
    !                                                 !! expansion (SOLER equation 3)

    real(dp), save :: dk = 2.0D0*pi/r_max             !! Spacing of k points

    real(dp), allocatable :: temp_array(:)            !! Temporary storage
    real(dp) :: temp_1, temp_2                        !! 

    allocate( temp_array(0:Nr_points) )

    D2 = 0

    temp_array = 0

    do r_i = 1, Nr_points - 1

       temp_1 = dble(r_i - (r_i - 1))/dble( (r_i + 1) - (r_i - 1) )
       temp_2 = temp_1 * D2(r_i-1) + 2.0D0

       D2(r_i) = (temp_1 - 1.0D0)/temp_2

       temp_array(r_i) = ( phi(r_i+1) - phi(r_i))/dble( dk*((r_i+1) - r_i) ) - &
            ( phi(r_i) - phi(r_i-1))/dble( dk*(r_i - (r_i-1)) )
       temp_array(r_i) = (6.0D0*temp_array(r_i)/dble( dk*((r_i+1) - (r_i-1)) )-temp_1*temp_array(r_i-1))/temp_2

    end do

    D2(Nr_points) = 0.0D0

    do  r_i = Nr_points-1, 0, -1

       D2(r_i) = D2(r_i)*D2(r_i+1) + temp_array(r_i)

    end do

    deallocate( temp_array )

  end subroutine set_up_splines


  !! ###########################################################################################################




  !! ###########################################################################################################
  !!                                          |             |
  !!                                          |  PHI_VALUE  |
  !!                                          |_____________|

  !! vv10 kernel phi


  real(dp) function phi_value(d1, d2)

    real(dp), intent(in) :: d1, d2                !! The point at which to evaluate the kernel.  Note that
    !                                             !! d1 = q1*r^2 and d2 = q2*r^2


     phi_value = - 24.0D0 / ( ( d1 + 1.0 ) * ( d2 + 1.0 ) * ( d1 + d2 + 2.0 ) ) 


    return

  end function phi_value


  !! ###########################################################################################################
  !!                                   |                           |
  !!                                   |  WRITE_KERNEL_TABLE_FILE  |
  !!                                   |___________________________|

  !! Subroutine to write out the vdW_kernel_table file.  All processors pass their data to processor 0 which
  !! is the one that actually does the writing.  This is the only communication in the entire program.

  
  subroutine write_kernel_table_file(phi, d2phi_dk2)
    
    real(dp), target :: phi(:,:), d2phi_dk2(:,:)      !! Each processor passes in its array of kernel values and second
    !                                                 !! derivative values for the q-pairs it calculated.  They are stored
    !                                                 !! as phi(index of function, function_values)

    integer :: proc_Nqs                               !! Number of calculated functions for a particular processor

    real(dp), pointer :: data(:,:)                    !! Pointer to point to the needed section of the phi and d2phi_dk2
    !                                                 !! arrays.  This is needed because some processors may have calculated
    !                                                 !! 1 extra function if the number of processors is not an even divisor
    !                                                 !! of (Nqs^2+Nqs)/2.  Processor 0 is guaranteed to be one of the ones
    !                                                 !! with an extra calculation (if there are any), so it can collect the
    !                                                 !! arrays from other processors and put it in its array.  Data then 
    !                                                 !! points to either the entire array (if the other processor also had
    !                                                 !! an extra calculation), or just the first proc_Nqs entries (which is
    !                                                 !! guaranteed to be at most 1 less than the proc_Nqs for processor 0.


    if (ionode) then
       
       !! Open the file for writing.  The file is written in binary to save space.
       !open(UNIT=21, FILE='vdW_kernel_table', status='replace', form='unformatted', action='write')
       open(UNIT=21, FILE='rVV10_kernel_table', status='replace', form='formatted', action='write') 
       !! Write the relevant header information that will be read in by the kernel_table module
       !! ---------------------------------------------------------------------------------------
       !write(*) "Writing headers..."

       write(21, '(2i5,f13.8)') Nqs, Nr_points
       write(21, double_format) r_max
       write(21, double_format) q_mesh
       !! ---------------------------------------------------------------------------------------



       !! Processor 0 writes its kernel functions first.  The subroutine "write_data" is defined
       !! below.
       !! ---------------------------------------------------------------------------------------
       !write(*) "Writing phi proc ", 0
       data => phi(:,:)
       call write_data(21, data)
 
       !! ---------------------------------------------------------------------------------------
       
    end if
    
    !! Now, loop over all other processors (if any) and collect their kernel functions in the phi
    !! array of processor 0, which is big enough to hold any of them.  Figure out how many functions
    !! should have been passed and make data point to just the right amount of the phi array.  Then
    !! write the data.
    !! -------------------------------------------------------------------------------------------

    do proc_i = 1, nproc-1
       
       call mp_get(phi, phi, mpime, 0, proc_i, 0, my_comm)
       
       if (ionode) then
          
          proc_Nqs = proc_indices(proc_i+1, 2) - proc_indices(proc_i+1,1) + 1
          
          !write(*) "Writing phi proc ", proc_i
          data => phi(:,1:proc_Nqs)
          call write_data(21, data)
          
       end if
       
    end do
    
    !! -------------------------------------------------------------------------------------------


    !! Here, we basically repeat the process exactly but for the second derivatives d2phi_dk2 
    !! instead of the kernel itself
    !! -------------------------------------------------------------------------------------------
    
    if (ionode) then
       
       !write(*) "Writing d2phi_dk2 proc ", 0
       data => d2phi_dk2(:,:)
       call write_data(21, data)

    end if
    
    
    do proc_i = 1, nproc-1
       
       call mp_get(d2phi_dk2, d2phi_dk2, mpime, 0, proc_i, 0, my_comm)
       
       if (mpime == 0) then
          
          proc_Nqs = proc_indices(proc_i+1,2) - proc_indices(proc_i+1,1) + 1
          
          !write(*) "Writing d2phi_dk2 proc ", proc_i 
          data => d2phi_dk2(:, 1:proc_Nqs)
          call write_data(21, data)
          
       end if
       
    end do
    
    !! -------------------------------------------------------------------------------------------

    if (ionode) then
       
       close(21)
       
    end if
    
  end subroutine write_kernel_table_file
  

  !! ###########################################################################################################
  !!                                      |              |
  !!                                      |  WRITE_DATA  |
  !!                                      !______________|

  !! Write matrix data held in the point "array" to the file with unit number "file".  Data is written
  !! in binary format.

  subroutine write_data(file, array)
    
    real(dp), pointer:: array(:,:)     !! Input pointer to the matrix data to be written

    integer, intent(in) :: file                    !! Unit number of file to write to

    integer :: index, ios                              !! Indexing variable
    

    do index = 1, size(array,2)
       
       ! write(file) array(:,index)
       write (file, double_format, err=100, iostat=ios) array(:,index)    
    
    end do

  100 call errore ('generate_vv10_sgd_kernel_table', 'Writing table file', abs (ios) )
  
  end subroutine write_data
  

  !! ###########################################################################################################
  !!                                            |              |
  !!                                            |  RADIAL_FFT  |
  !!                                            |______________|   

  !! This subroutine performs a radial Fourier transform on the real-space kernel functions.  Basically, this is
  !! just int( 4*pi*r^2*phi*sin(k*r)/(k*r))dr integrated from 0 to r_max.  That is, it is the kernel function phi
  !! integrated with the 0^th spherical Bessel function radially, with a 4*pi assumed from angular integration
  !! since we have spherical symmetry.  The spherical symmetry comes in because the kernel function depends only
  !! on the magnitude of the vector between two points.  The integration is done using the trapezoid rule.


  subroutine radial_fft(phi)

    real(dp), intent(inout) :: phi(0:Nr_points)      !! On input holds the real-space function phi_q1_q2(r)
    !                                                !! On output hold the reciprocal-space function phi_q1_q2(k)

    real(dp) :: phi_k(0:Nr_points)                   !! Temporary storage for phi_q1_q2(k)

    real(dp) :: dr = r_max/Nr_points                 !! Spacing between real-space sample points

    real(dp) :: dk = 2.0D0*pi/r_max                  !! Spacing between reciprocal space sample points

    integer :: k_i, r_i                              !! Indexing variables

    real(dp) :: r, k                                 !! The real and reciprocal space points

    phi_k = 0.0D0

    !! Handle the k=0 point separately
    !! -------------------------------------------------------------------------------------------------
    
    do r_i = 1, Nr_points

       r = r_i * dr

       phi_k(0) = phi_k(0) + phi(r_i)*r**2

    end do
    
    !! Subtract half of the last value of because of the trapezoid rule
    phi_k(0) = phi_k(0) - 0.5D0 * (Nr_points*dr)**2 * phi(Nr_points)

    !! -------------------------------------------------------------------------------------------------


    !! Integration for the rest of the k-points
    !! -------------------------------------------------------------------------------------------------
    
    do k_i = 1, Nr_points

       k = k_i * dk

       do r_i = 1, Nr_points

          r = r_i * dr

          phi_k(k_i) = phi_k(k_i) + phi(r_i) * r * sin(k*r) / k

       end do

       phi_k(k_i) = phi_k(k_i) - 0.5D0 * phi(Nr_points) * r *sin(k*r) / k

    end do
    
    !! Add in the 4*pi and the dr factor for the integration
    phi = 4.0D0 * pi * phi_k * dr

    !! -------------------------------------------------------------------------------------------------

  end subroutine radial_fft


  !! ###########################################################################################################

end subroutine generate_rVV10_kernel_table
