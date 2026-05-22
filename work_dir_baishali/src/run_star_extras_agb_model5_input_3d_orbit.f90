! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

         use star_lib
         use star_def
         use const_def
         use math_lib
         
         implicit none

         !arrays
         !real(dp), allocatable :: v(:), Ra(:), mdot_hl(:), fd_hl(:), edot_hl(:), Eorb(:), dEorb(:)                
         !real(dp), allocatable :: eps_rho(:), fd_mr15_ratio(:), mdot_mr15_ratio(:), mdot_mr15(:), fd_mr15(:)
         !real(dp), allocatable :: mdot_edd(:), mdot_hyper(:), mdot_arr(:), fd_arr(:), edot(:), f(:)

         ! removed loop from all the variables
         real(dp) :: v, Ra, mdot_hl, fd_hl, edot_hl, Eorb, dEorb
         real(dp) :: eps_rho, fd_mr15_ratio, mdot_mr15_ratio, mdot_mr15, fd_mr15
         real(dp) :: mdot_edd, mdot_hyper, mdot_arr, fd_arr, edot
         real(dp), save :: f = 0.0d0, f_old = 0.0d0

         ! one loop for allocating the kernel 
         real(dp), allocatable :: kernel(:)

         !extra smoothing variables
         real(dp), allocatable :: R_smooth(:), m_smooth(:), rho_boundary_smooth(:), scale_height_smooth(:), v_smooth(:)
         real(dp), allocatable :: Ra_smooth(:), eps_rho_smooth(:), mr15_ratio_smooth(:),mdot_hl_smooth(:), fd_hl_smooth(:)
         real(dp), allocatable :: fd_mr15_smooth(:), fd_array_smooth(:), dEorb_smooth(:),f_smooth(:), edot_smooth(:)

         !extra reading force variables
         real(dp), allocatable ,save:: t_position(:), r_position(:), t_power(:), f_power(:), t_velocity(:), f_velocity(:)

         !xtra variables - values for current step
         integer, parameter :: a_curr = 1, M_ns_curr = 2, M_acc_curr = 3, omega_env = 4                    !values for the current timestep
         integer, parameter :: omega_curr = 5, Q_curr = 6, Qmax_curr = 7, Qtb_curr = 8       !values for the current timestep
         integer, parameter :: e_curr = 9, aa_curr = 10, bb_curr = 11, mom_inert_curr = 12   !values for the current timestep
         integer, parameter :: strain_curr = 13                                              !values for the current timestep
         integer, parameter :: a1_curr = 14, b1_curr = 15, c1_curr = 16, d1_curr = 17
         integer, parameter :: a2_curr = 18, b2_curr = 19, c2_curr = 20, d2_curr = 21

         !values for next step
         real(dp) :: a_next, M_ns_next, M_acc_next, mdot_next, omega_next, Q_next
         real(dp) :: Qmax_next, Qtb_next, e_next, aa_next, bb_next, strain_next, mom_inert_next
         real(dp) :: a1_next, b1_next, c1_next, d1_next, a2_next, b2_next, c2_next, d2_next

         !input from inlists
         real(dp) :: M_ns_initial, omega_initial, e_initial, D, R0, r1, r2
         real(dp) :: op_const, eta, efactor, M_crust, n_poly, beta_sec                      
         real(dp) :: omega_env_factor, prescription  

         !other variables
         real(dp) :: temp, fd, decay_coeff, Req, Rbar, beta, ebind, eorb_change, mdot       
         real(dp) :: v_rel, Ra_br, mdot_br, fd_br, edot_br, mdot_hl_br
         real(dp) :: mdot_mr15_ratio_br, mdot_mr15_br, fd_mr15_ratio_br, fd_mr15_br, fd_hl_br
         integer :: azone

         ! variables for interpolation 
         integer :: n_profile_points
         real(dp), allocatable :: prof_mass_target(:),prof_mass_target_values(:),prof_mass_y2(:), prof_rho(:), prof_y2(:)
         logical :: relax_density = .true.
         logical :: relax_mass = .true.

         contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! FUNCTIONS NOT IN USE !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !DONE: To find the closest lower value to a target in an array
      subroutine find_closest_value_index(arr, tar, n, closest_index)
         implicit none
         real(dp), intent(in) :: arr(:)         ! Array to search
         real(dp), intent(in) :: tar            ! Target value
         integer, intent(in) :: n               ! Size of the array
         integer, intent(out) :: closest_index  ! Index of the closest value
   
         integer :: l
         real(dp) :: min_diff, diff, temp
         
         ! print *, 'find closest value index called'
         ! Initialize closest_index and min_diff
         closest_index = 1
         min_diff = abs(arr(1) - tar)
   
         ! Loop through the array to find the closest value
         do l = 2, n
            !  temp = arr(l) - tar
               diff = abs(arr(l) - tar)
               if (diff < min_diff) then
                  min_diff = diff
                  closest_index = l
               end if
         end do
      end subroutine find_closest_value_index

      !DONE: To evaluate the function F for the differential equation da/dt = F(a)
      subroutine linspace(n, from, to, array)
         real(dp), intent(in) :: from, to
         integer, intent(in) :: n
         real(dp), allocatable, intent(out) :: array(:)
         real(dp) :: range
         integer :: i
         range = to - from
   
         if (n == 0) return
         allocate(array(n))
         if (n == 1) then
               array(1) = from
               return
         end if
   
         do i=1, n
               array(i) = from + range * (i - 1) / (n - 1)
         end do
      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! OPERATIONAL FUNCTIONS !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !forward euler
      subroutine forward_euler(x_n, dxdt_n, dt_n, x_m)
         implicit none
         real(dp), intent(in) :: x_n, dxdt_n, dt_n
         real(dp), intent(out) :: x_m

         x_m = x_n + dxdt_n*dt_n
      
      end subroutine forward_euler

      !backward euler
      subroutine backward_euler(x_n, dxdt_m, dt_n, x_m)
         implicit none
         real(dp), intent(in) :: x_n, dxdt_m, dt_n
         real(dp), intent(out) :: x_m

         x_m = x_n + dxdt_m*dt_n
      
      end subroutine backward_euler

      !trapezoidal rule
      subroutine trapezoidal_rule(x_n, dxdt_n, dxdt_m, dt_n, x_m)
         implicit none
         real(dp), intent(in) :: x_n, dxdt_n, dxdt_m, dt_n
         real(dp), intent(out) :: x_m

         x_m = x_n + (dxdt_n + dxdt_m)*dt_n/2.0_dp
      
      end subroutine trapezoidal_rule

      !find zone
      subroutine find_zone(array, target, index)
         implicit none
         real(dp), intent(in) :: array(:)
         real(dp), intent(in) :: target
         integer, intent(out) :: index
         integer :: i
         real(dp) :: diff

         index = 1
         diff = abs(array(1) - target)
         do i = 2, size(array)
            if (diff > abs(array(i) - target)) then
               diff = abs(array(i) - target)
               index = i
            end if
         end do
      end subroutine find_zone

      ! subroutine find_azone(array, target, index)
      !    implicit none
      !    real(dp), intent(in) :: array(:)
      !    real(dp), intent(in) :: target
      !    integer, intent(out) :: index
      !    integer :: i
      !    real(dp) :: diff

      !    index = 1
      !    do i = 1, size(array)
      !       diff = array(i) - target
      !       if (diff .le. 0) then
      !          index = i
      !          return
      !       end if
      !    end do
      ! end subroutine find_azone

      !DONE: To append elements to an allocatable array
      subroutine AddToList(list, element)
         integer :: i, isize
         integer, intent(in) :: element
         integer, dimension(:), allocatable, intent(inout) :: list
         integer, dimension(:), allocatable :: clist

         if(allocated(list)) then
             isize = size(list)
             allocate(clist(isize+1))
             do i=1,isize          
             clist(i) = list(i)
             end do
             clist(isize+1) = element

             deallocate(list)
             call move_alloc(clist, list)

         else
             allocate(list(1))
             list(1) = element
         end if
      end subroutine AddToList
      
      !DONE: To find array indices where array values lie in a given range
      subroutine find_inbetween_value_index(arr, v1, v2, n, ind_btw)

         real(dp), intent(in) :: arr(:)
         real(dp), intent(in) :: v1, v2
         integer, intent(in) :: n
         integer, allocatable, intent(out) :: ind_btw(:)
         
         integer :: l

         ! print *, 'find inbetween value index called'
         do l = 1, n
            if (arr(l) >= v1 .and. arr(l) <= v2) then
               call AddToList(ind_btw, l)
            end if
         end do
         ! print *, '*16*'
      end subroutine find_inbetween_value_index

      !DONE: fourth order runge kutta solver
      !currently, this considers the orbital_evolution_function_holgado as the function to be solved
      subroutine rk4(id, ierr, a_in, t, dt, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         real(dp), intent(in) :: a_in, t, dt
         real(dp), intent(out) :: a_out
         real(dp) :: k1, k2, k3, k4
         procedure(orbital_evolution_function_holgado), pointer :: G
         
         !getting pointer to star
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! print *, 'rk4 called'

         !setting the function pointer
         G => orbital_evolution_function_holgado
         
         !evaluating the k values
         call G(id, ierr, t, a_in, k1)
         call G(id, ierr, t + dt/2, a_in + dt*k1/2, k2)
         call G(id, ierr, t + dt/2, a_in + dt*k2/2, k3)
         call G(id, ierr, t + dt, a_in + dt*k3, k4)
         
         !evaluating the next value of a
         a_out = a_in + dt*(k1 + 2*k2 + 2*k3 + k4)/6
      end subroutine rk4

      subroutine rk4_small_step(id, ierr, a_in, t, dt, a_out,R_smooth, f_smooth, n_smooth, v_smooth, fd_array_smooth, dEorb_smooth, m_smooth, rho_boundary_smooth)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         real(dp), intent(in) :: a_in, t, dt
         real(dp), intent(out) :: a_out
         real(dp) :: k1, k2, k3, k4
         real(dp) :: dt_small , t_curr
         real(dp) :: a_temp
         integer :: i
         integer, intent(in) :: n_smooth
         integer :: this_zone
         integer :: unit_log

         ! arrays for printing every small timestep
         real(dp), intent(in) :: R_smooth(:), f_smooth(:), v_smooth(:)
         real(dp), intent(in) :: fd_array_smooth(:), dEorb_smooth(:)
         real(dp), intent(in) :: m_smooth(:), rho_boundary_smooth(:)

         ! file variables
         logical :: file_exists
         integer :: filesize
         
         !getting pointer to star
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! print *, 'rk4 called'

         dt_small = dt/10
         a_temp = a_in
         t_curr = t
         

         ! check if file exists and has content
         inquire(file='orbit_substeps.data', exist=file_exists, size=filesize)

         unit_log = 255
         ! open the file in append mode
         open(unit=unit_log, file='orbit_substeps.data', status='unknown', position='append', action='write')

         ! only write header if file is new or empty
         if (.not. file_exists .or. filesize == 0) then
            write(unit_log,'(a)') 'model    step    time        a_temp       zone         R          f(k)        v         fd         dEorb        m         rho              problematic_value'
         end if
         
         !evaluating the k values  ! not using pointer to avoid extra arguments problem
         do i = 1,10
            call smooth_orbital_evolution_function_holgado(id, ierr, t_curr, a_temp, k1, R_smooth, f_smooth, n_smooth)
            call smooth_orbital_evolution_function_holgado(id, ierr, t_curr + dt_small/2, a_temp + dt_small*k1/2, k2, R_smooth, f_smooth, n_smooth)
            call smooth_orbital_evolution_function_holgado(id, ierr, t_curr + dt_small/2, a_temp + dt_small*k2/2, k3, R_smooth, f_smooth, n_smooth)
            call smooth_orbital_evolution_function_holgado(id, ierr, t_curr + dt_small, a_temp + dt_small*k3, k4, R_smooth, f_smooth, n_smooth)
            a_temp = a_temp + dt_small*(k1 + 2*k2 + 2*k3 + k4)/6
            t_curr = t_curr + dt_small

            ! Find the orbital zone corresponding to current separation
            call find_zone(R_smooth(1:n_smooth), a_temp, this_zone)

            
            ! Write quantities to file
            
            write(unit_log,'(I6,1X,I6,1X,ES12.4,1X,ES12.4,1X,I6,1X,8(1X,ES12.4))') &
               s%model_number, i, t_curr, a_temp, this_zone, &
               R_smooth(this_zone), f_smooth(this_zone), v_smooth(this_zone), &
               fd_array_smooth(this_zone), dEorb_smooth(this_zone), &
               m_smooth(this_zone), rho_boundary_smooth(this_zone), &
               m_smooth(this_zone)/R_smooth(this_zone) - 4.0_dp*pi*(R_smooth(this_zone)**2)*rho_boundary_smooth(this_zone)
         end do

         close(unit_log)

         a_out = a_temp

      end subroutine rk4_small_step


      !Cash-Karp method(rk5 with embedded rk4)
      subroutine integrate_orbit_adaptive2(id, ierr, a_in, t_start, dt_mesa, a_out)
         use const_def, only: dp
         implicit none
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         real(dp), intent(in) :: a_in, t_start, dt_mesa
         real(dp), intent(out) :: a_out
         type(star_info), pointer :: s

         ! RKCK / Cash-Karp coefficients
         real(dp), parameter :: c1=37.0d0/378.0d0, c3=250.0d0/621.0d0, c4=125.0d0/594.0d0, c6=512.0d0/1771.0d0
         real(dp), parameter :: d1= 2825.0d0/27648.0d0, d3=18575.0d0/48384.0d0, d4=13525.0d0/55296.0d0, d5=277.0d0/14336.0d0, d6=1.0d0/4.0d0

         ! Local variables
         real(dp) :: t, a, dt, dt_min, dt_max, rel_tol, safety, err_est, abs_tol, err_norm, tol
         real(dp) :: k1, k2, k3, k4, k5, k6, a_5th, a_4th
         integer :: nsteps
         real(dp), parameter :: abs_max_change = 0.05  ! max 5% change per substep
         real(dp), parameter :: min_increment = 1d-8 ! minimum absolute change
         procedure(orbital_evolution_function_holgado), pointer :: G

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         G => orbital_evolution_function_holgado

         ! ---- Adaptive RK settings ----
         rel_tol = 1d-10             ! relative error tolerance ! should i change it to 1d-5??
         nsteps = 10             ! minimum number of substeps per MESA dt
         dt_max = dt_mesa / nsteps
         dt_min = dt_mesa *1e-5_dp  ! prevent dt from being ridiculously small(lol)
         safety = 0.85           ! timestep safety factor
         t = t_start
         a = a_in
         dt = dt_max
         abs_tol = 1d-10
         tol = 1d-10

         ! ---- Integration loop ----
         do while (t < t_start + dt_mesa)

            if (t + dt > t_start + dt_mesa) dt = t_start + dt_mesa - t

            ! Compute RKCK k-values
            call G(id, ierr, t,a,k1)
            call G(id, ierr, t+dt/5, a+dt*k1/5,  k2)
            call G(id, ierr, t+3*dt/10, a+dt*(3*k1/40 + 9./40*k2), k3)
            call G(id, ierr, t+3*dt/5, a+dt*(3*k1/10 - 9*k2/10 + 6*k3/5), k4)
            call G(id, ierr, t+dt, a+dt*(-11*k1/54 + 5*k2/2 - 70*k3/27 + 35*k4/27), k5)
            call G(id, ierr, t+7*dt/8, a+dt*(1631*k1/55296 + 175*k2/512 + 575*k3/13824 + 44275*k4/110592 + 253*k5/4096), k6)

            ! 5th and 4th order solutions
            a_5th = a + dt*(37*k1/378 + 250*k3/621 + 125*k4/594 + 512*k6/1771)
            a_4th = a + dt*(2825*k1/27648 + 18575*k3/48384 + 13525*k4/55296 + 277*k5/14336 + k6/4)

            

            ! testing for absolute and relative error tolerances
            err_est = abs(a_5th - a_4th) / (abs_tol + rel_tol*a_5th)

            ! inside the integration loop, after computing k1..k6 and a_5th, a_4th:
            !if (a_5th < 26.44670755724950*Rsun .and. a_5th > 24.0*Rsun) then
               !print *, "DBG_ORBIT: model=", s%model_number, " t=", t, " dt=", dt, &
               !" a=", a, " a5th=", a_5th, " k1=", k1, " k2=", k2, &
               !" k3=", k3, " k4=", k4, " k5=", k5, " k6=", k6, &
               !" err=", err_est
            !end if


            !if (s%model_number >= 38) then
               !print *, "Model dt:", dt, err_est
            !end if


            ! Accept/reject step based on relative error 
            if (err_est <= 1d0) then
               ! accept the step
               t = t + dt
               a = a_5th


               if (err_est > 0d0) then  ! choose next timestep
                  dt = min(dt_max, dt * safety * (err_est)**(-0.2d0))
               else
                  dt = min(dt_max, dt * 10_dp)   !making timsstep 1 order bigger if error is zero
               end if
            
            else
               ! error too large, decrease timestep
               dt = max(dt_min, dt * safety * err_est**(-0.25d0))  !decreasing more than increasing

               if (dt <= dt_min) then
                  print *, "WARNING: Minimum timestep reached and error still above tolerance!"
                  print *, " t=", t, " a=", a, " err_est=", err_est, " dt=", dt
                  ! At smallest timestep: accept anyway to keep going
                  t = t + dt
                  a = a_5th       ! accept best guess even if inaccurate
                  dt = dt_min     ! stay at min dt or slowly increase later
                  dt = min(dt_max, dt * 2.0d0) ! grow for the next step
               end if
            end if

         end do

         a_out = a

      end subroutine integrate_orbit_adaptive2


      ! this routine assumes everything is increasing
      subroutine cubic_spline(x, y, n, x_new, y_new, n_new)
         implicit none
         integer, intent(in) :: n, n_new
         real(dp), intent(in) :: x(n), y(n), x_new(n_new)
         real(dp), intent(out) :: y_new(n_new)

         ! Local variables
         integer :: i, j, k
         real(dp), allocatable :: a(:), b(:), c(:), d(:), h(:), alpha(:), l(:), mu(:), z(:)

         ! Allocate arrays
         allocate(a(n), b(n-1), c(n), d(n-1), h(n-1), alpha(n-1), l(n), mu(n), z(n))

         ! Copy y into a
         a = y

         ! Step 1: Compute h and alpha
         do i = 1, n-1
            h(i) = x(i+1) - x(i)
         end do

         do i = 2, n-1
            alpha(i) = (3.0_dp/h(i))*(a(i+1)-a(i)) - (3.0_dp/h(i-1))*(a(i)-a(i-1))
         end do

         ! Step 2: Solve tridiagonal system for c
         l(1) = 1.0_dp
         mu(1) = 0.0_dp
         z(1) = 0.0_dp

         do i = 2, n-1
            l(i) = 2.0_dp*(x(i+1)-x(i-1)) - h(i-1)*mu(i-1)
            mu(i) = h(i)/l(i)
            z(i) = (alpha(i)-h(i-1)*z(i-1))/l(i)
         end do

         l(n) = 1.0_dp
         z(n) = 0.0_dp
         c(n) = 0.0_dp

         do j = n-1, 1, -1
            c(j) = z(j) - mu(j)*c(j+1)
            b(j) = (a(j+1)-a(j))/h(j) - h(j)*(c(j+1)+2.0_dp*c(j))/3.0_dp
            d(j) = (c(j+1)-c(j))/(3.0_dp*h(j))
         end do

         ! Step 3: Interpolate y_new at x_new
         do k = 1, n_new
            ! Find the interval
            j = 1
            do i = 1, n-1
               if (x_new(k) >= x(i) .and. x_new(k) <= x(i+1)) then
                  j = i
                  exit
               end if
            end do

            y_new(k) = a(j) + b(j)*(x_new(k)-x(j)) + c(j)*(x_new(k)-x(j))**2 + d(j)*(x_new(k)-x(j))**3
         end do

         ! Deallocate
         deallocate(a,b,c,d,h,alpha,l,mu,z)
      end subroutine cubic_spline




      ! subroutine is modified to cap the outer energy injection radius at the outer envelope of star
      subroutine find_azone_r1r2(array, a_in, r1_out, r2_out, index_out, rzones_out, nr_out)
         implicit none
         real(dp), intent(in) :: array(:)
         real(dp), intent(in) :: a_in
         real(dp), intent(out) :: r1_out, r2_out
         integer, intent(out) :: index_out, nr_out
         integer, allocatable, intent(out) :: rzones_out(:)

         !zone array - to find the zone closest to s%xtra(a_curr)
         call find_zone(array, a_in, index_out)

         !finding the radius range where energy is to be injected

         if (prescription == 1) then
            r1_out = max(array(index_out) - Ra, array(size(array)))      !cm
            r2_out = min(array(index_out) + Ra, array(1))      !cm
         else if (prescription == 2) then
            r1_out = array(index_out) - Ra_br       !cm
            r2_out = array(index_out) + Ra_br       !cm
         end if

         call find_inbetween_value_index(array, r1_out, r2_out, size(array), rzones_out)
         nr_out = size(rzones_out)
      end subroutine find_azone_r1r2

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! EVALUATION FUNCTIONS !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine other_cgrav(id, ierr)  ! cannot have any extra arguments in the call
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ! subroutine variables
         integer :: k , azone_temp

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         

         call find_zone(s%R(1:s%nz), s%xtra(a_curr), azone_temp)

         if (s%x_ctrl(20) == 0) then ! companion's gravitational effect disabled
            do k = 1, s%nz
               s%cgrav(k) = standard_cgrav
            end do
            
               
         else if (s%x_ctrl(20) == 1) then  ! companion's gravitational effect enabled
            do k = 1, s%nz
               if (k <= azone_temp) then   
                  s%cgrav(k) = standard_cgrav * (1.0 + s%xtra(M_ns_curr)/s%m(k))
               else
                  s%cgrav(k) = standard_cgrav  
               end if
            end do
         end if

         !call write_cgrav_profile(s)   

      end subroutine other_cgrav


      ! writes profiles for cgrav, to see if the values are correct
      subroutine write_cgrav_profile(s)
         use star_def, only: star_info
         use const_def, only: dp, standard_cgrav
         implicit none

         type (star_info), pointer :: s


         integer :: k, iounit
         logical :: dir_exists
         character(len=256) :: fname, dirname
         real(dp) :: factor

         if (s%model_number < 1 .or. s%nz <= 0) return

         ! directory for output (inside work/)
         dirname =  'cgrav_profiles'

         ! Check if directory exists; if not, create it via execute_command_line
         inquire(file=trim(dirname), exist=dir_exists)
         if (.not. dir_exists) call execute_command_line('mkdir -p ' // trim(dirname))
         
         ! create a unique file name for each model
         write(fname, '(A,I6.6,A)') 'cgrav_profiles/cgrav_profile_', s%model_number, '.dat'

         open(newunit=iounit, file=fname, status='replace', action='write', form='formatted')
         write(iounit,'(A)') 'k   radius(cm)   mass(g)   factor(1+M2/Mr)   cgrav(cgs)   standard_cgrav    M   g_analytic'

         do k = 1, s%nz
            factor = 1 + s%xtra(M_ns_curr)/s%m(k)
            write(iounit,'(I5,1X,7(1PE14.6,1X))') k, s%r(k), s%m(k), factor, s%cgrav(k), standard_cgrav, s%xtra(M_ns_curr), 1.0_dp + s%xtra(M_ns_curr)/s%m(k)
         end do

         close(iounit)

         
      end subroutine write_cgrav_profile



      ! both updates to a file and prints relevant variables to the output log
      subroutine holgado_prescription(id, ierr, M, a_in, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k
         integer :: my_unit
         real(dp), intent(in) :: M
         real(dp), intent(out) :: a_out
         integer :: azone_temp
         real(dp), intent(in) :: a_in
         integer :: filesize
         logical :: file_exists

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         ! Check if file exists and its size
         inquire(file='dEorb.data', exist=file_exists, size=filesize)

         ! open the file in append mode
         open(newunit=my_unit, file='dEorb.data', status='unknown', action='write', position='append')
         
         if (.not. file_exists .or. filesize == 0) then
            write(my_unit,'(a)') 'model   zone        R(cm)            m/R          4πR²ρ        dEorb           m               rho                 f(k)         fd_arr(k)       v(k)         Eorb          M'
         end if

         !printing dEorb and f(k) at the zone of the orbit 
         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         !mass is in g, radius is in cm, time is in s

         ! exact 2 body problem result for orbital velocity
         ! removing loop
         k=azone_temp
         v = SQRT(s%cgrav(k)*(M + s%m(k))/s%R(k))   !cm/s
         

         call hl_and_energy(id, ierr, M, s%xtra(a_curr))            !all cgs
         call mr15(id, ierr, s%xtra(a_curr))                        !all cgs
         call edd_and_hyper(id, ierr, M, s%xtra(a_curr))            !all cgs

         
         if (eta*mdot_mr15 >= mdot_hyper .or. eta*mdot_mr15 <= mdot_edd) then
            mdot_arr = eta*mdot_mr15       !g/s
         else
            mdot_arr = mdot_edd            !g/s
         end if
         fd_arr = eta*fd_mr15              !dyne (g cm s^-2)
         edot = fd_arr*v                !erg/s
         f = -fd_arr*v*1/dEorb     !cm s^-1(da/dt)
         
            
            
         write(my_unit,'(i6,1x,i6,1x,12(3x,es12.4))') &
         s%model_number, k, s%R(k), s%m(k)/s%R(k), &
         4.0_dp*pi*s%R(k)**2*s%rho_face(k), dEorb, s%m(k), s%rho_face(k), f, fd_arr, v, Eorb, M 
            

            
         ! writing in output log
         write(*,'(" Zone=",I6," R=",ES12.4," f(k)=",ES12.4," v=",ES12.4," fd=",ES12.4," dEorb=",ES12.4," m=",ES12.4," rho on boundary=",ES12.4," problematic_value=",ES12.4, "  M=",ES12.4)') &
         k, s%R(k), f, v, fd_arr, dEorb, s%m(k), s%rho_face(k), (s%m(k)/s%R(k)) - 4*pi*(s%R(k)**2)*s%rho_face(k), M


         close(my_unit)

         call rk4(id, ierr, s%xtra(a_curr), s%time, s%dt_next, a_out)  !cm
         !call integrate_orbit_adaptive2(id, ierr, s%xtra(a_curr), s%time, s%dt_next, a_out)

      end subroutine holgado_prescription

      subroutine interpolate_force(t, n, t_data, F_data, F_out)
         implicit none
         real(dp), intent(in)  :: t
         integer, intent(in) :: n
         real(dp), intent(in) :: t_data(n), F_data(n)
         real(dp), intent(out) :: F_out
         integer :: i

         if (t <= t_data(1)) then
            F_out = F_data(1)
            return
         end if

         if (t >= t_data(n)) then
            F_out = F_data(n)
            return
         end if

         do i = 1, n-1
            if (t >= t_data(i) .and. t <= t_data(i+1)) then
               ! piecewise linear interpolation
               F_out = F_data(i) + (F_data(i+1) - F_data(i)) * (t - t_data(i)) / (t_data(i+1) - t_data(i))
               return
            end if
         end do

      end subroutine interpolate_force



      subroutine holgado_simulation_force(id, ierr, M, a_out)
         ! star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(star_info), pointer :: s
         

         !subroutine variables
         integer :: k,j
         integer :: my_unit
         real(dp), intent(in) :: M
         real(dp), intent(out) :: a_out
         integer :: azone_temp
         real(dp) :: rho_bar, term1, term2, term3, term4, term5, term
         integer :: filesize
         logical :: file_exists
         logical, save :: r_loaded = .false.
         logical, save :: p_loaded = .false.
         logical, save :: v_loaded = .false.
         
         
         ! position and power reading variables
         integer :: ios, i
         real(dp) :: dummy
         integer :: n_position = 732727
         real(dp) :: r_current

         integer :: n_power = 1278
         real(dp) :: p_curr 

         integer :: n_velocity = 1278
         real(dp) :: v_curr

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         ! read position vs time data from file
         ! Read ONLY once
         if (.not. r_loaded) then

            open(newunit=my_unit, file='time_vs_separation.txt', &
               status='old', action='read', iostat=ios)
            if (ios /= 0) then
               write(*,*) 'Error: cannot open time_vs_separation.txt'
               ierr = 1
               return
            end if

            read(my_unit,*)  ! skip header

            do i = 1, n_position
               read(my_unit,*,iostat=ios) t_position(i), r_position(i)
               if (ios /= 0) then
                  write(*,*) 'Error reading position vs time table at line ', i
                  ierr = 1
                  close(my_unit)
                  return
               end if
            end do

            close(my_unit)
            r_loaded = .true.

            write(*,*) 'position vs time table loaded successfully.'
         end if

         ! read power data
         if (.not. p_loaded) then

            open(newunit=my_unit, file='Power_phi_vs_time.dat', &
               status='old', action='read', iostat=ios)
            if (ios /= 0) then
               write(*,*) 'Error: cannot open Power_phi_vs_time.dat'
               ierr = 1
               return
            end if

            read(my_unit,*)  ! skip header

            do i = 1, n_power
               read(my_unit,*,iostat=ios) dummy, t_power(i), f_power(i)
               if (ios /= 0) then
                  write(*,*) 'Error reading Power vs time table at line ', i
                  ierr = 1
                  close(my_unit)
                  return
               end if
            end do

            close(my_unit)
            p_loaded = .true.

            write(*,*) 'power vs time table loaded successfully.'
         end if

         ! read velocity data
         if (.not. v_loaded) then

            open(newunit=my_unit, file='vphi_vs_time.dat', &
               status='old', action='read', iostat=ios)
            if (ios /= 0) then
               write(*,*) 'Error: cannot open vphi_vs_time.dat'
               ierr = 1
               return
            end if

            read(my_unit,*)  ! skip header

            do i = 1, n_velocity
               read(my_unit,*,iostat=ios) dummy, t_velocity(i), f_velocity(i)
               if (ios /= 0) then
                  write(*,*) 'Error reading Power vs time table at line ', i
                  ierr = 1
                  close(my_unit)
                  return
               end if
            end do

            close(my_unit)
            v_loaded = .true.

            write(*,*) 'velocity vs time table loaded successfully.'
         end if


         ! writing dEorb.data
         ! Check if file exists and its size
         inquire(file='dEorb.data', exist=file_exists, size=filesize)

         ! open the file in append mode
         open(newunit=my_unit, file='dEorb.data', status='unknown', action='write', position='append')
         
         if (.not. file_exists .or. filesize == 0) then
            write(my_unit,'(a)') 'model   zone        R(cm)     edot              M            v_curr'
         end if

         !printing dEorb and f(k) at the zone of the orbit 
         

         !mass is in g, radius is in cm, time is in s

         ! exact 2 body problem result for orbital velocity
         ! removing loop
         v = 0.0d0   !cm/s
         

         !call hl_and_energy(id, ierr, M, s%xtra(a_curr))            !all cgs
         !call mr15(id, ierr, s%xtra(a_curr))                        !all cgs
         !call edd_and_hyper(id, ierr, M, s%xtra(a_curr))            !all cgs

         
         !if (eta*mdot_mr15 >= mdot_hyper .or. eta*mdot_mr15 <= mdot_edd) then
            !mdot_arr = eta*mdot_mr15       !g/s
         !else
            !mdot_arr = mdot_edd            !g/s
         !end if

         !no timeshift

         call interpolate_force(s%time , n_position, t_position, r_position, r_current)
         call interpolate_force(s%time , n_power, t_power, f_power, p_curr)
         call interpolate_force(s%time , n_velocity, t_velocity, f_velocity, v_curr)

         call find_zone(s%R(1:s%nz), r_current , azone_temp)

         k=azone_temp

         Ra = 2*s%cgrav(k)*M/(v_curr**2 + s%csound(k)**2 )        

         fd_arr =  0.0d0           !dyne (g cm s^-2) ! since the original phi data is negative, and we will be accounting for it in dEorb
         edot = p_curr      
         
         ! new terms in energy
         !rho_bar = 3*s%m(k)/(4*pi*s%R(k)**3)  ! mean density of enclosed mass  
         !term1 =  2*s%R(k)**2/(s%cgrav(k)*M*s%m(k))
         !term2 =  1/(1+ (3*s%rho_face(k)/rho_bar))

         !term3 = 0.0d0
         !do j = 1, k
            !term = s%rho(j)*s%u(j)*(s%R(j) - s%R(j+1))
            !term3 = term3 + term   
         !end do
           
         !term4 =  4*pi*s%cgrav(k)*M*( s%R(k)* s%rho(k) * s%u(k)/2  - term3)
         !term5 = edot + term4

         !f =  (f_curr_phi*v + term5)/((1/(term1*term2) - f_curr_r) )    !cm s^-1(da/dt)  ! already negative, no need for another negative sign  
            
         write(my_unit,'(i6,1x,i6,1x,18(3x,es12.4))') &
         s%model_number, k, s%R(k), edot , M, v_curr
            

            
         ! writing in output log
         write(*,'(" Zone=",I6," R=",ES12.4, "edot = " , ES12.4, "  M=",ES12.4)') &
         k, s%R(k), edot, M


         close(my_unit)

         a_out = r_current
         !call integrate_orbit_adaptive2(id, ierr, s%xtra(a_curr), s%time, s%dt_next, a_out)

      end subroutine holgado_simulation_force 




      ! remove loop if smoothing is required
      subroutine smooth_orbital_evolution_function_holgado(id, ierr, t, a_in, fa, R_smooth, f_smooth, n_smooth) 
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         real(dp), intent(in) :: a_in, t
         real(dp), intent(out) :: fa
         integer :: k, current_zone
         real(dp), intent(in) :: R_smooth(:), f_smooth(:)
         integer, intent(in) :: n_smooth
         
         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return 
         
         ! print *, 'orbital evolution function called'

         !finding the zone where the function is to be evaluated
         call find_zone(R_smooth(1:n_smooth), a_in, current_zone)  ! current_zone is on the smoother grid, used to not confuse

         !evaluating the function
         fa = f_smooth(current_zone)

      end subroutine smooth_orbital_evolution_function_holgado


      subroutine smooth_profile_holgado(id, ierr, M, a_in, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         !subroutine inputs and outputs
         real(dp), intent(in)  :: M       
         real(dp), intent(in)  :: a_in    
         real(dp), intent(out) :: a_out   

         ! Star pointer
         type(star_info), pointer :: s

         ! subroutine local variables
         integer :: k,i, p, q, rs_idx
         integer :: s_nz, n_smooth, this_zone
         real(dp), allocatable :: R_increasing(:), rho_increasing(:), m_increasing(:), h_increasing(:)
         real(dp):: da_dt, step
         real(dp) :: f1,f2,f3
         integer :: my_unit
         logical :: file_exists
         integer :: filesize
         integer :: n_refine

         ! variables of writing interpolated quantities into file
         integer :: prof_unit
         character(len=128) :: prof_name
         integer :: model_no
         character(len=20) :: temp_model, temp_nsmooth
         character(len=256) :: folder_name
         folder_name = 'smooth_profiles/'   ! my subfolder where i want the smoothed files to go to

         ! assigning the star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         model_no = s%model_number
         write(prof_name,'(A,I5.5,A)') trim(folder_name)//'smooth_profile', model_no, '.data'

         allocate(R_increasing(s%nz))
         allocate(rho_increasing(s%nz))
         allocate(m_increasing(s%nz))
         allocate(h_increasing(s%nz))

         s_nz = s%nz
         f1 = 1.91791946
         f2 = -1.52814698
         f3 = 0.75992092
         
         ! number of points in new grid
         n_refine = 10
         n_smooth = (s%nz-1)*n_refine + 1
         

         ! Check if file exists and its size
         inquire(file='dEorb.data', exist=file_exists, size=filesize)

         ! open the file in append mode
         open(newunit=my_unit, file='dEorb.data', status='unknown', action='write', position='append')
         
         if (.not. file_exists .or. filesize == 0) then
            write(my_unit,'(a)') 'model   zone        R(cm)            m/R          4πR²ρ        dEorb           m               rho                 f(k)         fd_arr(k)       v(k)'
         end if


         open(newunit=prof_unit, file=trim(prof_name), status='replace', action='write')

         write(temp_model,'(i0)') model_no
         write(temp_nsmooth,'(i0)') n_smooth

         ! Write header to file
         write(prof_unit,'(a)') '# smooth_profile_holgado output'
         write(prof_unit,'(a)') '# model_number = '//trim(temp_model)
         write(prof_unit,'(a)') '# n_smooth = '//trim(temp_nsmooth)
         write(prof_unit,'(a)') '# Columns: zone R(cm) rho(g/cm^3) m(g) scale_height(cm) v(cm/s) fd dyn dEorb f'
         

         ! this array is increasing, R_smooth(1)= R(s_nz) and R_smooth(s_nz) = R(1)
         ! 1 to s_nz, radius is increasing
         ! inserts 9 points inside each interval, divides each interval by 10
         rs_idx = 1
         do p = s%nz, 2, -1
            step = (s%R(p-1) - s%R(p)) / n_refine
            do q = 0, n_refine - 1
               R_smooth(rs_idx) = s%R(p) + q*step
               rs_idx = rs_idx + 1
            end do
         end do
         R_smooth(rs_idx) = s%R(1)  ! last point = surface
         

         do i = 1, s%nz
            R_increasing(i)   = s%R(s%nz - i + 1)        ! now R_increasing(1) = center
            rho_increasing(i) = s%rho_face(s%nz - i + 1)
            m_increasing(i)   = s%m(s%nz - i + 1)
            h_increasing(i)   = s%scale_height(s%nz - i + 1)
         end do


         ! reverse the array while interpolating so that radius always increases(inside to outside)

         !Interpolate density on the new grid as a function of radius(using rho_face instead of rho)
         ! cubic spline is breaking
         call cubic_spline(R_increasing(1:s%nz), rho_increasing(1:s%nz), s_nz, R_smooth, rho_boundary_smooth, n_smooth)

         ! interpolate mass also on the smoother grid(questionable??)  
         call cubic_spline(R_increasing(1:s%nz), m_increasing(1:s%nz), s_nz, R_smooth, m_smooth, n_smooth)

         !interpolate scale_height( scale height probably uses rho, not rho_face, option: calculate the values in the boundary and store them in an array and pass them here)
         call cubic_spline(R_increasing(1:s%nz), h_increasing(1:s%nz), s_nz, R_smooth, scale_height_smooth, n_smooth)

         ! now reverse all of the smooth arrays to match the original mesa ordering
         R_smooth = R_smooth(n_smooth:1:-1)
         rho_boundary_smooth = rho_boundary_smooth(n_smooth:1:-1)
         m_smooth = m_smooth(n_smooth:1:-1)
         scale_height_smooth = scale_height_smooth(n_smooth:1:-1)

         !compute velocity on this new grid using mass and radius
         do k= 1, n_smooth
            v_smooth(k) = SQRT(s%cgrav(k)*(M + m_smooth(k))/( R_smooth(k)))
         end do



         !comptuting necessary values for calculating drag force and dEorb  
         !every formula is same as original, only on a new grid

         
         do k = 1, n_smooth
            Ra_smooth(k) = (2*s%cgrav(k)*M)/(v_smooth(k)**2)
            eps_rho_smooth(k) = Ra_smooth(k)/(scale_height_smooth(k))
            mr15_ratio_smooth(k)= f1 + f2*eps_rho_smooth(k) + f3*(eps_rho_smooth(k)**2)
            mdot_hl_smooth(k) = pi*(Ra_smooth(k)**2)*rho_boundary_smooth(k)*v_smooth(k)
            fd_hl_smooth(k) = mdot_hl_smooth(k)*v_smooth(k)
            fd_mr15_smooth(k) = mr15_ratio_smooth(k)*fd_hl_smooth(k)
            fd_array_smooth(k) = eta*fd_mr15_smooth(k)
            edot_smooth(k) = fd_array_smooth(k)*v_smooth(k)
            dEorb_smooth(k) = ((s%cgrav(k)*M)/(2*R_smooth(k))) * ((m_smooth(k)/R_smooth(k)) - 4*pi*(R_smooth(k)**2)*rho_boundary_smooth(k))
            f_smooth(k) = - ((fd_array_smooth(k)*v_smooth(k))/(dEorb_smooth(k)))
         end do

         !printing dEorb and f(k) at the zone of the orbit 
         call find_zone(R_smooth(1:n_smooth), a_in, this_zone)

         do k=1,n_smooth
            if (k == this_zone) then
               write(my_unit,'(i6,1x,i6,1x,10(3x,es12.4))') &
                  s%model_number, this_zone, R_smooth(k), m_smooth(k)/R_smooth(k), &
                  4.0_dp*pi*R_smooth(k)**2*rho_boundary_smooth(k), dEorb_smooth(k), m_smooth(k), rho_boundary_smooth(k), f_smooth(k), fd_array_smooth(k), v_smooth(k)
            end if

            if (k == this_zone) then
               ! writing in output log
               write(*,'(" Zone=",I6," R=",ES12.4," f(k)=",ES12.4," v=",ES12.4," fd=",ES12.4," dEorb=",ES12.4," m=",ES12.4," rho=",ES12.4," problematic_value=",ES12.4)') &
               k, R_smooth(k), f_smooth(k), v_smooth(k), fd_array_smooth(k), dEorb_smooth(k), m_smooth(k), rho_boundary_smooth(k), m_smooth(k)/R_smooth(k) - 4*pi*(R_smooth(k)**2)*rho_boundary_smooth(k)
            end if
         end do

         ! loop over all zones to write the data in the smooth profile file
         do k = 1, n_smooth
            write(prof_unit,'(i8,1x,9(es20.10))') k, &
            R_smooth(k), &
            rho_boundary_smooth(k), &
            m_smooth(k), &
            scale_height_smooth(k), &
            v_smooth(k), &
            fd_array_smooth(k), &
            dEorb_smooth(k), &
            f_smooth(k)
         end do

         close(prof_unit)

         close(my_unit)

         !advance orbit
         call rk4_small_step(id, ierr, s%xtra(a_curr), s%time, s%dt_next, a_out, R_smooth, f_smooth, n_smooth, v_smooth, fd_array_smooth, dEorb_smooth, m_smooth, rho_boundary_smooth)

      end subroutine smooth_profile_holgado


      subroutine bronner_prescription(id, ierr, M, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k, azone_temp
         real(dp), intent(in) :: M
         real(dp), intent(out) :: a_out
         real(dp) :: a_in, f1, f2, f3, mu1, mu2, mu3, mu4, eps_rho_br, mach_br

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is just parameters for mr15 part?? eventho its meant to be given by the kim2010
         f1 = 1.91791946
         f2 = -1.52814698
         f3 = 0.75992092

         mu1 = -2.14034214
         mu2 = 1.94694764
         mu3 = 1.19007536
         mu4 = 1.05762477

         a_in = SQRT(( s%xtra(a1_curr) - s%xtra(c1_curr) )**2 + (s%xtra(b1_curr) - s%xtra(d1_curr))**2)
         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         if (a_in > 100) then
            v_rel = SQRT( (s%xtra(c2_curr) - s%xtra(a2_curr) + (s%xtra(omega_env)*(s%xtra(d1_curr) - s%xtra(b1_curr))) )**2 + (s%xtra(d2_curr) - s%xtra(b2_curr) - (s%xtra(omega_env)*(s%xtra(c1_curr) - s%xtra(a1_curr))) )**2 )
         else 
            v_rel = SQRT( (s%xtra(c2_curr) - s%xtra(a2_curr) + (s%xtra(omega_env)*(s%xtra(d1_curr) - s%xtra(b1_curr))) - (s%u(azone)*(s%xtra(c1_curr) - s%xtra(a1_curr)))/a_in )**2 + (s%xtra(d2_curr) - s%xtra(b2_curr) - (s%xtra(omega_env)*(s%xtra(c1_curr) - s%xtra(a1_curr))) - (s%u(azone)*(s%xtra(d1_curr) - s%xtra(b1_curr)))/a_in )**2 )
         end if

         Ra_br = 2*s% cgrav(azone_temp)*M/(v_rel**2)                 !cm  !different than the one in bronner
         mdot_hl_br = pi*(Ra_br**2)*(s%rho(azone_temp))*v_rel           !g/s
         fd_hl_br = mdot_hl_br*v_rel                           !dyne (g cm s^-2)
         eps_rho_br = Ra_br/s%scale_height(azone_temp)                                                !dimensionless
         mdot_mr15_ratio_br = (10)**(mu1 + mu2/(1 + mu3*eps_rho_br + mu4*(eps_rho_br**2)))   !dimensionless
         mdot_mr15_br = mdot_mr15_ratio_br*mdot_hl_br                                        !g/s
         fd_mr15_ratio_br = f1 + f2*eps_rho_br + f3*(eps_rho_br**2)
         fd_mr15_br = fd_mr15_ratio_br*fd_hl_br                                              !dyne (g cm s^-2)

         call edd_and_hyper(id, ierr, M, s%xtra(a_curr))            !all cgs

         if (eta*mdot_mr15_br >= mdot_hyper .or. eta*mdot_mr15_br <= mdot_edd) then
            mdot_br = eta*mdot_mr15_br       !g/s
         else
            mdot_br = mdot_edd           !g/s
         end if

         mach_br = v_rel / s% csound(azone_temp)

         call kim2010(id, ierr, mach_br, a_in, Ra_br, azone_temp, v_rel, M, fd_br)  !all cgs

         edot_br = fd_br*v_rel                !erg/s

         call orbital_evolution_bronner(id, ierr, s%time, s%dt_next, a_in, azone_temp, fd_br, s%xtra(omega_env), M, v_rel, a_out)  !cm
      end subroutine bronner_prescription

      subroutine orbital_evolution_bronner(id, ierr, t, dt, a_in, azone_in, fd_in, omega_in, M_in, v_in, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp), intent(in) :: t, dt, a_in, fd_in, omega_in, M_in, v_in
         integer, intent(in) :: azone_in
         real(dp), intent(out) :: a_out
         integer :: azone_temp

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         a1_next = s%xtra(a1_curr) + s%xtra(a2_curr)*dt
         b1_next = s%xtra(b1_curr) + s%xtra(b2_curr)*dt
         c1_next = s%xtra(c1_curr) + s%xtra(c2_curr)*dt
         d1_next = s%xtra(d1_curr) + s%xtra(d2_curr)*dt

         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         a2_next = a2_curr + dt*(s%cgrav(azone_temp)*s%xtra(M_ns_curr)*(s%xtra(c1_curr) - s%xtra(a1_curr))/(a_in**3))
         b2_next = b2_curr + dt*(s%cgrav(azone_temp)*s%xtra(M_ns_curr)*(s%xtra(d1_curr) - s%xtra(b1_curr))/(a_in**3))
         if (a_in > 100*Rsun) then
            c2_next = c2_curr + dt*(s%cgrav(azone_temp)*s%m(azone_in)*(s%xtra(c1_curr) - s%xtra(a1_curr))/(a_in**3) + fd_in*(s%xtra(c2_curr) - s%xtra(a2_curr) + omega_in*(s%xtra(d1_curr) - s%xtra(b1_curr)))/(M_in*v_in))
            d2_next = d2_curr + dt*(s%cgrav(azone_temp)*s%m(azone_in)*(s%xtra(d1_curr) - s%xtra(b1_curr))/(a_in**3) + fd_in*(s%xtra(d2_curr) - s%xtra(b2_curr) - omega_in*(s%xtra(c1_curr) - s%xtra(a1_curr)))/(M_in*v_in))
         else
            c2_next = c2_curr + dt*(s%cgrav(azone_temp)*s%m(azone_in)*(s%xtra(c1_curr) - s%xtra(a1_curr))/(a_in**3) + fd_in*(s%xtra(c2_curr) - s%xtra(a2_curr) + omega_in*(s%xtra(d1_curr) - s%xtra(b1_curr)) - (s%u(azone_in)*(s%xtra(c1_curr) - s%xtra(a1_curr)))/a_in)/(M_in*v_in))
            d2_next = d2_curr + dt*(s%cgrav(azone_temp)*s%m(azone_in)*(s%xtra(d1_curr) - s%xtra(b1_curr))/(a_in**3) + fd_in*(s%xtra(d2_curr) - s%xtra(b2_curr) - omega_in*(s%xtra(c1_curr) - s%xtra(a1_curr)) - (s%u(azone_in)*(s%xtra(d1_curr) - s%xtra(b1_curr)))/a_in)/(M_in*v_in))
         end if

         a_out = SQRT( (a1_next - c1_next)**2 + (b1_next - d1_next)**2 )

      end subroutine orbital_evolution_bronner

      subroutine kim2010(id, ierr, mach_in, a_in, Ra_in, azone_in, v_in, M_in, fd_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k
         real(dp), intent(in) :: M_in, mach_in, a_in, Ra_in, v_in
         integer, intent(in) :: azone_in
         real(dp), intent(out) :: fd_out
         real(dp) :: i_var, beta, eta_b, cd
         integer :: azone_temp

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         beta = s%cgrav(azone_temp)*M_in/(a_in*s%csound(azone_in)**2)   !dimensionless
         eta_b = beta/(mach_in**2 - 1)                          !dimensionless
         cd = 0.002

         if (mach_in .lt. 1.01) then
            i_var = 0.7706*LOG((1+mach_in)/(1.0004 - 0.9185*mach_in)) - 1.473*mach_in
         else if (mach_in >= 1.01 .and. mach_in .lt. 4.4) then
            i_var = LOG(330*(a_in*(mach_in-0.71)**5.72)/(Ra_in*mach_in**9.58))
         else
            i_var = LOG(a_in / (Ra_in*(0.11*mach_in + 1.65)))
         end if

         if (eta_b > 0.1 .and. mach_in > 1.01) then
            fd_out = -cd*0.7*4*pi*s%rho(azone_in)*(1 + 0.46*(beta**1.1)/(mach_in**2 - 1)**0.11)*(s%cgrav(azone_temp)*(M_in**2))/((v_in**2)*(eta_b**0.5))
         else
            fd_out = -cd*4*pi*s%rho(azone_in)*(s%cgrav(azone_temp)*(M_in**2))*i_var/(v_in**2)
         end if
      end subroutine kim2010

      subroutine hl_and_energy(id, ierr, M, a_in)
         ! star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         ! subroutine variables
         integer :: k
         integer :: j
         real(dp) :: dm
         real(dp), intent(in) :: M
         integer :: azone_temp
         real(dp), intent(in) :: a_in
         real(dp) :: E_outer

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return


         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         k = azone_temp

         ! mass is in g, radius is in cm, time is in s
         Ra = 2*s%cgrav(k)*M/(v**2 + s%csound(k)**2 )                 !cm
         mdot_hl = pi*(Ra**2)*(s%rho_face(k))*v           !g/s
         fd_hl = mdot_hl*v                           !dyne (g cm s^-2)
         edot_hl = fd_hl*v                          !erg/s

         ! try to skip the allocate deallocate part
         ! see how the sum is working and how we can only take the sum of the zones we need
         
         E_outer = 0.0d0

         ! loop only runs upto azone_temp
         do j = 2, k
            dm = s%m(j-1) - s%m(j)   ! positive shell mass
            E_outer = E_outer + dm / s%R(j-1)  ! outer radius
         end do

         ! compute total orbital energy including outer shell contribution
         Eorb =  s%cgrav(k)*M*(s%m(k))/(2.0d0*s%R(k))
         Eorb =  Eorb + s%cgrav(k)*M*E_outer  ! erg

         dEorb = (s%cgrav(k)*M/(2*s%R(k)))*((s%m(k)/s%R(k)) + 4*pi*(s%R(k)**2)*s%rho_face(k)) 

         
      end subroutine hl_and_energy


      subroutine mr15(id, ierr, a_in)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k, azone_temp
         real(dp) :: f1, f2, f3, mu1, mu2, mu3, mu4
         real(dp), intent(in) :: a_in

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         f1 = 1.91791946
         f2 = -1.52814698
         f3 = 0.75992092

         mu1 = -2.14034214
         mu2 = 1.94694764
         mu3 = 1.19007536
         mu4 = 1.05762477

         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         k=azone_temp

         eps_rho = Ra/s%scale_height(k)                                                !dimensionless
         fd_mr15_ratio = f1 + f2*eps_rho + f3*(eps_rho**2)                          !dimensionless
         mdot_mr15_ratio = (10)**(mu1 + mu2/(1 + mu3*eps_rho + mu4*(eps_rho**2)))   !dimensionless
         mdot_mr15 = mdot_mr15_ratio*mdot_hl                                        !g/s
         fd_mr15 = fd_mr15_ratio*fd_hl                                              !dyne (g cm s^-2)
      end subroutine mr15

      subroutine edd_and_hyper(id, ierr, M, a_in)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp), intent(in) :: M
         real(dp), intent(in) :: a_in

         !subroutine variables
         integer :: k, azone_temp

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !call find_zone(s%R(1:s%nz), a_in, azone_temp)

         !mass is in g, radius is in cm, time is in s
         mdot_edd = 3.5*(1d-8)*(M/(1.33*Msun))*(0.34/op_const)*Msun/secyer      !gm/s
         mdot_hyper = 8.9*(1d-5)*((op_const/0.34)**(-0.73))*Msun/secyer                         !gm/s

      end subroutine edd_and_hyper


      ! isnt complete   
      subroutine kim2007(id, ierr, n, M)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, n
         real(dp) :: mach(n)
         real(dp), intent(in) :: M

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return 

         
         do k = 1, n
            mach(k) = v / s% csound(k)  
            if (mach(k) .lt. 1.0) then
               fd_arr = 0.7706*LOG((1+mach(k))/(1.0004 - 0.9185*mach(k))) - 1.473*mach(k)
            else if (mach(k) >= 1.0 .and. mach(k) .lt. 4.4) then
               fd_arr = LOG(330*(s%R(k)*(mach(k)-0.71)**5.72)/(R0*mach(k)**9.58))
            else
               fd_arr = LOG(s%R(k) / (R0*(0.11*mach(k) + 1.65)))
            end if
            fd_arr = fd_arr*(4*pi*s%rho(k)*((s%cgrav(k)*M)**2)/(v**2))
         end do
      end subroutine kim2007

      !DONE: to evaluate the function F for the differential equation da/dt = F(a)
      subroutine orbital_evolution_function_holgado(id, ierr, t, a_in, fa) 
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         real(dp), intent(in) :: a_in, t
         real(dp), intent(out) :: fa
         integer :: k, azone_temp
         
         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return 
         
         ! print *, 'orbital evolution function called'

         !finding the zone where the function is to be evaluated
         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         !evaluating the function
         fa = f   ! no loop now

      end subroutine orbital_evolution_function_holgado

      subroutine energy_injection_bronner(id, ierr, edot_in, nr_in, rzones_in, a_in, r1_in, r2_in)

         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer, intent(in) :: nr_in, rzones_in(:)
         real(dp), intent(in) :: a_in, r1_in, r2_in, edot_in
         real(dp) :: temp_heat
         real(dp) :: sigma
         integer :: i, k

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         print *, 'edot_in = ', edot_in

         temp_heat = 0d0
         kernel(1:s%nz) = 0d0
         
         ! CASE 1 : Companion inside envelope

         if (a_in >= r1_in .and. a_in <= r2_in) then

            do i = 1, nr_in
               k = rzones_in(i)
               ! inner side
               if (s%R(k) >= s%R(rzones_in(nr_in)) .and. s%R(k) < a_in) then
                  sigma = a_in - s%R(rzones_in(nr_in))
                  kernel(k) = exp(-((a_in - s%R(k))/sigma)**2)
               ! outer side
               else if (s%R(k) >= a_in .and. s%R(k) <= s%R(rzones_in(1))) then
                  sigma = s%R(rzones_in(1)) - a_in
                  kernel(k) = exp(-((s%R(k) - a_in)/sigma)**2)
               end if
               temp_heat = temp_heat + kernel(k)*s%dm(k)
            end do

   
         ! CASE 2 : Companion outside envelope
         ! Inject using inward half-Gaussian from outer edge
   
         else if (a_in > r2_in) then

            sigma = r2_in - r1_in
            do i = 1, nr_in
               k = rzones_in(i)
               ! Half-Gaussian centered at outer boundary r2_in
               if (s%R(k) >= r1_in .and. s%R(k) <= r2_in) then
                  kernel(k) = exp(-((r2_in - s%R(k))/sigma)**2)
               end if
               temp_heat = temp_heat + kernel(k)*s%dm(k)
            end do
         end if

         ! Normalize heating
         s%extra_heat(1:s%nz)%val = efactor * kernel(1:s%nz) * edot_in / temp_heat

         print *, 'nr_in = ', nr_in
         print *, 'temp_heat = ', temp_heat
         print *, 'extra_heat = ', SUM(s%extra_heat(1:s%nz)%val)

      end subroutine energy_injection_bronner

      subroutine binding_energy(id, ierr, zone_in, ebind_out)
         !star variables
         integer, intent(in) :: id, zone_in
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp), intent(out) :: ebind_out      !eorb_change_out
         integer :: i

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! print *, 'binding energy called'

         ebind_out = 0
         do i = 1, zone_in
            ebind_out = ebind_out + (s%energy(i) - s%cgrav(i)*s%m(i)/s%R(i))*s%dm(i)
         end do
         !eorb_change_out = Eorb(zone_in) - Eorb(1)

      end subroutine binding_energy

      !DONE: To find the equatorial radius and beta secular
      subroutine Req_and_beta(e1, Req1, Rbar1, beta1)
         real(dp), intent(in) :: e1
         real(dp), intent(out) :: Req1, Rbar1, beta1
         
         beta1 = 3*(1 - ((e1*sqrt(1-e1**2))/(asin(e1))))/(2*e1**2) - 1                             !dimensionless
         Rbar1 = R0*((asin(e1) * ((1-e1**2)**(1/6)) * (1-beta1)) / e1)**(-1 * n_poly / (3-n_poly)) !cm
         Req1 = Rbar1/((1-e1**2)**(1/6))                                                           !cm
      end subroutine Req_and_beta

      !DONE: To evaluate the omega function and solve for the inverse
      subroutine omega_function(e1, M, omega1, cgrav)
         !subroutine variables
         real(dp), intent(in) :: e1, M
         real(dp), intent(in) :: cgrav
         real(dp), intent(out) :: omega1
         real(dp) :: rho_bar, qn

         rho_bar = 3*M/(4*pi*(R0**3))      !gm/cm^3
         qn = (1-n_poly/5)                 !dimensionless

         
         omega1 = sqrt(2*pi*cgrav*rho_bar*( (sqrt(1-e1**2)*(3-2*e1**2)*asin(e1)/(e1**3)) - 3*(1-e1**2)/(e1**2) )/qn) !Hz
      end subroutine omega_function

      !DONE: To solve for the inverse of the omega function
      subroutine omega_func_solve_inverse(e_in, omega_val, tol, M, e_out, cgrav)
         !subroutine variables
         real(dp), intent(in) :: e_in, omega_val, tol, M
         real(dp), intent(out) :: e_out
         real(dp), intent(in) :: cgrav

         !local variables
         real(dp) :: e1, omega1
         integer :: iterations

         iterations = 0

         e1 = e_in !dimensionless
         call omega_function(e1, M, omega1, cgrav)

         do while (abs(omega1 - omega_val) >= tol .and. iterations < 100)
            e1 = e1*omega_val/omega1
            call omega_function(e1, M, omega1, cgrav)
            iterations = iterations + 1
         end do

         e_out = e1 !dimensionless
      end subroutine omega_func_solve_inverse

      !DONE: to evaluate the spin evolution and quadrupole moment evolution
      subroutine omega_and_q(id, ierr, M, omega_function , cgrav)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp), intent(in) :: M, omega_function, cgrav
         integer :: i


         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         Qmax_next = 4*1d39*(M_acc_next/M_crust)             !g cm^2
                           
         if (s%x_ctrl(14) /= s%x_ctrl(15)) then
            call forward_euler(s%xtra(omega_curr), omega_function, s%dt_next, omega_next)  !Hz
         else 
            call forward_euler(s%xtra(omega_curr), omega_function, 0d0, omega_next)          !Hz
         end if

         call omega_func_solve_inverse(s%xtra(e_curr), omega_next, 1d-9, M, e_next, cgrav)

         aa_next = R0*(1 + e_next/2)                        !cm
         bb_next = R0*(1 - e_next/2)                        !cm
         mom_inert_next = M*(aa_next**2 + bb_next**2)/5     !gm cm^2

         call Req_and_beta(e_next, Req, Rbar, beta)
         Qtb_next = sqrt((5*(clight**5)*mdot*sqrt(cgrav*M*Req))/(32*cgrav*(omega_next**5))) !g cm^2
         if (Qtb_next > Qmax_next) then
            Q_next = Qmax_next*exp(-1*s%time*mdot/M_crust)   !g cm^2
            decay_coeff = exp(-1*s%time*mdot/M_crust)        !dimensionless
         else
            Q_next = Qtb_next !g cm^2
            decay_coeff = exp(-1*s%time*mdot/M_crust) !dimensionless
         end if
      end subroutine omega_and_q

      subroutine inject_energy(id, ierr)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp) :: temp, temp_heat
         integer, allocatable :: rzones(:)
         integer :: i, j, nr, azone_temp, zone_curr
         integer :: n_smooth, n_refine
         integer :: n_position, n_power, n_velocity

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! INITIALISING AND ALLOCATING
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         n_refine = 10
         n_smooth = (s%nz-1)*n_refine + 1

         n_position = 732727
         n_power = 1278
         n_velocity = 1278
         
         s%x_ctrl(14) = s%model_number

         !Initialising
         M_ns_initial = s%x_ctrl(1)*Msun          !mass of the NS (gm)
         op_const = s%x_ctrl(2)                   !opacity constant (cm^2/g)
         eta = s%x_ctrl(3)                        !efficiency factor
         efactor = s%x_ctrl(4)                    !multiplication factor for injected energy
         M_crust = s%x_ctrl(5)*Msun               !mass of the crust (gm)
         omega_initial = 2*pi*s%x_ctrl(6)         !initial spin frequency (Hz)
         R0 = s%x_ctrl(7)                         !equatorial radius (cm)
         e_initial = s%x_ctrl(8)                  !initial ellipticity
         n_poly = s%x_ctrl(9)                     !polytropic index
         beta_sec = s%x_ctrl(10)                  !beta secular
         D = s%x_ctrl(11)*1d3*pc                  !distance to the source (cm)
         omega_env_factor = s%x_ctrl(12)          !omega_env_factor*initial_orbital_period = envelope omega
         prescription = s%x_ctrl(13)              !general prescription for the CE evolution

         !allocating relevant variables
         !if (allocated(fd_arr)) then 
            !deallocate(fd_arr, mdot_arr, f, edot)
            !deallocate(mdot_hl, fd_hl, edot_hl, v, Ra, Eorb, dEorb)
            !deallocate(mdot_mr15, fd_mr15, mdot_mr15_ratio, fd_mr15_ratio, eps_rho)
            !deallocate(mdot_edd, mdot_hyper)!, omega_arr, e_arr, beta_arr)!, rand)
            !deallocate(R_smooth, m_smooth, rho_boundary_smooth ,scale_height_smooth, v_smooth)
            !deallocate(Ra_smooth, eps_rho_smooth, mr15_ratio_smooth, mdot_hl_smooth, fd_hl_smooth)
            !deallocate(fd_mr15_smooth, fd_array_smooth, dEorb_smooth, f_smooth, edot_smooth)
         !end if
         !allocate(fd_arr(s%nz), mdot_arr(s%nz), f(s%nz), edot(s%nz))
         !allocate(mdot_hl(s% nz), fd_hl(s% nz), edot_hl(s% nz), v(s% nz), Ra(s% nz), Eorb(s%nz), dEorb(s%nz))
         !allocate(mdot_mr15(s%nz), fd_mr15(s%nz), mdot_mr15_ratio(s%nz), fd_mr15_ratio(s%nz), eps_rho(s%nz))
         !allocate(mdot_edd(s%nz), mdot_hyper(s%nz))
         !allocate(R_smooth(n_smooth), m_smooth(n_smooth), rho_boundary_smooth(n_smooth),scale_height_smooth(n_smooth), v_smooth(n_smooth))
         !allocate(Ra_smooth(n_smooth), eps_rho_smooth(n_smooth), mr15_ratio_smooth(n_smooth), mdot_hl_smooth(n_smooth), fd_hl_smooth(n_smooth))
         !allocate(fd_mr15_smooth(n_smooth), fd_array_smooth(n_smooth), dEorb_smooth(n_smooth), f_smooth(n_smooth), edot_smooth(n_smooth))
         if (.not. allocated(t_position)) then
            allocate(t_position(n_position))
         end if
         if (.not. allocated(r_position)) then
            allocate(r_position(n_position))
         end if
         if (.not. allocated(t_power)) then
            allocate(t_power(n_power))
         end if
         if (.not. allocated(f_power)) then
            allocate(f_power(n_power))
         end if
         if (.not. allocated(t_velocity)) then
            allocate(t_velocity(n_velocity))
         end if
         if (.not. allocated(f_velocity)) then
            allocate(f_velocity(n_velocity))
         end if
         if (allocated(kernel)) then
            deallocate(kernel)
         end if
         allocate(kernel(s%nz))

         if (s%model_number == 1) then
            s%xtra(a_curr) = 124*Rsun !9.8d-1*s%R(1)          !cm
            s%xtra(M_ns_curr) = M_ns_initial                  !gm
            s%xtra(M_acc_curr) = 0                            !gm   
            s%xtra(omega_curr) = omega_initial                !spin frequency (Hz)
            s%xtra(e_curr) = e_initial                        !ellipticity
            s%xtra(Qmax_curr) = 4*1d39*(s%xtra(M_acc_curr)/M_crust)   !gm cm^2
            s%xtra(aa_curr) = R0*(1+s%xtra(e_curr)/2)                 !cm
            s%xtra(bb_curr) = R0*(1-s%xtra(e_curr)/2)                 !cm
            s%xtra(mom_inert_curr) = s%xtra(M_ns_curr)*(s%xtra(aa_curr)**2 + s%xtra(bb_curr)**2)/5                 !gm cm^2
            
            call find_zone(s%R(1:s%nz), s%xtra(a_curr), azone_temp)
 
            s%xtra(omega_env) = omega_env_factor*SQRT(s%cgrav(azone_temp)*(s%xtra(M_ns_curr) + s%m(azone_temp))/(s%R(azone_temp)**3))
            if (prescription == 2) then
               s%xtra(a1_curr) = 0
               s%xtra(b1_curr) = 0
               s%xtra(c1_curr) = s%R(azone_temp)
               s%xtra(d1_curr) = 0
               s%xtra(a2_curr) = 0
               s%xtra(b2_curr) = 0
               s%xtra(c2_curr) = 0
               s%xtra(d2_curr) = SQRT(s%cgrav(azone_temp)*(s%xtra(M_ns_curr) + s%m(azone_temp))/s%R(azone_temp))
            end if
            s%x_ctrl(15) = 2
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! EVALUATING CURRENT PARAMETERS
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !doing the orbital evolution
         if (prescription == 1) then  
            if (s%x_ctrl(19) == 1) then                         ! smoothing turned on
               call smooth_profile_holgado(id, ierr, s%xtra(M_ns_curr), s%xtra(a_curr), a_next)
            else if (s%x_ctrl(19) == 0) then                    ! smothing disabled, call normal holgado_prescription
               call holgado_prescription(id, ierr, s%xtra(M_ns_curr), s%xtra(a_curr), a_next)
            else if (s%x_ctrl(19)== 2) then
               call holgado_simulation_force(id, ierr, s%xtra(M_ns_curr), a_next)
            end if
         else if (prescription == 2) then
            call bronner_prescription(id, ierr, s%xtra(M_ns_curr), a_next)
         end if

         !injecting energy into the envelope
         call find_azone_r1r2(s%R(1:s%nz), s%xtra(a_curr), r1, r2, azone, rzones, nr)
         
         if (s%x_ctrl(16) == 1) then    !energy injection enabled
            if (prescription == 1) then
               call energy_injection_bronner(id, ierr, -edot, nr, rzones, s%xtra(a_curr), r1, r2) 
            else if (prescription == 2) then
               call energy_injection_bronner(id, ierr, -edot_br, nr, rzones, s%xtra(a_curr), r1, r2)
            end if
         else if (s%x_ctrl(16) == 0) then       !energy injection disabled
            s%extra_heat(:)%val = 0.0d0
         end if


         !calculating the binding energy of the envelope(not used in evaluating anything else)
         call binding_energy(id, ierr, azone, ebind)                                                 !ergs

         ! accretion switch
         if (s%x_ctrl(18) == 1) then ! accretion on

            if (prescription == 1) then
               mdot = mdot_arr
            else if (prescription == 2) then
               mdot = mdot_br
            end if 

            ! update mass here
            if (s%x_ctrl(14) /= s%x_ctrl(15)) then
               call forward_euler(s%xtra(M_acc_curr), mdot, s%dt_next, M_acc_next)     !gm
               call forward_euler(s%xtra(M_ns_curr), mdot, s%dt_next, M_ns_next)       !gm
            else 
               call forward_euler(s%xtra(M_acc_curr), mdot, 0d0, M_acc_next)     !gm
               call forward_euler(s%xtra(M_ns_curr), mdot, 0d0, M_ns_next)       !gm
            end if
         else if (s%x_ctrl(18) == 0) then ! accretion disabled
            mdot = 0d0
            M_acc_next = s%xtra(M_acc_curr)
            M_ns_next = s%xtra(M_ns_curr)
         end if 

         ! to disable gravitational waves, built the switch inside this routine
         call evaluate_strain(id, ierr)

         print *, 'model number = ', s% model_number
         print *, 'orbital separation = ', s%xtra(a_curr)/Rsun
         print *, 'omega = ', s%xtra(omega_curr)
         print *, 'Q = ', s%xtra(Q_curr)
         print *, 'strain = ', s%xtra(strain_curr)
         if (prescription == 1) then
            !print *, 'mach = ', v/s% csound(azone)
            print *, 'edot = ', edot
         else if (prescription == 2) then
            print *, 'mach = ', v_rel/s% csound(azone)
            print *, 'edot = ', -edot_br
            ! print *, 'fd_br = ', fd_br
            print *, 'a1_curr = ', s%xtra(a1_curr)
            print *, 'b1_curr = ', s%xtra(b1_curr)
            print *, 'c1_curr = ', s%xtra(c1_curr)
            print *, 'd1_curr = ', s%xtra(d1_curr)
            print *, 'a2_curr = ', s%xtra(a2_curr)
            print *, 'b2_curr = ', s%xtra(b2_curr)
            print *, 'c2_curr = ', s%xtra(c2_curr)
            print *, 'd2_curr = ', s%xtra(d2_curr)
            print *, 'total injected energy = ', efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz))*s%dt_next
         end if

         print *, '##############################################################################'
         print *, 'ONE TIMESTEP DONE'
         print *, '##############################################################################'


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! SAVING RELEVANT VARIABLES
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         s%xtra(a_curr) = a_next
         s%xtra(M_acc_curr) = M_acc_next
         s%xtra(M_ns_curr) = M_ns_next
         s%xtra(omega_curr) = omega_next
         s%xtra(e_curr) = e_next
         s%xtra(Qtb_curr) = Qtb_next
         s%xtra(Qmax_curr) = Qmax_next
         s%xtra(Q_curr) = Q_next
         s%xtra(aa_curr) = aa_next
         s%xtra(bb_curr) = bb_next
         s%xtra(mom_inert_curr) = mom_inert_next
         s%xtra(strain_curr) = strain_next

         s%xtra(a1_curr) = a1_next
         s%xtra(b1_curr) = b1_next
         s%xtra(c1_curr) = c1_next
         s%xtra(d1_curr) = d1_next
         s%xtra(a2_curr) = a2_next
         s%xtra(b2_curr) = b2_next
         s%xtra(c2_curr) = c2_next
         s%xtra(d2_curr) = d2_next

         s%x_ctrl(15) = s%model_number

      end subroutine inject_energy

      subroutine evaluate_strain(id, ierr)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp) :: omega_func_curr
         integer :: azone_temp
         real(dp) :: cgrav

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         call find_zone(s%R(1:s%nz), s%xtra(a_curr), azone_temp)

         cgrav = s%cgrav(azone_temp)

         if (s%x_ctrl(17) == 1) then   ! gravitational waves enabled
            if (s%model_number == 1) then
               call Req_and_beta(s%xtra(e_curr), Req, Rbar, beta)
               omega_next = s%xtra(omega_curr)                             !Hz
               s%xtra(Qtb_curr) = sqrt((5*(clight**5)*mdot*sqrt(cgrav*s%xtra(M_ns_curr)*Req))/(32*cgrav*((s%xtra(omega_curr))**5))) !g cm^2
               Qtb_next = s%xtra(Qtb_curr)                                 !g cm^2
               Qmax_next = s%xtra(Qmax_curr)                               !g cm^2
               s%xtra(Q_curr) = s%xtra(Qtb_curr)                           !g cm^2
               Q_next = s%xtra(Q_curr)                                     !g cm^2
               e_next = s%xtra(e_curr)                                     !dimensionless
               aa_next = s%xtra(aa_curr)                                   !cm
               bb_next = s%xtra(bb_curr)                                   !cm
               mom_inert_next = s%xtra(mom_inert_curr)                     !gm cm^2

               decay_coeff = 1                                             !dimensionless

               s%xtra(strain_curr) = 2*cgrav*((s%xtra(omega_curr))**2)*(s%xtra(Q_curr))/(D*(clight**4)) !dimensionless
               strain_next = s%xtra(strain_curr)    !dimensionless

            else 
               if (s%xtra(e_curr) > 0.817) then 
                  s%xtra(e_curr) = 0.817                               !dimensionless
               end if

               call Req_and_beta(s%xtra(e_curr), Req, Rbar, beta)
               omega_func_curr = (mdot*sqrt(cgrav*s%xtra(M_ns_curr)*Req) - (32*cgrav*(s%xtra(omega_curr)**5)*(s%xtra(Q_curr)**2))/(5*clight**5))/s%xtra(mom_inert_curr)
               call omega_and_q(id, ierr, M_ns_next, omega_func_curr, cgrav)
               strain_next = 2*cgrav*(omega_next**2)*Q_next/(D*(clight**4))  !dimensionless

            end if

         else if (s%x_ctrl(17) == 0) then ! gravitaional waves disabled
            Req = 0d0
            Rbar = 0d0
            beta = 0d0
            Qtb_next = 0d0
            Q_next = 0d0
            Qmax_next = 0d0
            strain_next = 0d0
            ! Keep geometrical / structural quantities unchanged
            e_next = s%xtra(e_curr)
            aa_next = s%xtra(aa_curr)
            bb_next = s%xtra(bb_curr)
            mom_inert_next = s%xtra(mom_inert_curr)
            omega_next = s%xtra(omega_curr)

         end if


      end subroutine evaluate_strain


      subroutine spline(x, y, n, yp1, ypn, y2)
         use const_def
         implicit none

         integer, intent(in) :: n
         real(dp), dimension(n), intent(in) :: x, y
         real(dp), intent(in) :: yp1, ypn
         real(dp), dimension(n), intent(out) :: y2
         integer :: i
         real(dp) :: p, sig

         ! Initialize output
         y2 = 0.0_dp

         ! Check for too small arrays
         if (n < 2) then
            return
         end if

         ! First boundary
         if (yp1 > 0.99d30) then
            y2(1) = 0.0_dp
         else
            y2(1) = -0.5_dp
         end if

         ! Loop over interior points safely
         do i = 2, n - 1
            if (abs(x(i+1) - x(i-1)) < 1.0d-30) then
               sig = 0.5d0
            else
               sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
            end if
            p = sig * y2(i-1) + 2.0d0
            y2(i) = (sig - 1.0d0) / p
         end do

         ! Last boundary
         if (ypn > 0.99d30) then
            y2(n) = 0.0_dp
         else
            y2(n) = -0.5_dp
         end if
      end subroutine spline


      subroutine splint(x, y, y2, n, xint, yint)
      use const_def
      implicit none
      
      integer, intent(in) :: n
      real(dp), dimension(n), intent(in) :: x, y, y2
      real(dp), intent(in) :: xint
      real(dp), intent(out) :: yint
      integer :: klo, khi, k
      real(dp) :: h, b, a

      klo = 1
      khi = n
      do while (khi - klo > 1)
         k = (khi + klo) / 2
         if (x(k) > xint) then
            khi = k
         else
            klo = k
         end if
      end do

      h = x(khi) - x(klo)
      if (h == 0.0d0) then
         print *, 'Error: xint is out of range.'
         stop
      end if

      a = (x(khi) - xint) / h
      b = (xint - x(klo)) / h
      yint = a * y(klo) + b * y(khi) + ((a**3 - a) * y2(klo) + (b**3 - b) * y2(khi)) * (h**2) / 6.0d0
      end subroutine splint


         
      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% other_energy => inject_energy
         s% other_cgrav => other_cgrav
         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns  

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls

      


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

      end subroutine extras_startup



      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if
         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         if (s% xtra(a_curr) .lt. s%R(s%nz)) then
            extras_check_model = terminate
            s% termination_code = t_xtra1
            termination_code_str(t_xtra1) = 'orbital separation less than core size'
            return
         end if

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (prescription == 1) then
            if (s%x_ctrl(19) == 1) then
               how_many_extra_history_columns = 37   ! smoothing on, the smooth_profile and dEorb.data act as profiles and history
            else if (s%x_ctrl(19) == 0) then
               how_many_extra_history_columns = 43    ! smoothing off, normal holgado
            else if (s%x_ctrl(19) == 2) then   
               how_many_extra_history_columns = 43    !smoothing on, holgado with 3d force data
            end if 
         else
            how_many_extra_history_columns = 43
         end if
      end function how_many_extra_history_columns
      
      
      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         integer :: k, zone_temp
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         

         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.

         names(1) = 'a_curr'
         names(2) = 'Injected_E_per_timestep'
         names(3) = 'mdot_ns'
         names(4) = 'azone'
         names(5) = 'r1'
         names(6) = 'r2'
         names(7) = 'M_ns'
         names(8) = 'M_acc'
         names(9) = '-ebind_curr'
         names(10) = 'Qmax'
         names(11) = 'Qtb'
         names(12) = 'Q'
         names(13) = 'omega'
         names(14) = 'Req'
         names(15) = 'Rbar'
         names(16) = 'beta'
         names(17) = 'e'
         names(18) = 'aa'
         names(19) = 'bb'
         names(20) = 'mom_inert'
         names(21) = 'strain'
         names(22) = 'mass_at_a'
         names(23) = 'R_azone'
         names(24) = 'decay_coeff'
         names(25) = 'v_esc'
         names(26) = 'u_k'
         names(27) = 'csound_ns'
         names(28) = 'lnT_ns'
         names(29) = 'rho_ns'
         names(30) = 'lnPgas_ns'
         names(31) = 'gamma1_ns'
         names(32) = 'scale_height_ns'
         names(33) = 'R99'
         names(34) = 'R95'
         names(35) = 'R90'
         names(36) = 'mdot_edd'
         if (s%x_ctrl(19) == 0) then  ! smoothing disabled
            names(37) = 'v_ns'
            names(38) = 'Ra'
            names(39) = 'mdot_hl'
            names(40) = 'mdot_MR15'
            names(41) = 'fd_ns'
            names(42) = 'fd_hl'
            names(43) = 'fd_MR15'
         end if
         if (s%x_ctrl(19) == 2) then  ! reading 3d force data
            names(37) = 'v_ns'
            names(38) = 'Ra'
            names(39) = 'mdot_hl'
            names(40) = 'mdot_MR15'
            names(41) = 'fd_ns'
            names(42) = 'fd_hl'
            names(43) = 'fd_MR15'
         end if

         vals(1) = s%xtra(a_curr)
         vals(2) = efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz))*s%dt_next
         vals(3) = mdot
         vals(4) = azone
         vals(5) = r1
         vals(6) = r2
         vals(7) = s%xtra(M_ns_curr)
         vals(8) = s%xtra(M_acc_curr)
         vals(9) = -ebind
         vals(10) = s%xtra(Qmax_curr)
         vals(11) = s%xtra(Qtb_curr)
         vals(12) = s%xtra(Q_curr)
         vals(13) = s%xtra(omega_curr)
         vals(14) = Req
         vals(15) = Rbar
         vals(16) = beta
         vals(17) = s%xtra(e_curr)
         vals(18) = s%xtra(aa_curr)
         vals(19) = s%xtra(bb_curr)
         vals(20) = s%xtra(mom_inert_curr)
         vals(21) = s%xtra(strain_curr)
         vals(22) = s%m(azone)
         vals(23) = s%R(azone)
         vals(24) = decay_coeff
         vals(25) = sqrt(2*s%cgrav(azone)*s% m(azone)/s% R(azone))
         vals(26) = s%u(azone)
         vals(27) = s%csound(azone)
         vals(28) = s%lnT(azone)
         vals(29) = s%rho(azone)
         vals(30) = s%lnPgas(azone)
         vals(31) = s%gamma1(azone)
         vals(32) = s%scale_height(azone) 
         call find_zone(s%m(1:s%nz), 0.99*s%m(1), zone_temp)
         vals(33) = s%R(zone_temp)

         call find_zone(s%m(1:s%nz), 0.95*s%m(1), zone_temp)
         vals(34) = s%R(zone_temp)

         call find_zone(s%m(1:s%nz), 0.90*s%m(1), zone_temp)
         vals(35) = s%R(zone_temp)
         vals(36) = mdot_edd

         if (s%x_ctrl(19) == 0) then  ! smoothing disabled
            if (prescription == 1) then
               vals(37) = v
               vals(38) = Ra
               vals(39) = mdot_hl
               vals(40) = mdot_mr15
               vals(41) = fd_arr
               vals(42) = fd_hl
               vals(43) = fd_mr15
            else if (prescription == 2) then
               vals(37) = v_rel
               vals(38) = Ra_br
               vals(39) = mdot_hl_br
               vals(40) = mdot_mr15_br
               vals(41) = fd_br
               vals(42) = fd_hl_br
               vals(43) = fd_mr15_br
            end if
         else if (s%x_ctrl(19) == 2) then
            vals(37) = v
            vals(38) = Ra
            vals(39) = mdot_hl
            vals(40) = mdot_mr15
            vals(41) = fd_arr
            vals(42) = fd_hl
            vals(43) = fd_mr15
         end if 

      end subroutine data_for_extra_history_columns

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (prescription == 1) then
            how_many_extra_profile_columns = 2
         else if (prescription == 2) then
            how_many_extra_profile_columns = 0
         end if
      end function how_many_extra_profile_columns
      
      
      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         if (prescription == 1) then
            names(1) = 'kernel'
            names(2) = 'edot_heat'
            !names(3) = 'mdot'
            !names(4) = 'fd_hl'
            !names(5) = 'fd_MR15'
            !names(6) = 'fd'
            !names(7) = 'edot_hl'
            !names(8) = 'edot'
            !names(9) = 'eps_rho'
            !names(10) = 'v_ns'
            !names(11) = 'Ra'
            !names(12) = 'Eorb'
            !names(13) = 'dEorb'
            !names(14) = 'f'

            do k = 1, s%nz
               vals(k,1) = kernel(k)
               vals(k,2) = s%extra_heat(k)%val
               !vals(k,3) = mdot_arr(k)
               !vals(k,4) = fd_hl(k)
               !vals(k,5) = fd_mr15(k)
               !vals(k,6) = fd_arr(k)
               !vals(k,7) = edot_hl(k)
               !vals(k,8) = edot(k)
               !vals(k,9) = eps_rho(k)
               !vals(k,10) = v(k)
               !vals(k,11) = Ra(k)
               !vals(k,12) = Eorb(k)
               !vals(k,13) = dEorb(k)
               !vals(k,14) = f(k)
            end do
         end if
      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

         names(1) = 'a_next'
         names(2) = 'mdot_ns'
         names(3) = 'azone'
         names(4) = 'r1'
         names(5) = 'r2'
         names(6) = 'M_acc'
         names(7) = 'Q'
         names(8) = 'omega'
         names(9) = 'v_ns'
         names(10) = 'fd_ns'

         vals(1) = s%xtra(a_curr)
         vals(2) = mdot
         vals(3) = azone
         vals(4) = r1
         vals(5) = r2
         vals(6) = s%xtra(M_acc_curr)
         vals(7) = s%xtra(Q_curr)
         vals(8) = s%xtra(omega_curr)

         if (prescription == 1) then 
            vals(9) = v
            vals(10) = fd_arr
         else if (prescription == 2) then
            vals(9) = v_rel
            vals(10) = fd_br
         end if

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_finish_step = keep_going

         ! to save a profile, 
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
   
         end module run_star_extras
         
   
