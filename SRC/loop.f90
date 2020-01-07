module loop

  use global_parameter
  use mpi_mod
  use parquet_ini
  use parquet_chi
  use parquet_formfactors
  use parquet_util
  use parquet_plot
  use Parquet_kernel
  use Parquet_PhiR
  use parquet_EQMatrix
  use parquet_EQOuter
  use Parquet_equation
  use parquet_selfenergy
  use parquet_SDE
  use parquet_sus
  use parquet_check
  use preprocessing
  use hdf5_wrapper
  use hdf5

  !$use omp_lib

  implicit none

contains

  subroutine execute_loop(Converged, ite)

    !boolean that determines whether one is converged
    logical, intent(inout) :: Converged
    !iteration one is at
    integer, intent(inout) :: ite

    !for timing
    integer(4) :: count_rate, count_max
    !count_launch -> counts at program launch, count_iteration -> count at
    !beginning of iteration
    integer(4) :: count_launch, count_iteration, count_start, count_end

    !write timing message to screen after each iteration
    character(len=16)  :: it_str, time_str
    character(len=50) :: iteration_message

    !check variables
    real(dp) :: chi0_check

    if(id == master) write(*, *) 'Starting Loop'

    do while ((.NOT.Converged) .AND. ite < ite_max)

      !time each iteration
      call system_clock(count_iteration, count_rate, count_max)

      if (id == master) then
        write (*, *)
        write (*, *) "Starting Iteration : ", ite
        write (*, *)
      end if

      ! --- Green's function ---

      call system_clock(count_start, count_rate, count_max)

      !calculate  Green's functon
      call calc_Gkw_Grt(ite, Grt, .false.)

      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time to calculate Greens funcion: ', &
                                                dble(count_end - count_start)/dble(count_rate)
      
      ! --- eval Bethe Salpeter ---

      call system_clock(count_start, count_rate, count_max)

      !call BS equations in k-space
      call BSE(ite)

      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time spent on bethe-salpeter is: ', &
                                          dble(count_end - count_start)/dble(count_rate)

      if(old .or. SDE_old) then
      ! --- calculate kernel functions ---
      
      call system_clock(count_start, count_rate, count_max)
      call get_kernel_function(ite)
      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time spent on calculating kernels: ', &!t_end - t_start
                                          dble(count_end - count_start)/dble(count_rate)
      end if !anything old


      if(.not. old) then
      ! --- calculate PhiR ---

      !at first iteration calculate exponential for FT
      if(ite == 0) call calculate_expSqR()
      if(ite == 0) call calculate_expqSR()
      !call calculate PhiR
      call system_clock(count_start, count_rate, count_max)
      call calc_PhiR(ite)
      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time spent on calculating PhiR: ', &!t_end - t_start
                                          dble(count_end - count_start)/dble(count_rate)
      if(ite == ite_max - 1) then
      !  if(id == 0) call invTrafo_PhiR()
      !  call write_G(ite)
      end if!

      ! --- calculate Transformation matricies

      call system_clock(count_start, count_rate, count_max)
      if(ite == 0) then

        call calc_M_phbar_to_ph()
        call calc_M_pp_to_ph()
        call calc_M_ph_to_pp()
        call calc_M_phbar_to_pp()

      end if !ite == 0
      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time spent on calculating Transformation Tensors: ', &!t_end - t_start
                                          dble(count_end - count_start)/dble(count_rate)

      end if !.not. old

      !G_s_LL = dcmplx(0.0d0)
      !G_t_LL = dcmplx(0.0d0)
      !G_d_LL = dcmplx(0.0d0)
      !G_m_LL = dcmplx(0.0d0)

      ! --- solve parquet equations ---

      call system_clock(count_start, count_rate, count_max)
      !call PA
      if(old) call solve_parquet_equation(ite)
      if(.not. old) call calculate_PAE()
      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time spent on parquet is: ', &!t_end - t_start
                                          dble(count_end - count_start)/dble(count_rate)

      ! --- calculate self-energy ---

      call system_clock(count_start, count_rate, count_max)

      if (.NOT. allocated(SigmaOld)) allocate (SigmaOld(Nt))

      if(.not. oneshot) then
        !save old sigma for convergence check
        SigmaOld = Sigma
        if(SDE_old) call self_energy(ite, Grt)
        if(.not. SDE_old) call calc_SDE(ite, Grt)
      else
        if(ite == ite_max - 1) then
          SigmaOld = Sigma
          if(SDE_old) call self_energy(ite, Grt)
          if(.not. SDE_old) call calc_SDE(ite, Grt)
        end if !ite == ite_max - 1
      end if !.not. oneshot

      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time spent on selfenergy calculation is: ', &!t_end - t_start
                                          dble(count_end - count_start)/dble(count_rate)


      ! --- check convergence ---

      call system_clock(count_start, count_rate, count_max)

      !check convergence
      Converged = conv()

      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time spent on convergence check: ', &
                                                dble(count_end - count_start)/dble(count_rate)


      ! --- produce some intermediate output ---
      if(mod(ite, 5) == 0 .or. ite == 0 .or. ite == ite_max - 1 .or. Converged) then

        if(id == master) write(*, *)

        write (it_str, '(I5.5)') ite

        if(id == master) then
          write(*, *)
          write(*, *) "Intemdiate output at iteration ", ite
        end if !id == master

        call system_clock(count_start, count_rate, count_max)
    
        !calculate  chi0
        call calc_Gkw_Grt(ite, Grt, .true.)

        !calculate susceptibility
        call calc_sus

        !calculate susceptibilities in whole BZ
        call symmetrize_sus(ite, Converged)

        if(id == master) write(*, *) 'write selfenergy' 
        call write_sigma(ite, Converged)
        if(id == master) write(*, *) 'write susceptibilities' 
        call write_chis(ite, Converged)

        call system_clock(count_end, count_rate, count_max)
        if (id == master) write (*, "(a, f12.6)") ' Intemediate-output took: ', &
                                                  dble(count_end - count_start)/dble(count_rate)
        if(id == master) write(*, *)

      end if !ite%50 == 0


      ! --- update vertices ---

      !again put irreducible vertex into G's
      G_d_LL = F_d_LL - G_d_LL
      G_m_LL = F_m_LL - G_m_LL
      G_s_LL = F_s_LL - G_s_LL
      G_t_LL = F_t_LL - G_t_LL


      ! --- finish interation ---

      !for timing iteration      
      call system_clock(count_end, count_rate, count_max)

      if (id == master) then

        write (*, *)
        !two temporary strings for formatting output appropriately
        write (it_str, '(I5.5)') ite
        write (time_str, '(f12.2)') dble(count_end - count_iteration)/dble(count_rate)
        iteration_message = "iteration "//TRIM(it_str)//" took "//TRIM(time_str)//" s"
        write (*, *) iteration_message
        write (*, *) "-------------------------------------------------------------------------------"

      end if

      ite = ite + 1

    END DO!outer while

    if(id == master) write(*, *) 'End Loop'


  end subroutine execute_loop

end module loop
