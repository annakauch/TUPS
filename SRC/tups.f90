program test_main

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
  use loop
  use hdf5_wrapper
  use hdf5

  !$use omp_lib

  implicit none

  logical       :: Converged
  integer       :: ite

  !for timing
  integer(4) :: count_rate, count_max
  !count_launch -> counts at program launch, count_iteration -> count at
  !beginning of iteration
  integer(4) :: count_launch, count_iteration, count_start, count_end

  Converged = .FALSE.

  call hdf5_init()

  call parallel_start
  master = 0

  if(id == master) then
    write(*, *) "Starting TUPS with ", ntasks, " tasks"
  end if !id == master

  !call cpu_time(time_begin)
  call system_clock(count_launch, count_rate, count_max)

  ite = 0
  Converged = .false.

  !call cpu_time(t_start)
  call system_clock(count_start, count_rate, count_max)

  ! --- calculations to be done before the start of the iteration ---
  call do_preprocessing()

  ! --- execution of actual parquet loop ---
  call execute_loop(Converged, ite)

  ! --- output + a little postprocessing ---
  ite = ite_max

  !calculate chi0
  call calc_Gkw_Grt(ite, Grt, .true.)

  !actually wrong to calculate it at this place
  !call calc_sus()

  !dont recalculate since G (= Gamma) containes the wrong vertices
  !calculate susceptibilities in whole BZ
  !call symmetrize_sus(ite, Converged)

  if(id == master) write(*, *) 'writing program output to ', file_main
  !call write_sigma(ite, Converged)
  !call write_chis(ite, Converged)
  
  if(save_F) then 
    call write_F(ite)
  end if

  call system_clock(count_end, count_rate, count_max)
  if (id == master) write (*, "(a, f12.6)") 'TUPS took:', &!t_end - time_begin
                                                  dble(count_end - count_launch)/dble(count_rate)

  call parallel_end

  ! closing the hdf5 interface
  if(id == master) write(*, *) 'finalizing hdf5'
  call hdf5_finalize()

end program test_main
