module preprocessing

  use parquet_check
  use global_parameter
  use mpi_mod
  use parquet_ini
  use parquet_chi
  use parquet_formfactors
  use parquet_util
  use parquet_plot
  use hdf5_wrapper
  use hdf5
  

  implicit none

contains

  subroutine do_preprocessing
  
  
    implicit none
  
    !for initialization loop
    integer :: idx_q, q_0
    !1/ff(0)/ff(0)
    real(dp) :: aux
  
  
    !for timing
    integer(4) :: count_rate, count_max
    !beginning of iteration
    integer(4) :: count_start, count_end
  
    !for plotting
    character(len=100) :: label, kind_str_temp, sus_str
  
    !format to write integers into strings + parameter names
    character(len=8) :: fmt, FF_num, mom_num, Nf_str, B_str, U_str, grain_str, nParticle_str
  
    !to label plots
    integer :: beta_int, U_int
  
    !write timing message to screen after each iteration
    character(len=16)  :: it_str, time_str
    character(len=50) :: iteration_message
  
    !specify format for label string
    !integer with 3 digits
    fmt = '(I4.4)'
  
    !set basic paramerters
    call set_parameters
  
    !readin data + some allocations
    call Initialization
  
    !calculate array with formfactors
    call FF_array
    !calculate s- and d-wave formfactor seperately
    call fill_FF_sdwave
    !check that s- and d-wave are calculated correctly
    call FF_sdwave_check 
  
    !some strings for naming
    write(FF_num, fmt) Nl
    write(mom_num, fmt) Nx * Ny
    write(Nf_str, fmt) Nf
    beta_int = INT(beta)
    write(B_str, '(I2.2)') beta_int
    U_int = INT(10 * xU)
    write(U_str, fmt) U_int
    write(grain_str, fmt) Ngrain
    write(nParticle_str, fmt) INT(nParticle * 100)
  
    kind_str = 'NxN' // trim(mom_num) // '_Nl' // trim(FF_num) // '_Nf' &
                // trim(Nf_str) // '_U' // trim(U_str) // '_B' // trim(B_str) &
                // '_Gr' // trim(grain_str) // '_n' // trim(nParticle_str)
  
    !to distiguish PA from DGA calculations
    if(.not. dga) kind_str = trim(kind_str) // "_pa"
    if(use_U2) kind_str = trim(kind_str) // "_U2"
    if(oneshot) kind_str = trim(kind_str) // "_one"
  
    if(old) kind_str = trim(kind_str) // "_old"
    if(.not. old) kind_str = trim(kind_str) // "_new"
    if(SDE_old) kind_str = trim(kind_str) // "_SDEold"
    if(.not. SDE_old) kind_str = trim(kind_str) // "_SDEnew"

    !kind_str = trim(kind_str) // "2"

    if(id == master) write(*, *) "label on output: ", kind_str
  
    !check that Nx*Ny*Nf/2 is a multiple of ntasks
    if( ntasks * (int(Nred/ntasks)) .NE. Nred ) then
      write(*, *) 'ERROR - choose number of tasks such that it divides Nred'
      stop
    end if
  
    aux = abs(1.0d0 / (FF(1, 1, 1) * FF(1, 1, 1)))
  
    do idx_q = 1, Nb
  
      q_0 = index_bosonic_IBZ(id * Nb + idx_q)%iw
  
      F_d_LL(1:Nf, 1:Nf, idx_q) = L_d(1:Nf, 1:Nf, q_0) * aux 
      F_m_LL(1:Nf, 1:Nf, idx_q) = L_m(1:Nf, 1:Nf, q_0) * aux
      F_s_LL(1:Nf, 1:Nf, idx_q) = L_s(1:Nf, 1:Nf, q_0) * aux
      F_t_LL(1:Nf, 1:Nf, idx_q) = L_t(1:Nf, 1:Nf, q_0) * aux
  
      G_d_LL(1:Nf, 1:Nf, idx_q) = L_d(1:Nf, 1:Nf, q_0) * aux 
      G_m_LL(1:Nf, 1:Nf, idx_q) = L_m(1:Nf, 1:Nf, q_0) * aux
      G_s_LL(1:Nf, 1:Nf, idx_q) = L_s(1:Nf, 1:Nf, q_0) * aux
      G_t_LL(1:Nf, 1:Nf, idx_q) = L_t(1:Nf, 1:Nf, q_0) * aux
  
    end do
  
    !allocate real time/space green's function
    if(.not. allocated(Grt)) allocate(Grt(Nx * Ngrain, Ny * Ngrain, Nf))
  
    !call cpu_time(t_end)
    call system_clock(count_end, count_rate, count_max)
    if (id == master) write (*, "(a, f12.6)") ' time spent on initialization process: ', &
                                              dble(count_end - count_start)/dble(count_rate)
  
    !set NR to something appropriate
    if(.not. old) then
      if(Nl == 1) NR = 1
      if(Nl == 5) NR = 4
      if(Nl == 9) NR = 6
    else
      NR = 1
    end if
    !NR = NIBZ
  
    !write parameters 
    if (id == master) then
      write (*, *) 'Parameters : '
      write (*, *) ' Beta : ', beta
      write (*, *) ' U : ', xU
      write (*, *) ' Nx : ', Nx
      write (*, *) ' Nf : ', Nf
      write (*, *) ' Nl : ', Nl
      write (*, *) ' NR : ', NR
      write (*, *) ' ite_max : ', ite_max
      write (*, *) 'DGA : ', DGA
  
      write (*, *) 
      write (*, *) 'fixed : '
      write (*, *) 'Ngrain : ', Ngrain
      write (*, *) 'f_range : ', f_range 
      write (*, *) 'f_damping : ', f_damping
      write (*, *) 'nParticle : ', nParticle
      write (*, *) 'f_sigma : ', f_sigma
      write (*, *) 'f_vert : ', f_vert
      write (*, *) 'f_kernel : ', f_kernel
    end if
  
  
    !setup hdf5
    file_main = 'data/TUPS_' // trim(kind_str) // '.hdf5'
    if(id == master) write(*, *) 'writing program output to ', file_main
    if(id == master) call setup_hdf5()
  
  end subroutine do_preprocessing

end module preprocessing
