module parquet_EQOuter

  use mpi_mod
  use parquet_ini
  use parquet_util
  use parquet_kernel
  use parquet_PhiR
  use parquet_EQContribs


  implicit none

contains


  subroutine calculate_PAE

    !loop variable to add fully irreducible vertex
    integer :: idxQ
    !loop variables for parquet loop
    integer :: channel, task, w2

    !frequency corresponding to idxQ
    integer :: q0

    !normalization of local quantities
    real(dp) :: aux
  
    !temporary matrix to be sent around and worked with
    complex(dp), allocatable :: PhiRWork(:, :, :, :)

    !for timing
    integer(4) :: count_rate, count_max
    !count_launch -> counts at program launch, count_iteration -> count at
    !beginning of iteration
    integer(4) :: count_launch, count_iteration, count_start, count_end


    if(.not. allocated(PhiRWork)) allocate(PhiRWork(Nz, Nz, wperTask, NR))


    ! --- reducible vertex + fully irreducible vertex ---
    F_d_LL = G_d_LL
    F_m_LL = G_m_LL
    F_s_LL = G_s_LL
    F_t_LL = G_t_LL

    aux = abs( 1.0d0/FF(1, 1, 1) * 1.0d0/FF(1, 1, 1) )

    do idxQ = 1, Nb

      q0 = index_bosonic_IBZ(id * Nb + idxQ)%iw

      F_d_LL(1:Nf, 1:Nf, idxQ) = F_d_LL(1:Nf, 1:Nf, idxQ) + L_d(1:Nf, 1:Nf, q0) * aux 
      F_m_LL(1:Nf, 1:Nf, idxQ) = F_m_LL(1:Nf, 1:Nf, idxQ) + L_m(1:Nf, 1:Nf, q0) * aux
      F_s_LL(1:Nf, 1:Nf, idxQ) = F_s_LL(1:Nf, 1:Nf, idxQ) + L_s(1:Nf, 1:Nf, q0) * aux
      F_t_LL(1:Nf, 1:Nf, idxQ) = F_t_LL(1:Nf, 1:Nf, idxQ) + L_t(1:Nf, 1:Nf, q0) * aux

    end do !idxQ

    !for kernel treatment
    call fillpickNu()

    do channel = 1, 4
      do w2 = 1, Nf/2

        if(id == getwTask(w2)) then
          select case(channel)
            case(1)
            PhiRWork = reshape(PhiRd, (/Nz, Nz, wperTask, NR/))
            case(2)
            PhiRWork = reshape(PhiRm, (/Nz, Nz, wperTask, NR/))
            case(3)
            PhiRWork = reshape(PhiRs, (/Nz, Nz, wperTask, NR/))
            case(4)
            PhiRWork = reshape(PhiRt, (/Nz, Nz, wperTask, NR/))
          end select 
        end if !id == getwTask(w2)

        !send working-data to all tasks
        call MPI_BCAST(PhiRWork, Nz * Nz * wperTask * NR, MPI_DOUBLE_COMPLEX, getwTask(w2), MPI_COMM_WORLD, rc)

        !actual evaluation of projected parquet equations
        select case(channel)
        case(1)
          call calc_phbar_to_ph(PhiRWork, w2, -0.5d0, -0.5d0)
          call calc_ph_to_pp(PhiRWork, w2, 0.5d0, 0.5d0)
          call calc_phbar_to_pp(PhiRWork, w2, 0.5d0, -0.5d0)
          if(w2 > 1) then
            call calc_phbar_to_ph_TR(PhiRWork, w2, -0.5d0, -0.5d0)
            call calc_ph_to_pp_TR(PhiRWork, w2, 0.5d0, 0.5d0)
            call calc_phbar_to_pp_TR(PhiRWork, w2, 0.5d0, -0.5d0)
          end if !w1 > 1
        case(2)
          call calc_phbar_to_ph(PhiRWork, w2, -1.5d0, 0.5d0)
          call calc_ph_to_pp(PhiRWork, w2, -1.5d0, 0.5d0)
          call calc_phbar_to_pp(PhiRWork, w2, -1.5d0, -0.5d0)
          if(w2 > 1) then
            call calc_phbar_to_ph_TR(PhiRWork, w2, -1.5d0, 0.5d0)
            call calc_ph_to_pp_TR(PhiRWork, w2, -1.5d0, 0.5d0)
            call calc_phbar_to_pp_TR(PhiRWork, w2, -1.5d0, -0.5d0)
          end if !w1 > 1
        case(3)
          call calc_pp_to_ph(PhiRWork, w2, 0.5d0, -0.5d0)
          if(w2 > 1) then
            call calc_pp_to_ph_TR(PhiRWork, w2, 0.5d0, -0.5d0)
          end if
        case(4)
          call calc_pp_to_ph(PhiRWork, w2, 1.5d0, 0.5d0)
          if(w2 > 1) then
            call calc_pp_to_ph_TR(PhiRWork, w2, 1.5d0, 0.5d0)
          end if
        end select

      end do !w2
    end do !channel

  end subroutine calculate_PAE


end module parquet_EQOuter

