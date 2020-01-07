module MPI_mod
  use mpi_f08
  !use mpi
  use global_parameter

implicit none

contains
  !-------------------------------------------------------------------------------------------
  subroutine parallel_start

    call MPI_INIT(rc)
    if (rc /= MPI_SUCCESS) then
       print*, "MPI initialization failed"
       stop
    end if
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, rc)
    call MPI_COMM_RANK(MPI_COMM_WORLD, id, rc)

    !allocate(stat(MPI_STATUS_SIZE))
    !allocate(send_request(0:ntasks-1))
    !allocate(recv_request(0:ntasks-1))

  end subroutine parallel_start

  !-------------------------------------------------------------------------------------------
  subroutine parallel_end

    !if (allocated(stat))  deallocate(stat)
    !if (allocated(send_request)) deallocate(send_request)
    !if (allocated(recv_request)) deallocate(recv_request)

    call MPI_FINALIZE(rc)

  end subroutine parallel_end

end module MPI_mod
