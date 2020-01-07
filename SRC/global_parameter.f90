module global_parameter

  use mpi_f08  
  !use mpi

  implicit none

  !specifies precision of floats
  integer,  parameter :: dp = 8
  real(dp), parameter :: pi = acos(-1.d0)

  !numbers that still lurk somewhere in the code
  real(dp), parameter :: Zero = 0.d0, Half = 0.5d0, One = 1.d0, Two = 2.d0  
  complex(dp), parameter :: xi = dcmplx(Zero, One)
  complex(dp), parameter :: Zero_c = dcmplx(Zero, Zero), One_c = dcmplx(One, Zero)
  complex(dp), parameter :: Half_c = dcmplx(Half, Zero), Two_c  = dcmplx(Two, Zero)

  !--- mpi ---
  !rc - error message
  integer  :: ntasks, master, id, rc
  !integer, allocatable :: stat(:)
  !integer, allocatable :: send_request(:), recv_request(:)
  !mpi_f08
  type(MPI_Request) :: send_request
  type(MPI_Request) :: recv_request
  type(MPI_Status) :: stat

  integer, parameter :: nSpin   = 2

  integer, parameter :: nTau    = 200


end module

