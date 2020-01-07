module parquet_PhiR

  use mpi_mod
  use parquet_ini
  use parquet_util
  use hdf5_wrapper
  use hdf5



  implicit none

  private

  !declare arrays here such that they are not reallocated each time
  !save symmetrization of PhiQ
  complex(dp), allocatable :: PhiQ_sym(:, :, :)

  !temporary memory to save Phi's in q one single node
  complex(dp), allocatable :: PhiQ(:, :, :)

  !buffers - without zgemm semms to segfault - dont know why
  complex(dp), allocatable :: dummy1(:, :), dummy2(:, :)

  public :: calc_PhiR
  public :: getWTask
  public :: task_has_w
  public :: symop_arr
  public :: calculate_expSqR
  public :: calculate_expqSR

contains

!solver for parquet equations in l-space
  subroutine calc_PhiR(ite)
    integer, intent(in) :: ite

    !loop varibales
    integer :: channel, w1


    !for hdf5 output
    integer(hid_t) :: file_ident
    character(len = 30) :: label

    !for timing
    integer(4) :: count_rate, count_max
    !count_launch -> counts at program launch, count_iteration -> count at
    !beginning of iteration
    integer(4) :: count_launch, count_iteration, count_start, count_end

    !allocate memory on nodes where Phi_R is saved
    if(task_has_w()) then

      if(.not. allocated(PhiRd)) allocate(PhiRd(Nz * Nz, wperTask, NR))
      if(.not. allocated(PhiRm)) allocate(PhiRm(Nz * Nz, wperTask, NR))
      if(.not. allocated(PhiRs)) allocate(PhiRs(Nz * Nz, wperTask, NR))
      if(.not. allocated(PhiRt)) allocate(PhiRt(Nz * Nz, wperTask, NR))

      PhiRd = 0.0d0
      PhiRm = 0.0d0
      PhiRs = 0.0d0
      PhiRt = 0.0d0

      if(.not. allocated(PhiQ)) allocate(PhiQ(Nz, Nz, NIBZ * wperTask))
      PhiQ = 0.0d0

    end if !task_has_w

    do channel = 1, 4


      call system_clock(count_start, count_rate, count_max)

      !!get PhiQ for specific Q on one node
      call communicatePhis(channel, PhiQ)

      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time for communicating Phis: ', &!t_end - t_start
                                          dble(count_end - count_start)/dble(count_rate)


      call system_clock(count_start, count_rate, count_max)

      do w1 = 1, Nf/2
        !Only transform stuff on nodes that actually have something to transform
        if(getwTask(w1) == id) then
          !transform to R-space
          call FTPhis(channel, w1, PhiQ)
        end if

      end do !w1

      call system_clock(count_end, count_rate, count_max)
      if (id == master) write (*, "(a, f12.6)") ' time for computing FT: ', &!t_end - t_start
                                          dble(count_end - count_start)/dble(count_rate)

    end do !channel

    !if(ite == ite_max - 1 .and. id == master) then
    !  ! produce output for comparison
    !  call hdf5_open_file(file_main, file_ident)
    !  label = 'checkPhiR/ZPhiRd'
    !  call hdf5_write_data(file_ident, label, reshape(PhiRd, (/Nz, Nz, wperTask * NR/)))
    !  call hdf5_close_file(file_ident)
    !end if


  end subroutine calc_PhiR


  ! ------
  !subroutine to collect Phi(Q) for one frequency on a task
  subroutine communicatePhis(channel, PhiQ)

    !channel for which to do communication
    integer, intent(in) :: channel
    !place to save communicated q-values
    complex(dp), dimension(Nz, Nz, NIBZ * wperTask), intent(out) :: PhiQ

    !buffer to be sent data
    !complex(dp), allocatable :: SendData(:, :, :)
    complex(dp), asynchronous, allocatable :: SendData(:, :, :)

    !loop variables
    integer :: idxQ, Qx, Qy, w1, idxQ_global
    !id of task to send to / receive from
    integer :: target_id, source_id
    type(Indxmap) :: mapQ, mapQ1


    type(MPI_Request), dimension(Nb) :: send_requestQ

    if(.not. allocated(SendData)) allocate(SendData(Nz, Nz, Nb))

    select case(channel)
    case(1)
      SendData = G_d_LL
    case(2)
      SendData = G_m_LL
    case(3)
      SendData = G_s_LL
    case(4)
      SendData = G_t_LL
    end select !channel 

    !Send-Loop
    do idxQ = 1, Nb

      mapQ = index_bosonic_IBZ(id * Nb + idxQ)
      target_id = getwTask(mapQ%iw)
      if(target_id .ne. id) then
        call MPI_ISEND(SendData(:, :, idxQ), Nz * Nz, MPI_DOUBLE_COMPLEX, target_id, id * Nb + idxQ, MPI_COMM_WORLD, send_requestQ(idxQ), rc)
      else
        PhiQ(:, :, indexR_from_indexQ(id * Nb + idxQ)) = SendData(:, :, idxQ)
      end if !frequency not already on my task

    end do !idxQ

    !Receive Data - if one has a correcponding freuquency
    do idxQ_global = 1, NIBZ * Nf/2
      mapQ1 = index_bosonic_IBZ(idxQ_global)
      if(getwTask(mapQ1%iw) .ne. id) cycle
      source_id = get_id(idxQ_global)
      if(source_id .ne. id) then
        call MPI_RECV(PhiQ(:, :, indexR_from_indexQ(idxQ_global)), Nz * Nz, MPI_DOUBLE_COMPLEX, source_id, idxQ_global, MPI_COMM_WORLD, stat)
      end if !source_id .ne. id
    end do !idxQ_global 

    !Loop for waiting for correct sends
    do idxQ = 1, Nb
      mapQ = index_bosonic_IBZ(id * Nb + idxQ)
      target_id = getwTask(mapQ%iw)
      if(target_id .ne. id) then
        call MPI_Wait(send_requestQ(idxQ), stat, rc)
      end if
    end do


    if(allocated(SendData)) deallocate(SendData)

  end subroutine communicatePhis

  !determine on which task corresponding qargument lays
  integer pure function get_id(idxQ) result(id)

    integer, intent(in) :: idxQ

    id = int((idxQ - 1)/Nb)

  end function get_id  


  ! ------
  !calculate R-space representation of Phi for a specific Phi
  subroutine FTPhis(channel, w1, PhiQ)

    !determines channel to transform
    integer, intent(in) :: channel
    !determines bosonic frequency to transform
    integer, intent(in) :: w1
    !place to save communicated q-values
    complex(dp), dimension(Nz, Nz, NIBZ * wperTask), intent(in) :: PhiQ

    !loop variables
    integer :: sym, ll, R, q



    if(.not. allocated(PhiQ_sym)) allocate(PhiQ_sym(Nz * Nz, wperTask, NIBZ))
    if(.not. allocated(dummy1)) allocate(dummy1(Nz * Nz, NIBZ))
    if(.not. allocated(dummy2)) allocate(dummy2(Nz * Nz, NR))
    PhiQ_sym = 0.0d0
    dummy1 = 0.0d0
    dummy2 = 0.0d0

    !now transform Phi_q to get Phi_R
    if(getwTask(w1) == id) then

      do sym = 1, Ns
  
        !act with one symmetry operation on all entries
        call symop_arr(reshape(PhiQ, (/Nz * Nz, wperTask, NIBZ/)), PhiQ_sym, sym, mod(w1 - 1, wperTask) + 1)

        dummy1 = PhiQ_sym(:, mod(w1 - 1, wperTask) + 1, :)
        call ZGEMM('N', 'N', Nz * Nz, NR, NIBZ, dcmplx(1.0d0/Nx, 0.0d0), dummy1(:, :), Nz * Nz, &
        expSqR(:, :, sym), NIBZ, dcmplx(1.0d0, 0.0d0), dummy2(:, :), Nz * Nz)

        select case(channel)

          case(1)
          PhiRd(:, mod(w1 - 1, wpertask) + 1, :) = dummy2
          case(2)
          PhiRm(:, mod(w1 - 1, wpertask) + 1, :) = dummy2
          case(3)
          PhiRs(:, mod(w1 - 1, wpertask) + 1, :) = dummy2
          case(4)
          PhiRt(:, mod(w1 - 1, wpertask) + 1, :) = dummy2

          !case(1)
          !call ZGEMM('N', 'N', Nz * Nz, NR, NIBZ, dcmplx(1.0d0/Nx, 0.0d0), PhiQ_sym(:, mod(w1 - 1, wpertask) + 1, :), Nz * Nz * wperTask, &
          !expSqR(:, :, sym), NIBZ, dcmplx(0.0d0, 0.0d0), PhiRd(:, mod(w1 - 1, wpertask) + 1, :), Nz * Nz * wperTask)
          !case(2)
          !call ZGEMM('N', 'N', Nz * Nz, NR, NIBZ, dcmplx(1.0d0, 0.0d0), PhiQ_sym(:, w1, :), Nz * Nz * wperTask, &
          !expSqR(:, :, sym), NIBZ, dcmplx(1.0d0, 0.0d0), PhiRm(:, w1, :), Nz * Nz * wperTask)
          !case(3)
          !call ZGEMM('N', 'N', Nz * Nz, NR, NIBZ, dcmplx(1.0d0, 0.0d0), PhiQ_sym(:, w1, :), Nz * Nz * wperTask, &
          !expSqR(:, :, sym), NIBZ, dcmplx(1.0d0, 0.0d0), PhiRs(:, w1, :), Nz * Nz * wperTask)
          !case(4)
          !call ZGEMM('N', 'N', Nz * Nz, NR, NIBZ, dcmplx(1.0d0, 0.0d0), PhiQ_sym(:, w1, :), Nz * Nz * wperTask, &
          !expSqR(:, :, sym), NIBZ, dcmplx(1.0d0, 0.0d0), PhiRt(:, w1, :), Nz * Nz * wperTask)
  
        end select

      end do !sym
    end if!getwTask == id

  end subroutine FTPhis



  ! ------ some utitily stuff ------


  !determine if I (current task) save a frequency argument for PhiR
  logical pure function task_has_w result(has_w)

    integer :: w

    has_w = .false.
    do  w = 1, Nf/2
      if(getwTask(w) == id) has_w = .true.
    end do !w
  end function task_has_w

  !for given frequency return the task that it recides on
  integer pure function getwTask(w) result(ident)
    
    integer, intent(in) :: w
    integer :: step

    !if one has few processes just fill up with appropriate number of frequencies
    if(ntasks .le. Nf/2) then
       ident = int((w-1)/wperTask)
       return
    end if !ntasks < Nf/2
    
    !if one has many task distribute memory over different nodes as good as possible
    step = ntasks / (Nf/2)
    ident = (w - 1) * step

  end function getwTask 

  !function to get an index for PhiR from a usual q-index in IBZ
  integer pure function indexR_from_indexQ(idxQ) result (idxR)

    !idxQ is a global index - not local on a particular node
    integer, intent(in) :: idxQ
    type(indxmap) :: mapQ

    mapQ = index_bosonic_IBZ(idxQ)
    !now idxR hold only IRlat-part - 1 ...
    idxR = (idxQ - mapQ%iw) / (Nf/2)
    !... again add frequencies in case more than one frequency recides on each task
    idxR = idxR * wperTask + mod(mapQ%iw - 1, wperTask) + 1

  end function indexR_from_indexQ

  
  !act with a symmetry operation an a whole array
  subroutine symop_arr(tosym, retsym, sym, w_outer)

    !last argument are momenta in IBZ
    !array that is supposed to be symmetrized
    complex(dp), dimension(Nz * Nz, wpertask, NIBZ), intent(in) :: tosym
    !readily symmterized return array
    complex(dp), dimension(Nz * Nz, wpertask, NIBZ), intent(out) :: retsym
    !symmetry operation to be performed
    integer, intent(in) :: sym
    !bosonic frequency argument on current node
    integer, intent(in) :: w_outer

    ! loop variable ...
    integer :: idxQ

    !prefactor after symmetry operation
    real(dp) :: prefacl1, prefacl2
    !new index after symmetry operation
    integer :: L1, L2, L1_0, L2_0, L1sym, L2sym
    integer :: Lidx_sym, Lidx_nosym

    retsym = 0.0d0

    do idxQ = 1, NIBZ
      do L2 = 1, Nl
        do L1 = 1, Nl
          call sym_op_L(sym, L2, L2sym, prefacL2)
          call sym_op_L(sym, L1, L1sym, prefacL1)
          
          do L2_0 = 1, Nf
            do L1_0 = 1, Nf
  
              Lidx_sym = ((L2sym - 1) * Nf + L2_0 - 1) * Nz + (L1sym - 1) * Nf + L1_0
              Lidx_nosym = ((L2 - 1) * Nf + L2_0 - 1) * Nz + (L1 - 1) * Nf + L1_0

              retsym(Lidx_nosym, w_outer, idxQ) = prefacL1 * prefacL2 * &
                                           tosym(Lidx_sym, w_outer, idxQ)


            end do !L1_0
          end do !L2_0
        end do !L1
      end do !L2
    end do !idxQ

  end subroutine symop_arr

  !save trafo matrix for FT
  !behaves differently to below function - contains zeros
  subroutine calculate_expSqR()

    !loop variables and indices
    integer :: idxQ, idxR, sym
    !actual arguments corresponding to indices
    type(Indxmap) :: mapQ, mapQ_sym, mapR
    !intermediate variable to save exponent 
    real(dp) :: QdotR

    !save whether q-point needs to be symmetrized
    logical, dimension(8) :: symList

    if(.not. allocated(expSqR)) allocate(expSqR(NIBZ, NR, Ns))
    expSqR = 0.0d0  

    do idxQ = 1, NIBZ

      mapQ = index_bosonic_IBZ((idxQ - 1) * Nf/2 + 1)
      !avoid counting points twice - leave their value at zero
      call list_symmetries(mapQ, symList) 

      do sym = 1, Ns
        if(.not. symList(sym)) cycle
        call symmetry_operation(sym, mapQ, mapQ_sym)
        do idxR = 1, NR !calculate full trafo matrix
          !just treat real space index like an IBZ index
          mapR = index_bosonic_IBZ((idxR - 1) * Nf/2 + 1)

          QdotR = (2.0d0 * pi / Nx) * (mapQ_sym%ix - 1) * (mapR%ix - 1) + (2.0d0 * pi / Ny) * (mapQ_sym%iy - 1) * (mapR%iy - 1)
          expSqR(idxQ, idxR, sym) = exp(dcmplx(0.0d0, QdotR))

        end do !idxR
      end do !sym
    end do !idxQ

  end subroutine calculate_expSqR


  !save trafo matrix for parquet equation
  !actually behaves differently than above function - no zeros contained
  subroutine calculate_expqSR()

    !loop variables and indices
    integer :: kx, ky, idxK, idxR, sym
    !actual arguments corresponding to indices
    type(Indxmap) :: mapK, mapR_sym, mapR
    !intermediate variable to save exponent 
    real(dp) :: KdotSR

    !save whether q-point needs to be symmetrized
    logical, dimension(8) :: symList

    if(.not. allocated(expqSR)) allocate(expqSR(Nx * Ny, NR, Ns))
    expqSR = 0.0d0  

    do idxR = 1, NR !calculate full trafo matrix
      !just treat real space index like an IBZ index
      mapR = index_bosonic_IBZ((idxR - 1) * Nf/2 + 1)
      call list_symmetries(mapR, symList) 
      do sym = 1, Ns
        if(.not. symList(sym)) cycle
        call symmetry_operation(sym, mapR, mapR_sym)
        do ky = 1, Ny
        do kx = 1, Nx
          mapK = Indxmap(kx, ky, 1)
          idxK = (ky - 1) * Nx + kx

          KdotSR = (2.0d0 * pi / Nx) * (mapK%ix - 1) * (mapR_sym%ix - 1) + (2.0d0 * pi / Ny) * (mapK%iy - 1) * (mapR_sym%iy - 1)
          expqSR(idxK, idxR, sym) = exp(dcmplx(0.0d0, - KdotSR))

        end do !ky
        end do !kx
      end do !sym
    end do !idxR

  end subroutine calculate_expqSR


end module parquet_PhiR

