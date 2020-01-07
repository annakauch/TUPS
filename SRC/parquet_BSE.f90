module Parquet_kernel

  ! Purpose
  ! ========
  !   determine the reducible vertex function in particle-hole and particle-particle channels

  use parquet_ini
  use parquet_util
  use parquet_chi
  use mpi_mod
  use hdf5_wrapper
  use hdf5

  implicit none

contains
  !-------------------------------------------------------------------------------------------

  subroutine BSE(ite)

    !iteration count of selfconsistency loop
    integer, intent(in) :: ite

    !bosonic loop variable
    integer :: idx_q
    !bosonic argument corresponding to loop variable
    type(Indxmap) :: map_q

    complex(dp), allocatable :: chi_ph_LL(:, :)
    complex(dp), allocatable :: chi_pp_LL(:, :)

    complex(dp), allocatable :: chiG(:, :)
    complex(dp), allocatable :: Gamtemp(:, :)

    !for output
    character(len = 20) :: label
    integer(hid_t) :: file_ident
    character(len = 5) :: qstr

    integer :: idx1, idx2, idx3


!      if(id == master) then
!        write(qstr, '(I2.2)') idx_q
!        label = 'chi0ph' // trim(qstr)
!        call hdf5_open_file(file_main, file_ident)
!        call hdf5_write_data(file_ident, label, reshape(chi_ph_LL, (/ Nl * Nf, Nl * Nf /)))
!        call hdf5_close_file(file_ident)
!
!      end if !id == master


    if (.not. allocated(chi_ph_LL)) allocate (chi_ph_LL(Nz, Nz))
    if (.not. allocated(chi_pp_LL)) allocate (chi_pp_LL(Nz, Nz))

    if (.not. allocated(chiG)) allocate (chiG(Nz, Nz))
    if (.not. allocated(Gamtemp)) allocate (Gamtemp(Nz, Nz))
    chiG = 0.0d0
    Gamtemp = 0.0d0

    !loop over bosonic argument
    do idx_q = 1, Nb
    
      map_q = Index_Bosonic_IBZ(idx_q + id * Nb)

      !calculate bubbles for specific q
      call calc_chi(map_q, chi_ph_LL, chi_pp_LL)

      ! --- density channel ---

      !first ordering ...
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(1.0d0/(beta * Nx * Ny), 0.0d0), G_d_LL(1:Nz, 1:Nz, idx_q), Nz, chi_ph_LL, Nz, dcmplx(0.0d0, 0.0d0), chiG, Nz)
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(0.5d0 * f_damping, 0.0d0), chiG, Nz, F_d_LL(1:Nz, 1:Nz, idx_q), Nz, dcmplx(0.0d0, 0.0d0), Gamtemp, Nz)


      ! ... and second ordering
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(1.0d0/(beta * Nx * Ny), 0.0d0), chi_ph_LL(1:Nz, 1:Nz), Nz, G_d_LL(1:Nz, 1:Nz, idx_q), Nz, dcmplx(0.0d0, 0.0d0), chiG, Nz)
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(0.5d0 * f_damping, 0.0d0), F_d_LL(1:Nz, 1:Nz, idx_q), Nz, chiG(1:Nz, 1:Nz), Nz, dcmplx(1.0d0, 0.0d0), Gamtemp, Nz)

      !Mixing
      G_d_LL(1:Nz, 1:Nz, idx_q) = Gamtemp + &
                                  (1.0d0 - f_damping) * (F_d_LL(1:Nz,1:Nz, idx_q) - G_d_LL(1:Nz, 1:Nz, idx_q))


      ! --- magnetic channel ---

      !first ordering ...
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(1.0d0/(beta * Nx * Ny), 0.0d0), G_m_LL(1:Nz, 1:Nz, idx_q), Nz, chi_ph_LL(1:Nz, 1:Nz), Nz, dcmplx(0.0d0, 0.0d0), chiG, Nz)
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(0.5d0 * f_damping, 0.0d0), chiG(1:Nz, 1:Nz), Nz, F_m_LL(1:Nz, 1:Nz, idx_q), Nz, dcmplx(0.0d0, 0.0d0), Gamtemp, Nz)

      ! ... and second ordering
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(1.0d0/(beta * Nx * Ny), 0.0d0), chi_ph_LL(1:Nz, 1:Nz), Nz, G_m_LL(1:Nz, 1:Nz, idx_q), Nz, dcmplx(0.0d0, 0.0d0), chiG, Nz)
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(0.5d0 * f_damping, 0.0d0), F_m_LL(1:Nz, 1:Nz, idx_q), Nz, chiG(1:Nz, 1:Nz), Nz, dcmplx(1.0d0, 0.0d0), Gamtemp, Nz)


      !Mixing
      G_m_LL(1:Nz, 1:Nz, idx_q) = Gamtemp + &
                                  (1.0d0 - f_damping) * (F_m_LL(1:Nz,1:Nz, idx_q) - G_m_LL(1:Nz, 1:Nz, idx_q))


      ! --- singlet channel ---

      !first ordering ...
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(1.0d0/(beta * Nx * Ny), 0.0d0), G_s_LL(1:Nz, 1:Nz, idx_q), Nz, chi_pp_LL(1:Nz, 1:Nz), Nz, dcmplx(0.0d0, 0.0d0), chiG, Nz)
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(- 0.25d0 * f_damping, 0.0d0), chiG(1:Nz, 1:Nz), Nz, F_s_LL(1:Nz, 1:Nz, idx_q), Nz, dcmplx(0.0d0, 0.0d0), Gamtemp, Nz)

      ! ... and second ordering
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(1.0d0/(beta * Nx * Ny), 0.0d0), chi_pp_LL(1:Nz, 1:Nz), Nz, G_s_LL(1:Nz, 1:Nz, idx_q), Nz, dcmplx(0.0d0, 0.0d0), chiG, Nz)
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(- 0.25d0 * f_damping, 0.0d0), F_s_LL(1:Nz, 1:Nz, idx_q), Nz, chiG(1:Nz, 1:Nz), Nz, dcmplx(1.0d0, 0.0d0), Gamtemp, Nz)

      !Mixing
      G_s_LL(1:Nz, 1:Nz, idx_q) = Gamtemp + &
                                  (1.0d0 - f_damping) * (F_s_LL(1:Nz,1:Nz, idx_q) - G_s_LL(1:Nz, 1:Nz, idx_q))


      ! --- triplet channel ---

      !first ordering ...
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(1.0d0/(beta * Nx * Ny), 0.0d0), G_t_LL(1:Nz, 1:Nz, idx_q), Nz, chi_pp_LL(1:Nz, 1:Nz), Nz, dcmplx(0.0d0, 0.0d0), chiG, Nz)
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx( 0.25d0 * f_damping, 0.0d0), chiG(1:Nz, 1:Nz), Nz, F_t_LL(1:Nz, 1:Nz, idx_q), Nz, dcmplx(0.0d0, 0.0d0), Gamtemp, Nz)

      ! ... and second ordering
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx(1.0d0/(beta * Nx * Ny), 0.0d0), chi_pp_LL(1:Nz, 1:Nz), Nz, G_t_LL(1:Nz, 1:Nz, idx_q), Nz, dcmplx(0.0d0, 0.0d0), chiG, Nz)
      call ZGEMM('N', 'N', Nz, Nz, Nz, dcmplx( 0.25d0 * f_damping, 0.0d0), F_t_LL(1:Nz, 1:Nz, idx_q), Nz, chiG(1:Nz, 1:Nz), Nz, dcmplx(1.0d0, 0.0d0), Gamtemp, Nz)

      !Mixing
      G_t_LL(1:Nz, 1:Nz, idx_q) = Gamtemp + &
                                  (1.0d0 - f_damping) * (F_t_LL(1:Nz,1:Nz, idx_q) - G_t_LL(1:Nz, 1:Nz, idx_q))



    end do !q


  end subroutine BSE



  subroutine get_kernel_function(ite)
    !
    ! Purpose
    ! =======
    !   From the exact reducible vertex function calculated in routine "reducible_Vertex" to determine the kernel functions
    ! by scanning the edge of it.
    !
    integer, intent(in) :: ite

    ! ... local vars ...
    integer     :: ichannel, i
    integer     :: l1, l2c, q
    complex(dp), allocatable :: dummy1(:, :, :), dummy2(:, :, :)

    type(MPI_Request) :: send_request1, send_request2
    type(MPI_Request) :: recv_request1, recv_request2

    !allocate kernel2 functions for all channels
    if (.NOT. allocated(K2_d1)) allocate (K2_d1(Nz, Nl, Nred))
    if (.NOT. allocated(K2_m1)) allocate (K2_m1(Nz, Nl, Nred))
    if (.NOT. allocated(K2_s1)) allocate (K2_s1(Nz, Nl, Nred))
    if (.NOT. allocated(K2_t1)) allocate (K2_t1(Nz, Nl, Nred))
    if (.NOT. allocated(K2_d2)) allocate (K2_d2(Nl, Nz, Nred))
    if (.NOT. allocated(K2_m2)) allocate (K2_m2(Nl, Nz, Nred))
    if (.NOT. allocated(K2_s2)) allocate (K2_s2(Nl, Nz, Nred))
    if (.NOT. allocated(K2_t2)) allocate (K2_t2(Nl, Nz, Nred))

    ! assign trivial initial values to kernel functions
    if (ite == 1) then
      K2_d1 = 0.0d0
      K2_m1 = 0.0d0
      K2_s1 = 0.0d0
      K2_t1 = 0.0d0
      K2_d2 = 0.0d0
      K2_m2 = 0.0d0
      K2_s2 = 0.0d0
      K2_t2 = 0.0d0
    end if

    ! determine K2 - one frequency argument less
    if (.NOT. allocated(dummy1)) allocate (dummy1(Nz, Nl, Nb))
    if (.NOT. allocated(dummy2)) allocate (dummy2(Nl, Nz, Nb))

    dummy1 = 0.0d0
    dummy2 = 0.0d0

    do ichannel = 1, 4
      do q = 1, Nb
        ! here we will simply scan the edge of the reducible vertex function to get kernel
        do l1 = 1, Nz
          do l2c = 1, Nl !argument which's frequency is out of the box
            select case (ichannel)
            case (1)
              dummy1(l1, l2c, q) = G_d_LL(l1, (l2c - 1) * Nf + Nf, q)
              dummy2(l2c, l1, q) = G_d_LL((l2c - 1) * Nf + Nf, l1, q)
            case (2)
              dummy1(l1, l2c, q) = G_m_LL(l1, (l2c - 1) * Nf + Nf, q)
              dummy2(l2c, l1, q) = G_m_LL((l2c - 1) * Nf + Nf, l1, q)
            case (3)
              dummy1(l1, l2c, q) = G_s_LL(l1, (l2c - 1) * Nf + 1, q)
              dummy2(l2c, l1, q) = G_s_LL((l2c - 1) * Nf + 1, l1, q)
            case (4)
              dummy1(l1, l2c, q) = G_t_LL(l1, (l2c - 1) * Nf + 1, q)
              dummy2(l2c, l1, q) = G_t_LL((l2c - 1) * Nf + 1, l1, q)
            end select
          end do !l1
        end do !l2c

      end do !q

      !now master assembles the complete kernel function
      if (id == master) then
        select case (ichannel)
        case (1)
          K2_d1(1:Nz, 1:Nl, 1:Nb) = dummy1(1:Nz, 1:Nl, 1:Nb)
          K2_d2(1:Nl, 1:Nz, 1:Nb) = dummy2(1:Nl, 1:Nz, 1:Nb)
        case (2)
          K2_m1(1:Nz, 1:Nl, 1:Nb) = dummy1(1:Nz, 1:Nl, 1:Nb)
          K2_m2(1:Nl, 1:Nz, 1:Nb) = dummy2(1:Nl, 1:Nz, 1:Nb)
        case (3)
          K2_s1(1:Nz, 1:Nl, 1:Nb) = dummy1(1:Nz, 1:Nl, 1:Nb)
          K2_s2(1:Nl, 1:Nz, 1:Nb) = dummy2(1:Nl, 1:Nz, 1:Nb)
        case (4)
          K2_t1(1:Nz, 1:Nl, 1:Nb) = dummy1(1:Nz, 1:Nl, 1:Nb)
          K2_t2(1:Nl, 1:Nz, 1:Nb) = dummy2(1:Nl, 1:Nz, 1:Nb)
        end select

        !master recieves other dummies ...
        do i = 1, ntasks - 1

          call MPI_IRECV(dummy1, Nz * Nl * Nb, MPI_DOUBLE_COMPLEX, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_request1, rc)
          call MPI_IRECV(dummy2, Nz * Nl * Nb, MPI_DOUBLE_COMPLEX, i, MPI_ANY_TAG, MPI_COMM_WORLD, recv_request2, rc)
          call MPI_Wait(recv_request1, stat, rc)
          call MPI_Wait(recv_request2, stat, rc)

          ! ... and also writes them into kernel functions

          select case (ichannel)
          case (1)
            K2_d1(1:Nz, 1:Nl, (i * Nb + 1):(i + 1) * Nb) = dummy1(1:Nz, 1:Nl, 1:Nb)
            K2_d2(1:Nl, 1:Nz, (i * Nb + 1):(i + 1) * Nb) = dummy2(1:Nl, 1:Nz, 1:Nb)
          case (2)
            K2_m1(1:Nz, 1:Nl, (i * Nb + 1):(i + 1) * Nb) = dummy1(1:Nz, 1:Nl, 1:Nb)
            K2_m2(1:Nl, 1:Nz, (i * Nb + 1):(i + 1) * Nb) = dummy2(1:Nl, 1:Nz, 1:Nb)
          case (3)
            K2_s1(1:Nz, 1:Nl, (i * Nb + 1):(i + 1) * Nb) = dummy1(1:Nz, 1:Nl, 1:Nb)
            K2_s2(1:Nl, 1:Nz, (i * Nb + 1):(i + 1) * Nb) = dummy2(1:Nl, 1:Nz, 1:Nb)
          case (4)
            K2_t1(1:Nz, 1:Nl, (i * Nb + 1):(i + 1) * Nb) = dummy1(1:Nz, 1:Nl, 1:Nb)
            K2_t2(1:Nl, 1:Nz, (i * Nb + 1):(i + 1) * Nb) = dummy2(1:Nl, 1:Nz, 1:Nb)
          end select
        end do !i

      else !other tasks (not master) send information to master

        call MPI_ISEND(dummy1, Nz * Nl * Nb, MPI_DOUBLE_COMPLEX, master, id, MPI_COMM_WORLD, send_request1, rc)
        call MPI_ISEND(dummy2, Nz * Nl * Nb, MPI_DOUBLE_COMPLEX, master, id, MPI_COMM_WORLD, send_request2, rc)
        call MPI_Wait(send_request1, stat, rc)
        call MPI_Wait(send_request2, stat, rc)

      end if !id == master

    end do !ichannel

    !now send readily assembled kernel functions to all processes
    call MPI_BCAST(K2_d1, Nl * Nz * Nred, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_d2, Nl * Nz * Nred, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_m1, Nl * Nz * Nred, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_m2, Nl * Nz * Nred, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_s1, Nl * Nz * Nred, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_s2, Nl * Nz * Nred, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_t1, Nl * Nz * Nred, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)
    call MPI_BCAST(K2_t2, Nl * Nz * Nred, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, rc)

  end subroutine get_kernel_function

  ! --- function that actually returns the calculated kernel functions ---
  !-------------------------------------------------------------------------------------------------

  complex(dp) pure function kernel(channel, l1, l2, q)

    !for given l1, l2, q = vertex arguments this function return the kernels
    !this function is only to be called if one of the frequency arguments runs out of the box

    !channel for which kernel should be calculated
    character(len=1), intent(in) :: channel
    !input maps in l- ...
    type(Indxmap_L), intent(in) :: l1, l2
    ! ... and k-space
    type(Indxmap), intent(in) :: q

    !out of the box largument
    integer :: l1c, l2c
    !frequency arguments of the input
    integer :: w1, w2, qw

    !in case both frequencies run out of the box
    complex(dp) :: background

    !return value
    kernel = 0.0d0

    !formfactor index without frequency
    l1c = l1%il
    l2c = l2%il

    w1 = l1%iw
    w2 = l2%iw
    qw = q%iw

    !no kernels for running out of q-box
    if (qw > Nf/2 .or. qw < 1) return

    select case (channel)

    case ('d')

      background = K2_d1((l1%il - 1) * Nf + Nf, l2c, List_index_IBZ(q))

      !should only be true in susceptibility
      if (w1 > Nf .or. w1 < 1) then

        if (w2 > Nf .or. w2 < 1) then

          kernel = background

        else

          kernel = K2_d2(l1c, list_index_L(l2), List_Index_IBZ(q))

        end if !w2 out of box

      !only this case needed in PAE and SDE
      else if (w2 > Nf .or. w2 < 1) then

          kernel = K2_d1(list_index_L(l1), l2c, List_Index_IBZ(q))

      else

          !in this case the kernels should not have been called

      end if !w1 out of box

    case ('m')

      background = K2_m1((l1%il - 1) * Nf + Nf, l2c, List_index_IBZ(q))

      !should only be true in susceptibility
      if (w1 > Nf .or. w1 < 1) then

        if (w2 > Nf .or. w2 < 1) then

          kernel = background

        else

          kernel = K2_m2(l1c, list_index_L(l2), List_Index_IBZ(q))

        end if !w2 out of box

      !only this case needed in PAE and SDE
      else if (w2 > Nf .or. w2 < 1) then

          kernel = K2_m1(list_index_L(l1), l2c, List_Index_IBZ(q))

      else

          !in this case the kernels should not have been called

      end if !w1 out of box

    case ('s')

      background = K2_s1((l1%il - 1) * Nf + 1, l2c, List_index_IBZ(q))

      !should only be true in susceptibility
      if (w1 > Nf .or. w1 < 1) then

        if (w2 > Nf .or. w2 < 1) then

          kernel = background

        else

          kernel = K2_s2(l1c, list_index_L(l2), List_Index_IBZ(q))

        end if !w2 out of box

      !only this case needed in PAE and SDE
      else if (w2 > Nf .or. w2 < 1) then

          kernel = K2_s1(list_index_L(l1), l2c, List_Index_IBZ(q))

      else

          !in this case the kernels should not have been called

      end if !w1 out of box

    case ('t')

      !should only be true in susceptibility
      background = K2_t1((l1%il - 1) * Nf + 1, l2c, List_index_IBZ(q))

      if (w1 > Nf .or. w1 < 1) then

        if (w2 > Nf .or. w2 < 1) then

          kernel = background

        else

          kernel = K2_t2(l1c, list_index_L(l2), List_Index_IBZ(q))

        end if !w2 out of box

      !only this case needed in PAE and SDE
      else if (w2 > Nf .or. w2 < 1) then

          kernel = K2_t1(list_index_L(l1), l2c, List_Index_IBZ(q))

      else

          !in this case the kernels should not have been called

      end if !w1 out of box

    end select

  end function kernel



end module Parquet_kernel
