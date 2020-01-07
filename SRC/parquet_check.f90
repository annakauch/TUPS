module parquet_check

  use parquet_ini
  use parquet_util
  use hdf5_wrapper
  use hdf5
  use parquet_PhiR    


  implicit none

contains

  !checks convergence as in the old victory code - relative difference of
  !averaged Sigma
  !mainly for comparision with k-space
  logical function conv_old() result(Converged)

    !loop variable
    integer :: idx_k

    relative_err = 0.0d0
    Converged = .false.

    do idx_k = 1, Nt
      relative_err = relative_err + abs(1.0d0 - SigmaOld(idx_k)/Sigma(idx_k))/Nt
    end do

    if(relative_err < f_sigma) Converged = .true.

    if(id == master) then
      write (*, *)
      write (*, '(a, f18.12)') ' Convergence : ---------------- relative error in self-energy is:', relative_err
    end if

  end function conv_old

  !more sophisticated convergence check
  logical function conv() result(Converged)

    !error on average and at maximum
    real(dp) :: rel_sig_av, rel_sig_max
    integer :: k_max(1)
    !for labels
    integer :: kx, ky, nu

    !indx of max_value in F
    integer :: max_pos(3)
    !check convergence in vertices
    real(dp) :: F_d_max_t, F_m_max_t, F_s_max_t, F_t_max_t
    real(dp) :: F_d_max_tot, F_m_max_tot, F_s_max_tot, F_t_max_tot
    real(dp) :: rel_err_d, rel_err_m, rel_err_s, rel_err_t

    ! --- convergence in selfenergy ---

    if(id == master) then
      write(*, *)
      write (*, '(a, f18.12)') ' Convergence :'
    end if


    rel_sig_av = 0.0d0
    rel_sig_max = 0.0d0

    if(id == master) then

      !calculate error via intrinsics
      SigmaOld = abs(SigmaOld - Sigma)/abs(Sigma)
      k_max = maxloc(abs(SigmaOld))
      rel_sig_av = abs(sum(SigmaOld)/size(SigmaOld))
      rel_sig_max = maxval(abs(SigmaOld)) 
  
      !just for convenien labelling
      kx = Index_Fermionic(k_max(1))%ix
      ky = Index_Fermionic(k_max(1))%iy
      nu = Index_Fermionic(k_max(1))%iw

      write (*, *)
      write (*, '(a, f18.12)') ' Sigma averaged error   : ', rel_sig_av
      write (*, '(a, f18.12)') ' Max difference         : ', rel_sig_max
      write (*, '(a, I2.2)', advance='no') ' largest discrepancy at (kx, ky, nu) = (', kx
      write (*, '(a, I2.2)', advance='no') ' ,', ky
      write (*, '(a, I2.2)', advance='no') ' ,', nu
      write (*, *) ')'
    end if

    ! --- convergence in vertices --- only at largest point heuristically 

    !find largest values in F
    max_pos = maxloc(abs(F_d_LL(:, :, :)))
    F_d_max_t = abs(F_d_ll(max_pos(1), max_pos(2), max_pos(3)))

    max_pos = maxloc(abs(F_m_LL(:, :, :)))
    F_m_max_t = abs(F_m_ll(max_pos(1), max_pos(2), max_pos(3)))

    max_pos = maxloc(abs(F_s_LL(:, :, :)))
    F_s_max_t = abs(F_s_ll(max_pos(1), max_pos(2), max_pos(3)))

    max_pos = maxloc(abs(F_t_LL(:, :, :)))
    F_t_max_t = abs(F_t_ll(max_pos(1), max_pos(2), max_pos(3)))

    !calculate relative error for this specific point
    rel_err_d = abs( 1.0d0 - abs(F_d_max_t)/abs(F_d_max) )
    rel_err_m = abs( 1.0d0 - abs(F_m_max_t)/abs(F_m_max) )
    rel_err_s = abs( 1.0d0 - abs(F_s_max_t)/abs(F_s_max) ) 
    rel_err_t = abs( 1.0d0 - abs(F_t_max_t)/abs(F_t_max) )

    !reduce over all MPI processes - to find total maximum
    call MPI_Reduce(rel_err_d, MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_Comm_world, rc)
    call MPI_Reduce(rel_err_m, MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_Comm_world, rc)
    call MPI_Reduce(rel_err_s, MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_Comm_world, rc)
    call MPI_Reduce(rel_err_t, MPI_IN_PLACE, 1, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_Comm_world, rc)

    call MPI_Reduce(F_d_max_t, F_d_max_tot, 1, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_Comm_world, rc)
    call MPI_Reduce(F_m_max_t, F_m_max_tot, 1, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_Comm_world, rc)
    call MPI_Reduce(F_s_max_t, F_s_max_tot, 1, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_Comm_world, rc)
    call MPI_Reduce(F_t_max_t, F_t_max_tot, 1, MPI_DOUBLE_PRECISION, MPI_MAX, master, MPI_Comm_world, rc)

    if(id == master) then

      write(*, *)
  
      write (*, '(a, f18.12)') ' Rel. Err. Max(F_d) : ', rel_err_d 
      write (*, '(a, f18.12)') ' Rel. Err. Max(F_m) : ', rel_err_m 
      write (*, '(a, f18.12)') ' Rel. Err. Max(F_s) : ', rel_err_s 
      write (*, '(a, f18.12)') ' Rel. Err. Max(F_t) : ', rel_err_t 

      write (*, '(a, f18.6)') ' Max(F_d) : ', abs(F_d_max_tot)
      write (*, '(a, f18.6)') ' Max(F_m) : ', abs(F_m_max_tot)
      write (*, '(a, f18.6)') ' Max(F_s) : ', abs(F_s_max_tot)
      write (*, '(a, f18.6)') ' Max(F_t) : ', abs(F_t_max_tot)

    end if 
    
    !set convergence logical
    Converged = (rel_sig_max < f_sigma)
    Converged = Converged .and. (rel_err_d < f_vert) .and. (rel_err_m < f_vert) .and. (rel_err_s < f_vert) .and. (rel_err_t < f_vert)

    !send value of converged to all
    call MPI_BCAST(Converged, 1, MPI_LOGICAL, master, MPI_COMM_WORLD, rc)

    !set saved values for future checks
    SigmaOld = Sigma

    F_d_max = F_d_max_t
    F_m_max = F_m_max_t
    F_s_max = F_s_max_t
    F_t_max = F_t_max_t

  end function conv

  !function to check whether to FT chi0 looks similar to one summed by hand
  !atm this is just for one point
  complex(dp) function check_chi0(box_extension) result(abs_diff)

    integer, intent(in) :: box_extension

    !arguments
    type(Indxmap) :: map_q, map_k, map_kq
    integer :: kx, ky, nu

    !greens functions
    complex(dp) :: Gk, Gkq

    !chi0(pi, pi, 0) summed by hand
    complex(dp) :: chi0_summed

    !calculate chi0(pi,pi,0) by 'hand'
    map_q = Indxmap(Nx/2, Ny/2, 1)
    chi0_summed = 0.0d0

    do kx = 1, Nx
      do ky = 1, Ny
        do nu = - box_extension * Nf + 1, (box_extension + 1) * Nf

          map_k = Indxmap(kx, ky, nu)
          call Index_FaddB(map_k, map_q, map_kq)
          Gk = get_green(map_k)
          Gkq = get_green(map_kq)

          chi0_summed = chi0_summed + Gk * Gkq * 1.0d0 / (Nx * Ny * beta)

        end do !nu
      end do !ky
    end do !kx          

    abs_diff = abs(chi0_summed - chi0_ph(List_Index_IBZ(map_q)))


  end function check_chi0


! --- --- --- --- --- --- ---

  subroutine invTrafo_PhiR()

    !bosonic loop varibales - with and without frequencies included
    integer :: idxQ, idxQexp, q
    !bosonic argument
    type(Indxmap) :: mapQ

    !array to be outputted
    complex(dp), allocatable :: outputBuffer(:)
    !symmetrized version of PhiRd in lattice
    complex(dp), allocatable :: PhiR_sym(:, :, :)

    !for labelling output
    character(len = 5) :: qxstr, qystr, wstr
    character(len = 30) :: label

    !for hdf5 output
    integer(hid_t) :: file_ident

    !loop variables
    integer :: sym, ll, R

    if(.not. allocated(outputBuffer)) allocate(outputBuffer(Nz * Nz))
    if(.not. allocated(PhiR_sym)) allocate(PhiR_sym(Nz * Nz, wpertask, NR))

    do idxQ = 1, Nb

      mapQ = index_bosonic_IBZ(id * Nb + idxQ)
      idxQexp = (list_index_IBZ(mapQ) - mapQ%iw)/(Nf/2) + 1

      !if I have corresponding frequency in PhiR ...
      if(id == getwTask(mapQ%iw)) then

        outputBuffer = 0.0d0

        ! ... perform trasform to qspace
        do sym = 1, Ns
          !act with one symmetry operation on all entries
          call symop_arr(PhiRd, PhiR_sym, sym, mod(mapQ%iw - 1, wperTask) + 1)
          do R = 1, NR
            do ll = 1, Nz * Nz
              outputBuffer(ll) = outputBuffer(ll) + &
                             (1.0d0 / Nx) * PhiR_sym(ll, mod(mapQ%iw - 1, wperTask) + 1, R) * &
                             expqSR((mapQ%iy - 1) * Nx + mapQ%ix, R, sym)
            end do !ll
          end do !R
        end do !sym

        write(qxstr, '(I3.3)') mapQ%ix
        write(qystr, '(I3.3)') mapQ%iy
        write(wstr, '(I3.3)') mapQ%iw
        
        ! produce output for comparison
        call hdf5_open_file(file_main, file_ident)

        label = 'checkPhiR/R_' // trim(qxstr) // '_' // trim(qystr) // '_' // trim(wstr) 
        call hdf5_write_data(file_ident, label, reshape(outputBuffer, (/Nz, Nz/)))

        label = 'checkPhiQ/Q_' // trim(qxstr) // '_' // trim(qystr) // '_' // trim(wstr)
        call hdf5_write_data(file_ident, label, G_d_LL(:, :, idxQ))

        call hdf5_close_file(file_ident)

      end if !id == getwTask(mapQ%iw)
    end do !idxQ

  end subroutine invTrafo_PhiR

  !check if s- and d-wave formfactors were calcualted correctly
  subroutine FF_sdwave_check()

    real(dp) :: check_zero

    check_zero = 0.0d0

    if(readinpick) then 

      write(*, *) 'FF_sdwave_check cannot be performed due to readin of FFs'

    else

      !check whether two ways to calculate this agree  
      check_zero = sum(abs(FF(:, :, 1) - 1.0d0/Nx * FF_sdwave(:, :, 1)))
      if(Nl > 1) then 
        check_zero = check_zero +  sum(abs(FF(:, :, 3) - 1.0d0/Nx * FF_sdwave(:, :, 2)))
        check_zero = check_zero +  sum(abs(FF(:, :, 4) - sqrt(2.0d0)/Nx * FF_sdwave(:, :, 3)))
      end if
  
      if (check_zero < 1e-12) then
        if(id == master) write(*, *) 'Check for s- and d-wave form factor --- SUCCESS!!!'
      else
        if(id == master) write(*, *) 'Check for s- and d-wave form factor --- FAILURE!!!'
        if(id == master) write(*, *) 'Difference: ', check_zero
      end if

    end if !readinpick

  end subroutine FF_sdwave_check

end module parquet_check

