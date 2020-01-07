module parquet_EQMatrix

  use mpi_mod
  use parquet_ini
  use parquet_util
  use parquet_kernel


  implicit none

  !tranformation tensors
  complex(dp), allocatable :: M_phbar_to_ph(:, :, :, :, :, :)
  complex(dp), allocatable :: M_phbar_to_ph_long(:, :, :, :, :, :)
  complex(dp), allocatable :: M_pp_to_ph(:, :, :, :, :, :)
  complex(dp), allocatable :: M_ph_to_pp(:, :, :, :, :, :)
  complex(dp), allocatable :: M_phbar_to_pp(:, :, :, :, :, :)


contains

  subroutine calc_M_phbar_to_ph()

    integer :: sym, kx, ky, idxQ
    integer :: l1, l2, l3, l4
    integer :: R2
    type(Indxmap) :: mapQ

    !upper bound to number of 'bosoic momenta residing on each task'
    integer :: numQ

    !intermediate results - to be summed over symmetries
    complex(dp), allocatable :: M1(:, :, :, :, :, :)
    complex(dp), allocatable :: M2(:, :, :)

    numQ = int((Nb + Nf/2 - 1)/(Nf/2)) + 1

    if(.not. allocated(M1)) allocate(M1(Nl, Nl, Nl, Nb, NR, Ns))
    if(.not. allocated(M2)) allocate(M2(Nl, NR, Ns))
    if(.not. allocated(M_phbar_to_ph)) allocate(M_phbar_to_ph(Nl, Nl, Nl, Nl, Nb, NR))
    M1 = 0.0d0
    M2 = 0.0d0
    M_phbar_to_ph = 0.0d0
    
    !calculate M2
    do sym = 1, Ns
      do R2 = 1, NR
        do l2 = 1, Nl
          !sum over k
          do ky = 1, Ny
          do kx = 1, Nx

          M2(l2, R2, sym) = M2(l2, R2, sym) + &
                            expqSR((ky - 1) * Nx + kx, R2, sym) * &
                            FF(kx, ky, l2)

          end do !kx
          end do !ky

        end do !l2
      end do !R2
    end do !sym

    !calculate M1
    do sym = 1, Ns
      do R2 = 1, NR
        do idxQ = 1, Nb
          mapQ = index_bosonic_IBZ(id * Nb + idxQ)
          do l4 = 1, Nl
          do l3 = 1, Nl
          do l1 = 1, Nl
            !sum over k
            do ky = 1, Ny
            do kx = 1, Nx

            M1(l1, l3, l4, idxQ, R2, sym) = M1(l1, l3, l4, idxQ, R2, sym) + &
                              conjg(expqSR((ky - 1) * Nx + kx, R2, sym)) * &
                              FF(kx, ky, l1) * &
                              FF_inv(sym, kx, ky, l3) * &
                              FF_inv(sym, mod(kx + mapQ%ix - 2, Nx) + 1, mod(ky + mapQ%iy - 2, Ny) + 1, l4)

            end do !kx
            end do !ky

          end do !l1
          end do !l3
          end do !l4
        end do !idxQ
      end do !R2
    end do !sym
 
    !now calculate transformaton matrix
    do R2 = 1, NR
      do idxQ = 1, Nb
        do l4 = 1, Nl
        do l3 = 1, Nl
        do l2 = 1, Nl
        do l1 = 1, Nl
          !summing symmetries
          do sym = 1, Ns
    
            M_phbar_to_ph(l1, l2, l3, l4, idxQ, R2) = M_phbar_to_ph(l1, l2, l3, l4, idxQ, R2) + &
                                                      M1(l1, l3, l4, idxQ, R2, sym) * M2(l2, R2, sym)

          end do !sym
        end do !l1
        end do !l2
        end do !l3
        end do !l4
      end do !idxQ  
    end do !R2

    !write(*, *)
    !write(*, *)
    !!check that correct number of matrix entries is zero
    !do R2 = 1, NR
    !  
    !  write(*, *) 'R2 = ',R2
    !  write(*, *) 'M_phbar_to_ph(..., R2) = ', M_phbar_to_ph(1, 1, 1, 1, 3, R2)

    !end do !R2

 
  end subroutine calc_M_phbar_to_ph


  subroutine calc_M_pp_to_ph()

    integer :: sym, kx, ky, idxQ
    integer :: l1, l2, l3, l4
    integer :: R2
    type(Indxmap) :: mapQ

    !intermediate results - to be summed over symmetries
    complex(dp), allocatable :: M1(:, :, :, :)

    if(.not. allocated(M1)) allocate(M1(Nl, Nl, NR, Ns))
    if(.not. allocated(M_pp_to_ph)) allocate(M_pp_to_ph(Nl, Nl, Nl, Nl, Nb, NR))
    M1 = 0.0d0
    M_pp_to_ph = 0.0d0

    do sym = 1, Ns
      do R2 = 1, NR
        do l3 = 1, Nl   
        do l1 = 1, Nl   

          do ky = 1, Ny
          do kx = 1, Nx

          M1(l1, l3, R2, sym) = M1(l1, l3, R2, sym) + & 
                                FF(kx, ky, l1) * FF_inv(sym, kx, ky, l3) * &
                                expqSR((ky - 1) * Nx + kx, R2, sym)
          end do !kx
          end do !ky

        end do !l1
        end do !l3
      end do !R2
    end do !sym

    do R2 = 1, NR
      do idxQ = 1, Nb
        mapQ = index_bosonic_IBZ(id * Nb + idxQ)
        do l4 = 1, Nl
        do l3 = 1, Nl
        do l2 = 1, Nl
        do l1 = 1, Nl

          do sym = 1, Ns

          M_pp_to_ph(l1, l2, l3, l4, idxQ, R2) = & 
              M_pp_to_ph(l1, l2, l3, l4, idxQ, R2) + &
              expqSR((mapQ%iy - 1) * Nx + mapQ%ix, R2, sym) * &
              M1(l1, l3, R2, sym) * M1(l2, l4, R2, sym)

          end do !sym

        end do !l1
        end do !l2
        end do !l3
        end do !l4
      end do !idxQ
    end do !R2

    !write(*, *)
    !write(*, *)
    !!check that correct number of matrix entries is zero
    !do R2 = 1, NR
    !  
    !  write(*, *) 'R2 = ',R2
    !  write(*, *) 'M_pp_to_ph(..., R2) = ', M_pp_to_ph(1, 1, 1, 1, 3, R2)

    !end do !R2

  end subroutine calc_M_pp_to_ph



  subroutine calc_M_ph_to_pp()

    integer :: sym, kx, ky, idxQ
    integer :: l1, l2, l3, l4
    integer :: R2
    type(Indxmap) :: mapQ

    !upper bound to number of 'bosoic momenta residing on each task'
    integer :: numQ

    !intermediate results - to be summed over symmetries
    complex(dp), allocatable :: M1(:, :, :, :)

    numQ = int((Nb + Nf/2 - 1)/(Nf/2)) + 1

    if(.not. allocated(M1)) allocate(M1(Nl, Nl, NR, Ns))
    if(.not. allocated(M_ph_to_pp)) allocate(M_ph_to_pp(Nl, Nl, Nl, Nl, Nb, NR))
    M1 = 0.0d0
    M_ph_to_pp = 0.0d0

    do sym = 1, Ns
      do R2 = 1, NR
        do l3 = 1, Nl   
        do l1 = 1, Nl   

          do ky = 1, Ny
          do kx = 1, Nx

          M1(l1, l3, R2, sym) = M1(l1, l3, R2, sym) + & 
                                FF(kx, ky, l1) * FF_inv(sym, kx, ky, l3) * &
                                conjg(expqSR((ky - 1) * Nx + kx, R2, sym))
          end do !kx
          end do !ky

        end do !l1
        end do !l3
      end do !R2
    end do !sym

    do R2 = 1, NR
      do idxQ = 1, Nb
        mapQ = index_bosonic_IBZ(id * Nb + idxQ)
        do l4 = 1, Nl
        do l3 = 1, Nl
        do l2 = 1, Nl
        do l1 = 1, Nl

          do sym = 1, Ns

          M_ph_to_pp(l1, l2, l3, l4, idxQ, R2) = & 
              M_ph_to_pp(l1, l2, l3, l4, idxQ, R2) + &
              expqSR((mapQ%iy - 1) * Nx + mapQ%ix, R2, sym) * &
              M1(l1, l3, R2, sym) * M1(l2, l4, R2, sym)

          end do !sym

        end do !l1
        end do !l2
        end do !l3
        end do !l4
      end do !idxQ
    end do !R2

    !write(*, *)
    !write(*, *)
    !!check that correct number of matrix entries is zero
    !do R2 = 1, NR
    !  
    !  write(*, *) 'R2 = ',R2
    !  write(*, *) 'M_ph_to_pp(..., R2) = ', M_ph_to_pp(1, 1, 1, 1, 3, R2)

    !end do !R2

  end subroutine calc_M_ph_to_pp


  subroutine calc_M_phbar_to_pp()

    integer :: sym, kx, ky, idxQ
    integer :: l1, l2, l3, l4
    integer :: R2
    type(Indxmap) :: mapQ

    !intermediate results - to be summed over symmetries
    complex(dp), allocatable :: M1(:, :, :, :), M2(:, :, :, :, :)

    if(.not. allocated(M1)) allocate(M1(Nl, Nl, NR, Ns))
    if(.not. allocated(M2)) allocate(M2(Nl, Nl, Nb, NR, Ns))
    if(.not. allocated(M_phbar_to_pp)) allocate(M_phbar_to_pp(Nl, Nl, Nl, Nl, Nb, NR))
    M1 = 0.0d0
    M2 = 0.0d0
    M_phbar_to_pp = 0.0d0

    !calculate M1
    do sym = 1, Ns
      do R2 = 1, NR
        do l3 = 1, Nl   
        do l1 = 1, Nl   

          do ky = 1, Ny
          do kx = 1, Nx

          M1(l1, l3, R2, sym) = M1(l1, l3, R2, sym) + & 
                                FF(kx, ky, l1) * FF_inv(sym, kx, ky, l3) * &
                                conjg(expqSR((ky - 1) * Nx + kx, R2, sym))
          end do !kx
          end do !ky

        end do !l1
        end do !l3
      end do !R2
    end do !sym

    !calculate M2
    do sym = 1, Ns
      do R2 = 1, NR
        do idxQ = 1, Nb
          mapQ = index_bosonic_IBZ(id * Nb + idxQ)
          do l4 = 1, Nl   
          do l2 = 1, Nl   
  
            do ky = 1, Ny
            do kx = 1, Nx
  
            M2(l2, l4, idxQ, R2, sym) = M2(l2, l4, idxQ, R2, sym) + & 
                                  FF_inv(sym, mod(mapQ%ix - kx + Nx, Nx) + 1, mod(mapQ%iy - ky + Ny, Ny) + 1, l4) * &
                                  FF(kx, ky, l2) * &
                                  expqSR((ky - 1) * Nx + kx, R2, sym)
            end do !kx
            end do !ky
  
          end do !l1
          end do !l3
        end do !idxQ
      end do !R2
    end do !sym


    do R2 = 1, NR
      do idxQ = 1, Nb
        mapQ = index_bosonic_IBZ(id * Nb + idxQ)
        do l4 = 1, Nl
        do l3 = 1, Nl
        do l2 = 1, Nl
        do l1 = 1, Nl

          do sym = 1, Ns

          M_phbar_to_pp(l1, l2, l3, l4, idxQ, R2) = & 
              M_phbar_to_pp(l1, l2, l3, l4, idxQ, R2) + &
              M1(l1, l3, R2, sym) * M2(l2, l4, idxQ, R2, sym)

          end do !sym

        end do !l1
        end do !l2
        end do !l3
        end do !l4
      end do !idxQ
    end do !R2


  end subroutine calc_M_phbar_to_pp



end module parquet_EQMatrix

