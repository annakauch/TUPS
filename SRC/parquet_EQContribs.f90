module parquet_EQContribs

  use mpi_mod
  use parquet_ini
  use parquet_util
  use parquet_kernel
  use parquet_EQMatrix

  implicit none

  integer, allocatable :: pickNu(:)
  integer, allocatable :: pickNu_kernel(:)


contains

  subroutine fillpickNu()

    integer :: nu

    !take Nf + 2*Nf on both sides - should be enough to cover box
    if(.not. allocated(pickNu)) allocate(pickNu(5 * Nf))
    if(.not. allocated(pickNu_kernel)) allocate(pickNu_kernel(5 * Nf))

    do nu = 1, 5 * Nf
  
      if(nu .le. 2 * Nf) then
        pickNu(nu) = 1
        pickNu_kernel(nu) = Nf
      else if(2 * Nf < nu .and. nu < 3 * Nf) then
        pickNu(nu) = nu - 2 * Nf
        pickNu_kernel(nu) = nu - 2 * Nf
      else 
        pickNu(nu) = Nf
        pickNu_kernel(nu) = Nf
      end if !nu

    end do !nu
  end subroutine fillpickNu




  subroutine calc_phbar_to_ph(PhiRWork, w2, fac1, fac2)

    !R-space reducible vertex to be worked with
    complex(dp), dimension(Nz, Nz, wperTask, NR), intent(in) :: PhiRWork
    !frequency which is supposed to be worked with atm
    integer, intent(in) :: w2
    !prefactors in PAE 
    real(dp), intent(in) :: fac1, fac2

    !loop variables
    integer :: idxQ1, R2, l1, l2, l3, l4, nu1
    !outer Bosonic arguments
    type(Indxmap) :: mapQ1 

    do idxQ1 = 1, Nb
      mapQ1 = index_bosonic_IBZ(id * Nb + idxQ1)
      do R2 = 1, NR
        do l4 = 1, Nl  
        do l3 = 1, Nl  
        do l2 = 1, Nl  
        do l1 = 1, Nl  

          do nu1 = 1, Nf

            !in case of too large frequency don't calculate value
            if(FpB(nu1, w2) > Nf) cycle
            !don't take asymptotics atm
            !if(FpB(nu1, mapQ1%iw) > Nf .or. FpB(nu1, mapQ1%iw) < 1) cycle

            F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) = &
            F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) + &
            PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + pickNu(FpB(nu1, mapQ1%iw) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2) * &
            M_phbar_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            fac1 * 1.0d0/Nx

            F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) = &
            F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) + &
            PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + pickNu(FpB(nu1, mapQ1%iw) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2) * &
            M_phbar_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            fac2 * 1.0d0/Nx

            !F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) = &
            !F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) + &
            !PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + pickNu_kernel(FpB(nu1, mapQ1%iw) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2) * &
            !M_phbar_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            !fac1 * 1.0d0/Nx

            !F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) = &
            !F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) + &
            !PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + pickNu_kernel(FpB(nu1, mapQ1%iw) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2) * &
            !M_phbar_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            !fac2 * 1.0d0/Nx
              
          end do !nu1

        end do !l1      
        end do !l2 
        end do !l3      
        end do !l4

      end do !R2
    end do !idxQ1      

  end subroutine calc_phbar_to_ph

  subroutine calc_phbar_to_ph_TR(PhiRWork, w2, fac1, fac2)

    !R-space reducible vertex to be worked with
    complex(dp), dimension(Nz, Nz, wperTask, NR), intent(in) :: PhiRWork
    !frequency which is supposed to be worked with atm
    integer, intent(in) :: w2
    !prefactors in PAE 
    real(dp), intent(in) :: fac1, fac2

    !loop variables
    integer :: idxQ1, R2, l1, l2, l3, l4, nu1
    !outer Bosonic arguments
    type(Indxmap) :: mapQ1 

    do idxQ1 = 1, Nb
      mapQ1 = index_bosonic_IBZ(id * Nb + idxQ1)
      do R2 = 1, NR
        do l4 = 1, Nl  
        do l3 = 1, Nl  
        do l2 = 1, Nl  
        do l1 = 1, Nl  

          do nu1 = 1, Nf

            if(FmB(nu1, w2) > Nf .or. FmB(nu1, w2) < 1) cycle
            !don't take asymptotics atm
            !if(mF(FpB(nu1, mapQ1%iw)) > Nf .or. mF(FpB(nu1, mapQ1%iw)) < 1) cycle

            F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) = &
            F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) + &
            signs(l3) * signs(l4) * & 
            conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + pickNu(mF(FpB(nu1, mapQ1%iw)) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2)) * &
            M_phbar_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            fac1 * 1.0d0/Nx

            F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) = &
            F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) + &
            signs(l3) * signs(l4) * & 
            conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + pickNu(mF(FpB(nu1, mapQ1%iw)) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2)) * &
            M_phbar_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            fac2 * 1.0d0/Nx

            !F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) = &
            !F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) + &
            !signs(l3) * signs(l4) * & 
            !conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + pickNu_kernel(mF(FpB(nu1, mapQ1%iw)) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2)) * &
            !M_phbar_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            !fac1 * 1.0d0/Nx

            !F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) = &
            !F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) + &
            !signs(l3) * signs(l4) * & 
            !conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + pickNu_kernel(mF(FpB(nu1, mapQ1%iw)) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2)) * &
            !M_phbar_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            !fac2 * 1.0d0/Nx

          end do !nu1

        end do !l1      
        end do !l2 
        end do !l3      
        end do !l4

      end do !R2
    end do !idxQ1      

  end subroutine calc_phbar_to_ph_TR



  subroutine calc_pp_to_ph(PhiRWork, w2, fac1, fac2)

    !R-space reducible vertex to be worked with
    complex(dp), dimension(Nz, Nz, wperTask, NR), intent(in) :: PhiRWork
    !frequency which is supposed to be worked with atm
    integer, intent(in) :: w2
    !prefactors in PAE 
    real(dp), intent(in) :: fac1, fac2

    !loop variables
    integer :: idxQ1, R2, l1, l2, l3, l4, nu1
    !outer Bosonic arguments
    type(Indxmap) :: mapQ1 

    do idxQ1 = 1, Nb
      mapQ1 = index_bosonic_IBZ(id * Nb + idxQ1)
      do R2 = 1, NR
        do l4 = 1, Nl  
        do l3 = 1, Nl  
        do l2 = 1, Nl  
        do l1 = 1, Nl  

          do nu1 = 1, Nf
            !in case of too large frequency don't calculate value
            if(BmBmF(w2, mapQ1%iw, nu1) > Nf .or. BmBmF(w2, mapQ1%iw, nu1) < 1) cycle

            F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BmBmF(w2, mapQ1%iw, nu1), idxQ1) = &
            F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BmBmF(w2, mapQ1%iw, nu1), idxQ1) + &
            PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + BmBmF(w2, mapQ1%iw, nu1), mod(w2 - 1, wperTask) + 1, R2) * &
            M_pp_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            fac1 * 1.0d0/Nx

            F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BmBmF(w2, mapQ1%iw, nu1), idxQ1) = &
            F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BmBmF(w2, mapQ1%iw, nu1), idxQ1) + &
            PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + BmBmF(w2, mapQ1%iw, nu1), mod(w2 - 1, wperTask) + 1, R2) * &
            M_pp_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            fac2 * 1.0d0/Nx
              
          end do !nu1

        end do !l1      
        end do !l2 
        end do !l3      
        end do !l4

      end do !R2
    end do !idxQ1      

  end subroutine calc_pp_to_ph



  subroutine calc_pp_to_ph_TR(PhiRWork, w2, fac1, fac2)

    !R-space reducible vertex to be worked with
    complex(dp), dimension(Nz, Nz, wperTask, NR), intent(in) :: PhiRWork
    !frequency which is supposed to be worked with atm
    integer, intent(in) :: w2
    !prefactors in PAE 
    real(dp), intent(in) :: fac1, fac2

    !loop variables
    integer :: idxQ1, R2, l1, l2, l3, l4, nu1
    !outer Bosonic arguments
    type(Indxmap) :: mapQ1 

    do idxQ1 = 1, Nb
      mapQ1 = index_bosonic_IBZ(id * Nb + idxQ1)
      do R2 = 1, NR
        do l4 = 1, Nl  
        do l3 = 1, Nl  
        do l2 = 1, Nl  
        do l1 = 1, Nl  

          do nu1 = 1, Nf
            !in case of too large frequency don't calculate value
            if(FpB(nu1, mapQ1%iw + w2 - 1) > Nf .or. FpB(nu1, mapQ1%iw + w2 - 1) < 1) cycle

            F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + mF(FpB(nu1, mapQ1%iw + w2 - 1)), idxQ1) = &
            F_d_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + mF(FpB(nu1, mapQ1%iw + w2 - 1)), idxQ1) + &
            signs(l3) * signs(l4) * &
            conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + FpB(nu1, mapQ1%iw + w2 - 1), mod(w2 - 1, wperTask) + 1, R2)) * &
            M_pp_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            fac1 * 1.0d0/Nx

            F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + mF(FpB(nu1, mapQ1%iw + w2 - 1)), idxQ1) = &
            F_m_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + mF(FpB(nu1, mapQ1%iw + w2 - 1)), idxQ1) + &
            signs(l3) * signs(l4) * &
            conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + FpB(nu1, mapQ1%iw + w2 - 1), mod(w2 - 1, wperTask) + 1, R2)) * &
            M_pp_to_ph(l1, l2, l3, l4, idxQ1, R2) * &
            fac2 * 1.0d0/Nx
              
          end do !nu1

        end do !l1      
        end do !l2 
        end do !l3      
        end do !l4

      end do !R2
    end do !idxQ1      

  end subroutine calc_pp_to_ph_TR


  subroutine calc_ph_to_pp(PhiRWork, w2, fac1, fac2)

    !R-space reducible vertex to be worked with
    complex(dp), dimension(Nz, Nz, wperTask, NR), intent(in) :: PhiRWork
    !frequency which is supposed to be worked with atm
    integer, intent(in) :: w2
    !prefactors in PAE 
    real(dp), intent(in) :: fac1, fac2

    !loop variables
    integer :: idxQ1, R2, l1, l2, l3, l4, nu1
    !outer Bosonic arguments
    type(Indxmap) :: mapQ1 

    do idxQ1 = 1, Nb
      mapQ1 = index_bosonic_IBZ(id * Nb + idxQ1)
      do R2 = 1, NR
        do l4 = 1, Nl  
        do l3 = 1, Nl  
        do l2 = 1, Nl  
        do l1 = 1, Nl  

          do nu1 = 1, Nf
            !in case of too large frequency don't calculate value
            if(BmBmF(mapQ1%iw, w2, nu1) > Nf .or. BmBmF(mapQ1%iw, w2, nu1) < 1) cycle

            F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BmBmF(mapQ1%iw, w2, nu1), idxQ1) = &
            F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BmBmF(mapQ1%iw, w2, nu1), idxQ1) + &
            PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + BmBmF(mapQ1%iw, w2, nu1), mod(w2 - 1, wperTask) + 1, R2) * &
            M_ph_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            fac1 * 1.0d0/Nx
              
            F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BmBmF(mapQ1%iw, w2, nu1), idxQ1) = &
            F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BmBmF(mapQ1%iw, w2, nu1), idxQ1) + &
            PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + BmBmF(mapQ1%iw, w2, nu1), mod(w2 - 1, wperTask) + 1, R2) * &
            M_ph_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            fac2 * 1.0d0/Nx

          end do !nu1

        end do !l1      
        end do !l2 
        end do !l3      
        end do !l4

      end do !R2
    end do !idxQ1      

  end subroutine calc_ph_to_pp



  subroutine calc_ph_to_pp_TR(PhiRWork, w2, fac1, fac2)

    !R-space reducible vertex to be worked with
    complex(dp), dimension(Nz, Nz, wperTask, NR), intent(in) :: PhiRWork
    !frequency which is supposed to be worked with atm
    integer, intent(in) :: w2
    !prefactors in PAE 
    real(dp), intent(in) :: fac1, fac2

    !loop variables
    integer :: idxQ1, R2, l1, l2, l3, l4, nu1
    !outer Bosonic arguments
    type(Indxmap) :: mapQ1 

    do idxQ1 = 1, Nb
      mapQ1 = index_bosonic_IBZ(id * Nb + idxQ1)
      do R2 = 1, NR
        do l4 = 1, Nl  
        do l3 = 1, Nl  
        do l2 = 1, Nl  
        do l1 = 1, Nl  

          do nu1 = 1, Nf
            !in fact check that bot w2 and w4 are in box
            if(BpBmF(mapQ1%iw, w2, nu1) > Nf .or. BpBmF(mapQ1%iw, w2, nu1) < 1) cycle

            F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BpBmF(mapQ1%iw, w2, nu1), idxQ1) = &
            F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BpBmF(mapQ1%iw, w2, nu1), idxQ1) + &
            signs(l3) * signs(l4) * &
            conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + FmB(nu1, mapQ1%iw + w2 - 1), mod(w2 - 1, wperTask) + 1, R2)) * &
            M_ph_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            fac1 * 1.0d0/Nx
              
            F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BpBmF(mapQ1%iw, w2, nu1), idxQ1) = &
            F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + BpBmF(mapQ1%iw, w2, nu1), idxQ1) + &
            signs(l3) * signs(l4) * &
            conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + FmB(nu1, mapQ1%iw + w2 - 1), mod(w2 - 1, wperTask) + 1, R2)) * &
            M_ph_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            fac2 * 1.0d0/Nx

          end do !nu1

        end do !l1      
        end do !l2 
        end do !l3      
        end do !l4

      end do !R2
    end do !idxQ1      

  end subroutine calc_ph_to_pp_TR



  subroutine calc_phbar_to_pp(PhiRWork, w2, fac1, fac2)

    !R-space reducible vertex to be worked with
    complex(dp), dimension(Nz, Nz, wperTask, NR), intent(in) :: PhiRWork
    !frequency which is supposed to be worked with atm
    integer, intent(in) :: w2
    !prefactors in PAE 
    real(dp), intent(in) :: fac1, fac2

    !loop variables
    integer :: idxQ1, R2, l1, l2, l3, l4, nu1
    !outer Bosonic arguments
    type(Indxmap) :: mapQ1 

    do idxQ1 = 1, Nb
      mapQ1 = index_bosonic_IBZ(id * Nb + idxQ1)
      do R2 = 1, NR
        do l4 = 1, Nl  
        do l3 = 1, Nl  
        do l2 = 1, Nl  
        do l1 = 1, Nl  

          do nu1 = 1, Nf
            !in case of too large frequency don't calculate value
            if(FpB(w2, nu1) > Nf) cycle
            !turn off asymptotics
            !if(BmBmF(mapQ1%iw, w2, nu1) > Nf .or. BmBmF(mapQ1%iw, w2, nu1) < 1) cycle

            F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) = &
            F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) + &
            PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + pickNu(BmBmF(mapQ1%iw, w2, nu1) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2) * &
            M_phbar_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            fac1 * 1.0d0/Nx
              
            F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) = &
            F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) + &
            PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + pickNu(BmBmF(mapQ1%iw, w2, nu1) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2) * &
            M_phbar_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            fac2 * 1.0d0/Nx

            !F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) = &
            !F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) + &
            !PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + pickNu_kernel(BmBmF(mapQ1%iw, w2, nu1) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2) * &
            !M_phbar_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            !fac1 * 1.0d0/Nx
            !  
            !F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) = &
            !F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FpB(nu1, w2), idxQ1) + &
            !PhiRWork((l3 - 1) * Nf + nu1, (l4 - 1) * Nf + pickNu_kernel(BmBmF(mapQ1%iw, w2, nu1) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2) * &
            !M_phbar_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            !fac2 * 1.0d0/Nx


          end do !nu1

        end do !l1      
        end do !l2 
        end do !l3      
        end do !l4

      end do !R2
    end do !idxQ1      

  end subroutine calc_phbar_to_pp




  subroutine calc_phbar_to_pp_TR(PhiRWork, w2, fac1, fac2)

    !R-space reducible vertex to be worked with
    complex(dp), dimension(Nz, Nz, wperTask, NR), intent(in) :: PhiRWork
    !frequency which is supposed to be worked with atm
    integer, intent(in) :: w2
    !prefactors in PAE 
    real(dp), intent(in) :: fac1, fac2

    !loop variables
    integer :: idxQ1, R2, l1, l2, l3, l4, nu1
    !outer Bosonic arguments
    type(Indxmap) :: mapQ1 

    do idxQ1 = 1, Nb
      mapQ1 = index_bosonic_IBZ(id * Nb + idxQ1)
      do R2 = 1, NR
        do l4 = 1, Nl  
        do l3 = 1, Nl  
        do l2 = 1, Nl  
        do l1 = 1, Nl  

          do nu1 = 1, Nf
            !in case of too large frequency don't calculate value
            if(FmB(nu1, w2) > Nf .or. FmB(nu1, w2) < 1) cycle
            !turn off asymptotics
            !if(FmB(nu1, mapQ1%iw + w2 - 1) > Nf .or. FmB(nu1, mapQ1%iw + w2 - 1) < 1) cycle

            F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) = &
            F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) + &
            signs(l3) * signs(l4) * &
            conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + pickNu(FmB(nu1, mapQ1%iw + w2 - 1) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2)) * &
            M_phbar_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            fac1 * 1.0d0/Nx
             
            F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) = &
            F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) + &
            signs(l3) * signs(l4) * &
            conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + pickNu(FmB(nu1, mapQ1%iw + w2 - 1) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2)) * &
            M_phbar_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            fac2 * 1.0d0/Nx

            !F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) = &
            !F_s_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) + &
            !signs(l3) * signs(l4) * &
            !conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + pickNu_kernel(FmB(nu1, mapQ1%iw + w2 - 1) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2)) * &
            !M_phbar_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            !fac1 * 1.0d0/Nx
            !  
            !F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) = &
            !F_t_LL((l1 - 1) * Nf + nu1, (l2 - 1) * Nf + FmB(nu1, w2), idxQ1) + &
            !signs(l3) * signs(l4) * &
            !conjg(PhiRWork((l3 - 1) * Nf + mF(nu1), (l4 - 1) * Nf + pickNu_kernel(FmB(nu1, mapQ1%iw + w2 - 1) + 2 * Nf), mod(w2 - 1, wperTask) + 1, R2)) * &
            !M_phbar_to_pp(l1, l2, l3, l4, idxQ1, R2) * &
            !fac2 * 1.0d0/Nx

          end do !nu1

        end do !l1      
        end do !l2 
        end do !l3      
        end do !l4

      end do !R2
    end do !idxQ1      

  end subroutine calc_phbar_to_pp_TR




  !function to add two frequencies
  !Fermionic plus Bosonic
  integer pure function mF(nu) result(mnu)

    integer, intent(in) :: nu

    mnu = - nu + Nf + 1

  end function mF

  !function to add two frequencies
  !Fermionic plus Bosonic
  integer pure function FpB(nu, omega) result(nuplusomega)

    integer, intent(in) :: nu
    integer, intent(in) :: omega

    nuplusomega = nu + omega - 1    

  end function FpB


  !function to add two frequencies
  !Fermionic plus Bosonic
  integer pure function BpBmF(omega1, omega2, nu) result(O1pO2mnu)

    integer, intent(in) :: omega1
    integer, intent(in) :: omega2
    integer, intent(in) :: nu

    O1pO2mnu = omega1 + omega2 - 1 + (-nu + Nf + 1) - 1

  end function BpBmF


  !function to add two frequencies
  !Fermionic plus Bosonic
  integer pure function FmB(nu, omega) result(numinusomega)

    integer, intent(in) :: nu
    integer, intent(in) :: omega

    numinusomega = nu - omega + 1

  end function FmB


  !function to add two frequencies
  !Fermionic plus Bosonic
  integer pure function BmBmF(omega1, omega2, nu) result(omega1momega2mnu)

    integer, intent(in) :: nu
    integer, intent(in) :: omega1, omega2

    omega1momega2mnu = (omega1 - 1) - (nu + omega2 - 1) + Nf + 1

  end function BmBmF




end module parquet_EQContribs

