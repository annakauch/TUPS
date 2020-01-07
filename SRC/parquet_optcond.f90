module parquet_optcond

!In this module the optical conductivity is calculated

  use parquet_ini
  use parquet_util
  use parquet_kernel
  use parquet_sus

  implicit none

  contains


  !this calculates the complete susceptibility
  subroutine calc_opt_cond()

    complex(dp), allocatable :: opt_cond_L(:) 
    complex(dp), allocatable :: opt_cond_d(:)
    complex(dp), allocatable :: opt_cond_m(:)
    complex(dp), allocatable :: opt_cond_s(:)
    complex(dp), allocatable :: opt_cond_t(:)

    if(.not. allocated(opt_cond_L)) allocate(opt_cond_L(Nf/2))
    !if(.not. allocated(opt_cond_d)) allocate(opt_cond_d(Nf/2))
    !if(.not. allocated(opt_cond_m)) allocate(opt_cond_m(Nf/2))
    !if(.not. allocated(opt_cond_s)) allocate(opt_cond_s(Nf/2))
    !if(.not. allocated(opt_cond_t)) allocate(opt_cond_t(Nf/2))

    opt_cond_L = dcmplx(0.0d0, 0.0d0)
    !opt_cond_d = dcmplx(0.0d0, 0.0d0)
    !opt_cond_m = dcmplx(0.0d0, 0.0d0)
    !opt_cond_s = dcmplx(0.0d0, 0.0d0)
    !opt_cond_t = dcmplx(0.0d0, 0.0d0)


    call calc_opt_cond_L(opt_cond_L)
    !call calc_opt_cond_d_native(opt_cond_d)
    !call calc_opt_cond_d_nonnative(opt_cond_d)
    !call calc_opt_cond_m(opt_cond_m)
    !call calc_opt_cond_s(opt_cond_s)
    !call calc_opt_cond_t(opt_cond_t)

    if(id == master) then
      write(*, *) 'writing opt_cond', opt_cond_L(1), ' ', opt_cond_L(2) 
      write(*, *) 'writing sus_d', sus_d_IBZ(1, 3), ' ', sus_d_IBZ(2, 3) 
    end if

  end subroutine calc_opt_cond

! --------------------------------------------------

  subroutine calc_opt_cond_L(opt_cond_L_IBZ)
  
    complex(dp), intent(inout) :: opt_cond_L_IBZ(:)
  
    !loop variables of arguments and ...
    integer :: idx_q1, idx_l1, idx_l2, idx_l1_l, idx_l2_l
    integer :: l1, l2, l10, l20
  
    !holds product of 2 Green's functions and a formfactor
    complex(dp), allocatable :: GG(:, :, :)
   
    !maps to pass to kernel function
    type(indxmap_L) :: map_l1, map_l2
    type(indxmap) :: map_q1
  
    !save product of two Green's functions
    if(.not. allocated(GG)) allocate(GG(Nl * Nf * (2 * f_range + 1), Nb, 4))
    GG = 0.0d0
  
    call calc_GG_ph(GG)
    
    do idx_q1 = 1, Nb
  
      map_q1 = index_bosonic_IBZ(id * Nb + idx_q1)
      !only take spatially uniform part
      if(map_q1%ix .ne. 1 .or. map_q1%iy .ne. 1) cycle
  
      do l2 = 1, Nl
      do l20 = 1, Nf
  
        idx_l2 = (l2 - 1) * Nf + l20
        idx_l2_l = (l2 - 1) * (2 * f_range + 1) * Nf + l20 + f_range * Nf
  
        do l1 = 1, Nl
        do l10 = 1, Nf
  
          idx_l1 = (l1 - 1) * Nf + l10
          idx_l1_l = (l1 - 1) * (2 * f_range + 1) * Nf + l10 + f_range * Nf
  
          !add density and magnetic part
          opt_cond_L_IBZ(map_q1%iw) = opt_cond_L_IBZ(map_q1%iw) + &
                              (F_d_LL(idx_l1, idx_l2, idx_q1)) * &
                              !(L_d_LL(idx_l1, idx_l2, idx_q1)) * &
                              GG(idx_l1_l, idx_q1, 3) * GG(idx_l2_l, idx_q1, 3) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
          
        end do !l1
        end do !l10
      end do !l2
      end do !l20
    
    end do !idx_q1
  
    ! - add asymptotics -
  
    do idx_q1 = 1, Nb
  
      map_q1 = index_bosonic_IBZ(id * Nb + idx_q1)
      !only take spatially uniform part
      if(map_q1%ix .ne. 1 .or. map_q1%iy .ne. 1) cycle
  
      do l2 = 1, Nl
      do l20 = - f_range * Nf + 1, (f_range + 1) * Nf
  
        !idx_l2 = (l2 - 1) * Nf + l20
        idx_l2_l = (l2 - 1) * (2 * f_range + 1) * Nf + l20 + f_range * Nf
        map_l2 = indxmap_L(l2, l20)
  
        do l1 = 1, Nl
        do l10 = - f_range * Nf + 1, (f_range + 1) * Nf
  
          !if both frequencies are inside the box contribution has already been taken
          if((l10 > 0 .and. l10 < Nf + 1) .and. (l20 > 0 .and. l20 < Nf + 1)) cycle
  
          !idx_l1 = (l1 - 1) * Nf + l10
          idx_l1_l = (l1 - 1) * (2 * f_range + 1) * Nf + l10 + f_range * Nf
          map_l1 = indxmap_L(l1, l10)
  
  
          !add contribution from bare U
          if((l1 == 1) .and. (l2 == 1)) then
  
            !contribution from bare U
            opt_cond_L_IBZ(map_q1%iw) = opt_cond_L_IBZ(map_q1%iw) + &
                                xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                                GG(idx_l1_l, idx_q1, 3) * GG(idx_l2_l, idx_q1, 3) * &
                                1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
            
          end if ! l1 == l2 == 1         
  
  
          opt_cond_L_IBZ(idx_q1) = opt_cond_L_IBZ(idx_q1) + &
                              G_d_LL((l1 - 1) * Nf + pickNu(l10 + f_range * Nf), (l2 - 1) * Nf + pickNu(l20 + f_range * Nf), idx_q1) * &
                              GG(idx_l1_l, idx_q1, 3) * GG(idx_l2_l, idx_q1, 3) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
                              
  
        end do !l1
        end do !l10
      end do !l2
      end do !l20
    
    end do !idx_q1

  end subroutine calc_opt_cond_L


end module parquet_optcond
