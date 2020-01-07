module parquet_sus

!module that holds functions to calculate physical quantities

  use parquet_ini
  use parquet_util
  use parquet_kernel

  implicit none

  integer, allocatable :: pickNu(:)
  integer, allocatable :: pickNu_K(:)

  contains


  !this calculates the complete susceptibility
  subroutine calc_sus()

    if(.not. allocated(sus_d_IBZ))  allocate(sus_d_IBZ(Nred, 4))
    if(.not. allocated(sus_m_IBZ))  allocate(sus_m_IBZ(Nred, 4))
    if(.not. allocated(sus_s_IBZ))  allocate(sus_s_IBZ(Nred, 4))
    if(.not. allocated(sus_t_IBZ))  allocate(sus_t_IBZ(Nred, 4))

    sus_d_IBZ = 0.0d0
    sus_m_IBZ = 0.0d0
    sus_s_IBZ = 0.0d0
    sus_t_IBZ = 0.0d0

    call fillpickNu_range()

    !calculate d and m susceptibility
    call calc_sus_ph()


    !calculate pp susceptibility
    call calc_sus_pp()

  end subroutine calc_sus

! --------------------------------------------------

  !this calculates the susceptibitlity from the ph vertices
  !this also adds the bare vertices - since this is convenient at this point
  subroutine calc_sus_ph()
  
    !loop variables of arguments and ...
    integer :: idx_q1, idx_l1, idx_l2, idx_l1_l, idx_l2_l
    integer :: l1, l2, l10, l20

    !holds product of 2 Green's functions and a formfactor
    complex(dp), allocatable :: GG(:, :, :)
 
    !susceptibility for each task that will be put together at the end
    complex(dp), allocatable :: sus_aux_d(:, :), sus_aux_m(:, :)

    !maps to pass to kernel function
    type(indxmap_L) :: map_l1, map_l2
    type(indxmap) :: map_q1

    !external form-factor argument
    integer :: form_factor

    !auxillary susceptibility to assembel whole in the end on each task
    if(.not. allocated(sus_aux_d)) allocate(sus_aux_d(Nb, 4)) 
    if(.not. allocated(sus_aux_m)) allocate(sus_aux_m(Nb, 4)) 

    !save product of two Green's functions
    if(.not. allocated(GG)) allocate(GG(Nl * Nf * (2 * f_range + 1), Nb, 4))
    GG = 0.0d0

    call calc_GG_ph(GG)
  
    sus_aux_d = 0.0d0
    sus_aux_m = 0.0d0

    do form_factor = 1, 4
    do idx_q1 = 1, Nb

      map_q1 = index_bosonic_IBZ(id * Nb + idx_q1)
      !if(map_q1%iw .ne. 1) cycle

      do l2 = 1, Nl
      do l20 = 1, Nf

        idx_l2 = (l2 - 1) * Nf + l20
        idx_l2_l = (l2 - 1) * (2 * f_range + 1) * Nf + l20 + f_range * Nf

        do l1 = 1, Nl
        do l10 = 1, Nf

          idx_l1 = (l1 - 1) * Nf + l10
          idx_l1_l = (l1 - 1) * (2 * f_range + 1) * Nf + l10 + f_range * Nf

          !add density and magnetic part
          sus_aux_d(idx_q1, form_factor) = sus_aux_d(idx_q1, form_factor) + &
                              (F_d_LL(idx_l1, idx_l2, idx_q1)) * &
                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          !sus_aux_d(idx_q1, form_factor) = sus_aux_d(idx_q1, form_factor) - &
          !                    G_d_LL(idx_l1, idx_l2, idx_q1) * &  
          !                    GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
          !                    1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          !if(l1 == 1 .and. l2 == 1) then
          !sus_aux_d(idx_q1, form_factor) = sus_aux_d(idx_q1, form_factor) - &
          !                    L_d(l10, l20, map_q1%iw) *  1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
          !                    GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
          !                    1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
          !end if !l1 == l2 == 1
          
          sus_aux_m(idx_q1, form_factor) = sus_aux_m(idx_q1, form_factor) + &
                              (F_m_LL(idx_l1, idx_l2, idx_q1)) * &
                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          
        end do !l1
        end do !l10
      end do !l2
      end do !l20
  
    end do !idx_q1
    end do !form_factor

    ! - add asymptotics -

    do form_factor = 1, 4
    do idx_q1 = 1, Nb

      map_q1 = index_bosonic_IBZ(id * Nb + idx_q1)
      !if(map_q1%iw .ne. 1) cycle

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
            sus_aux_d(idx_q1, form_factor) = sus_aux_d(idx_q1, form_factor) + &
                                xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                                GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                                1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
            
            sus_aux_m(idx_q1, form_factor) = sus_aux_m(idx_q1, form_factor) - &
                                xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1))* &
                                GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                                1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          end if ! l1 == l2 == 1         


          sus_aux_d(idx_q1, form_factor) = sus_aux_d(idx_q1, form_factor) + &
                              G_d_LL((l1 - 1) * Nf + pickNu(l10 + f_range * Nf), (l2 - 1) * Nf + pickNu(l20 + f_range * Nf), idx_q1) * &
                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
                              
          sus_aux_m(idx_q1, form_factor) = sus_aux_m(idx_q1, form_factor) + &
                              G_m_LL((l1 - 1) * Nf + pickNu(l10 + f_range * Nf), (l2 - 1) * Nf + pickNu(l20 + f_range * Nf), idx_q1) * &
                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)


        end do !l1
        end do !l10
      end do !l2
      end do !l20
  
    end do !idx_q1
    end do !form_factor

    do form_factor = 1, 4
      !assembel the IBZ vertex from contributions of each task
      call assemble_IBZ(sus_aux_d(:, form_factor), sus_d_IBZ(:, form_factor))
      call assemble_IBZ(sus_aux_m(:, form_factor), sus_m_IBZ(:, form_factor))
    end do !form_factor

  end subroutine calc_sus_ph


  !this calculates the susceptibitlity from the ph vertices
  !this also adds the bare vertices - since this is convenient at this point
  subroutine calc_sus_pp()
  
    !loop variables of arguments and ...
    integer :: idx_q1, idx_l1, idx_l2, idx_l1_l, idx_l2_l
    integer :: l1, l2, l10, l20

    !holds product of 2 Green's functions and a formfactor
    complex(dp), allocatable :: GG(:, :, :)
 
    !susceptibility for each task that will be put together at the end
    complex(dp), allocatable :: sus_aux_s(:, :), sus_aux_t(:, :)

    !maps to pass to kernel function
    type(indxmap_L) :: map_l1, map_l2
    type(indxmap) :: map_q1

    !external form-factor argument
    integer :: form_factor

    !auxillary susceptibility to assembel whole in the end on each task
    if(.not. allocated(sus_aux_s)) allocate(sus_aux_s(Nb, 4)) 
    if(.not. allocated(sus_aux_t)) allocate(sus_aux_t(Nb, 4)) 

    !save product of two Green's functions
    if(.not. allocated(GG)) allocate(GG(Nl * Nf * (2 * f_range + 1), Nb, 4))
    GG = 0.0d0

    call calc_GG_pp(GG)
  
    sus_aux_s = 0.0d0
    sus_aux_t = 0.0d0

    do form_factor = 1, 4
    do idx_q1 = 1, Nb

      map_q1 = index_bosonic_IBZ(id * Nb + idx_q1)
      !if(map_q1%iw .ne. 1) cycle

      do l2 = 1, Nl
      do l20 = 1, Nf

        idx_l2 = (l2 - 1) * Nf + l20
        idx_l2_l = (l2 - 1) * (2 * f_range + 1) * Nf + l20 + f_range * Nf

        do l1 = 1, Nl
        do l10 = 1, Nf

          idx_l1 = (l1 - 1) * Nf + l10
          idx_l1_l = (l1 - 1) * (2 * f_range + 1) * Nf + l10 + f_range * Nf

          !add density and magnetic part
          sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) + &
                              F_s_LL(idx_l1, idx_l2, idx_q1) * &
                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          !sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) - &
          !                    G_s_LL(idx_l1, idx_l2, idx_q1) * &  
          !                    GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
          !                    1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          !if(l1 == 1 .and. l2 == 1) then
          !sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) - &
          !                    L_s(l10, l20, map_q1%iw) *  1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
          !                    GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
          !                    1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
          !end if !l1 == l2 == 1
          
          sus_aux_t(idx_q1, form_factor) = sus_aux_t(idx_q1, form_factor) + &
                              (F_t_LL(idx_l1, idx_l2, idx_q1)) * &
                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          !sus_aux_t(idx_q1, form_factor) = sus_aux_t(idx_q1, form_factor) - &
          !                    G_t_LL(idx_l1, idx_l2, idx_q1) * &
          !                    GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
          !                    1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          !if(l1 == 1 .and. l2 == 1) then
          !sus_aux_t(idx_q1, form_factor) = sus_aux_t(idx_q1, form_factor) - &
          !                    L_t(l10, l20, map_q1%iw) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
          !                    GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
          !                    1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          !end if !l1 == l2 == 1
          
        end do !l1
        end do !l10
      end do !l2
      end do !l20
  
    end do !idx_q1
    end do !form_factor

    ! - add asymptotics -

    do form_factor = 1, 4
    do idx_q1 = 1, Nb

      map_q1 = index_bosonic_IBZ(id * Nb + idx_q1)
      !if(map_q1%iw .ne. 1) cycle

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
            sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) + &
                                2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                                GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                                1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
            
            sus_aux_t(idx_q1, form_factor) = sus_aux_t(idx_q1, form_factor) - &
                                0.0d0/(FF(1, 1, 1) * FF(1, 1, 1))* &
                                GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                                1.0d0/(Nx * Nx * Ny * Ny * beta * beta)

          end if ! l1 == l2 == 1         


          sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) + &
                              G_s_LL((l1 - 1) * Nf + pickNu(l10 + f_range * Nf), (l2 - 1) * Nf + pickNu(l20 + f_range * Nf), idx_q1) * &
                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
                              
          sus_aux_t(idx_q1, form_factor) = sus_aux_t(idx_q1, form_factor) + &
                              G_t_LL((l1 - 1) * Nf + pickNu(l10 + f_range * Nf), (l2 - 1) * Nf + pickNu(l20 + f_range * Nf), idx_q1) * &
                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)


        end do !l1
        end do !l10
      end do !l2
      end do !l20
  
    end do !idx_q1
    end do !form_factor

    do form_factor = 1, 4
      !assembel the IBZ vertex from contributions of each task
      call assemble_IBZ(sus_aux_s(:, form_factor), sus_s_IBZ(:, form_factor))
      call assemble_IBZ(sus_aux_t(:, form_factor), sus_t_IBZ(:, form_factor))
    end do !form_factor

  end subroutine calc_sus_pp



!! --------------------------------------------------
!
!  !this calculates the susceptibitlity from the ph vertices
!  !this also adds the bare vertices - since this is convenient at this point
!  subroutine calc_sus_pp()
!  
!    !loop variables of arguments and ...
!    integer :: idx_q1, idx_l1, idx_l2, idx_l1_l, idx_l2_l
!    integer :: l1, l2, l10, l20
!
!    !holds product of 2 Green's functions and a formfactor
!    complex(dp), allocatable :: GG(:, :, :)
! 
!    !susceptibility for each task that will be put together at the end
!    complex(dp), allocatable :: sus_aux_s(:, :), sus_aux_t(:, :)
!
!    !maps to pass to kernel function
!    type(indxmap_L) :: map_l1, map_l2
!    type(indxmap) :: map_q1
!
!    !external form-factor argument
!    integer :: form_factor
!
!
!    !auxillary susceptibility to assembel whole in the end on each task
!    if(.not. allocated(sus_aux_s)) allocate(sus_aux_s(Nb, 4)) 
!    if(.not. allocated(sus_aux_t)) allocate(sus_aux_t(Nb, 4)) 
!
!    !save product of two Green's functions
!    if(.not. allocated(GG)) allocate(GG(Nl * Nf * (2 * f_range + 1), Nb, 4))
!    GG = 0.0d0
!
!    call calc_GG_pp(GG)
! 
!    sus_aux_s = 0.0d0
!    sus_aux_t = 0.0d0
!
!    do form_factor = 1, 4
!    do idx_q1 = 1, Nb
!
!      map_q1 = index_bosonic_IBZ(id * Nb + idx_q1)
!      !if(map_q1%iw .ne. 1) cycle
!
!      do l2 = 1, Nl
!      do l20 = 1, Nf
!
!        idx_l2 = (l2 - 1) * Nf + l20
!        idx_l2_l = (l2 - 1) * (2 * f_range + 1) * Nf + l20 + f_range * Nf
!
!        do l1 = 1, Nl
!        do l10 = 1, Nf
!
!          idx_l1 = (l1 - 1) * Nf + l10
!          idx_l1_l = (l1 - 1) * (2 * f_range + 1) * Nf + l10 + f_range * Nf
!
!          !add density and magnetic part
!          sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) + &
!                              F_s_LL(idx_l1, idx_l2, idx_q1) * &
!                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
!                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
!
!          sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) - &
!                              G_s_LL(idx_l1, idx_l2, idx_q1) * &
!                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
!                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
!
!          if(l1 == 1 .and. l2 == 1) then
!
!          sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) - &
!                              L_s(l10, l20, map_q1%iw) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
!                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
!                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
!
!          end if !l1 == l2 == 1
!          
!          sus_aux_t(idx_q1, form_factor) = sus_aux_t(idx_q1, form_factor) + &
!                              F_t_LL(idx_l1, idx_l2, idx_q1) * &
!                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
!                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
!
!          
!        end do !l1
!        end do !l10
!      end do !l2
!      end do !l20
!  
!    end do !idx_q1
!    end do !form_factor
!
!    ! - add asymptotics -
!
!    do form_factor = 1, 4
!    do idx_q1 = 1, Nb
!
!      map_q1 = index_bosonic_IBZ(id * Nb + idx_q1)
!      !if(map_q1%iw .ne. 1) cycle
!
!      do l2 = 1, Nl
!      do l20 = - f_range * Nf + 1, (f_range + 1) * Nf
!
!        !idx_l2 = (l2 - 1) * Nf + l20
!        idx_l2_l = (l2 - 1) * (2 * f_range + 1) * Nf + l20 + f_range * Nf
!        map_l2 = indxmap_L(l2, l20)
!
!        do l1 = 1, Nl
!        do l10 = - f_range * Nf + 1, (f_range + 1) * Nf
!
!          !if both frequencies are inside the box contribution has already been taken
!          if((l10 > 0 .and. l10 < Nf + 1) .and. (l20 > 0 .and. l20 < Nf + 1)) cycle
!
!          !idx_l1 = (l1 - 1) * Nf + l10
!          idx_l1_l = (l1 - 1) * (2 * f_range + 1) * Nf + l10 + f_range * Nf
!          map_l1 = indxmap_L(l1, l10)
!
!
!          !add contribution from bare U
!          if((l1 == 1) .and. (l2 == 1)) then
!
!            !contribution from bare U
!            sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) + &
!                                xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
!                                GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
!                                1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
!            
!            sus_aux_t(idx_q1, form_factor) = sus_aux_t(idx_q1, form_factor) - &
!                                xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1))* &
!                                GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
!                                1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
!
!          end if ! l1 == l2 == 1         
!
!          sus_aux_s(idx_q1, form_factor) = sus_aux_s(idx_q1, form_factor) + &
!                              G_s_LL((l1 - 1) * Nf + pickNu(l10 + f_range * Nf), (l2 - 1) * Nf + pickNu(l20 + f_range * Nf), idx_q1) * &
!                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
!                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
!                              
!          sus_aux_t(idx_q1, form_factor) = sus_aux_t(idx_q1, form_factor) + &
!                              G_t_LL((l1 - 1) * Nf + pickNu(l10 + f_range * Nf), (l2 - 1) * Nf + pickNu(l20 + f_range * Nf), idx_q1) * &
!                              GG(idx_l1_l, idx_q1, form_factor) * GG(idx_l2_l, idx_q1, form_factor) * &
!                              1.0d0/(Nx * Nx * Ny * Ny * beta * beta)
!
!        end do !l1
!        end do !l10
!      end do !l2
!      end do !l20
!  
!    end do !idx_q1
!    end do !form_factor
!
!
!    do form_factor = 1, 4
!      !assembel the IBZ vertex from contributions of each task
!      call assemble_IBZ(sus_aux_s(:, form_factor), sus_s_IBZ(:, form_factor))
!      call assemble_IBZ(sus_aux_t(:, form_factor), sus_t_IBZ(:, form_factor))
!    end do !form_factor
!
!  end subroutine calc_sus_pp


  ! --------------------------------------

  !subroutine to calculate the product of two Green's functions
  subroutine calc_GG_ph(GG)

    !product of two Green's functions to be outputted
    complex(dp), intent(out), dimension(Nl * Nf * (2 * f_range + 1), Nb, 4) :: GG
    
    !loop variable for sum
    integer :: idx_q, idx_l, kx, ky
    integer :: l, l0
    integer :: kx_grain, ky_grain

    !arguments to put into FF
    type(indxmap) :: map_k, map_q, map_kq

    !Green's functions to be multiplied
    complex(dp) :: Gk, Gkq

    !which shape does the external FF-argument have
    integer :: form_factor

    do form_factor = 1, 4
    do idx_q = 1, Nb
      map_q = Index_bosonic_IBZ(id * Nb + idx_q)
      do l = 1, Nl
        do l0 = -f_range * Nf + 1, (f_range + 1) * Nf

        idx_l = (l - 1) * (2 * f_range + 1) * Nf + l0 + f_range * Nf

        do ky = 1, Ny
          do kx = 1, Nx

            map_k = indxmap(kx, ky, l0)
            call Index_FaddB(map_k, map_q, map_kq)

            do ky_grain = 1, Ngrain
              do kx_grain = 1, Ngrain

                Gk = get_green_coarse(kx_grain, ky_grain, map_k)
                Gkq = get_green_coarse(kx_grain, ky_grain, map_kq)

                GG(idx_l, idx_q, form_factor) = GG(idx_l, idx_q, form_factor) + &
                                                FF_sdwave(kx, ky, form_factor) * Gk * Gkq * FF(kx, ky, l) * &
                                                1.0d0/Ngrain/Ngrain

              end do !kx_grain
            end do !ky_grain


          end do !kx
        end do !ky    
      end do !l
      end do !l0
    end do !idx_q
    end do !form_factor


  end subroutine calc_GG_ph

  ! --------------------------------------

  !subroutine to calculate the product of two Green's functions
  subroutine calc_GG_pp(GG)

    !product of two Green's functions to be outputted
    complex(dp), intent(out), dimension(Nl * Nf * (2 * f_range + 1), Nb, 4) :: GG
    
    !loop variable for sum
    integer :: idx_q, idx_l, kx, ky
    integer :: l, l0
    integer :: kx_grain, ky_grain

    !arguments to put into FF
    type(indxmap) :: map_k, map_q, map_qmk, map_mk

    !Green's functions to be multiplied
    complex(dp) :: Gk, Gqmk

    !which shape does the external FF-argument have
    integer :: form_factor

    do form_factor = 1, 4
      do idx_q = 1, Nb
        map_q = Index_bosonic_IBZ(id * Nb + idx_q)
        do l = 1, Nl
          do l0 = -f_range * Nf + 1, (f_range + 1) * Nf

          idx_l = (l - 1) * (2 * f_range + 1) * Nf + l0 + f_range * Nf

          do ky = 1, Ny
            do kx = 1, Nx

              map_k = indxmap(kx, ky, l0)
              call index_minusF(map_k, map_mk) 
              call index_FaddB(map_mk, map_q, map_qmk) 

              do ky_grain = 1, Ngrain
                do kx_grain = 1, Ngrain

                  Gk = get_green_coarse(kx_grain, ky_grain, map_k)
                  Gqmk = get_green_coarse(Ngrain - kx_grain + 1, Ngrain - ky_grain + 1, map_qmk)

                  GG(idx_l, idx_q, form_factor) = GG(idx_l, idx_q, form_factor) + & 
                                   FF_sdwave(kx, ky, form_factor) * Gk * Gqmk * FF(kx, ky, l) * &
                                   1.0d0/Ngrain/Ngrain

                end do !kx_grain
              end do !ky_grain


            end do !kx
          end do !ky    
        end do !l
        end do !l0
      end do !idx_q
    end do !form_factor

  end subroutine calc_GG_pp


! ------------------------------

  !calculate double occupancie
  subroutine double_occupation()

    integer :: idxQ 

    double_occ_1P = 0.0d0
    double_occ_2P = 0.0d0

    double_occ_1P = 1.0d0/(beta * Nx * Ny) * sum(Sigma * Gkw)

    double_occ_2P = 0.5d0/(beta * Nx * Ny) * sum(sus_d - sus_m) + 1.0d0/beta * (NParticle/2.0d0) * (NParticle/2.0d0)

  end subroutine double_occupation


  ! ----------------------------------------------------


  !assembels whole IBZ vertrex from the contributions on individual nodes
  subroutine assemble_IBZ(in_node, out_IBZ)

    !bosonic arguments helt by each node
    complex(dp), dimension(Nb), intent(in) :: in_node
    !array in IBZ
    complex(dp), dimension(Nred), intent(inout) :: out_IBZ

    !loop variable to loop over tasks
    integer :: task
    !temporary in communication
    complex(dp), allocatable :: arr_temp(:)
   
    !assembel whole suszeptibility on master
    if(id .NE. master) then
  
      call MPI_SEND(in_node, Nb, MPI_DOUBLE_COMPLEX, master, id, MPI_COMM_WORLD, rc)
  
    else !if id == master
 
      !only needed in mpi communication
      if(.not. allocated(arr_temp)) allocate(arr_temp(Nb)) 

      !first master adds its own contribution
      out_IBZ(1 : Nb) = out_IBZ(1:Nb) + in_node(1 : Nb)

      !then it recieves information from other tasks
      do task = 1, ntasks - 1
    
        !recv susceptibility part from other tasks
        call MPI_RECV(arr_temp, Nb, MPI_DOUBLE_COMPLEX, task, task, MPI_COMM_WORLD, stat, rc)
      
        !add it to susceptibility with appropriate indexing
        out_IBZ(task * Nb + 1 : (task + 1) * Nb) = out_IBZ(task * Nb + 1 : (task + 1) * Nb) + arr_temp(1:Nb)
  
      end do !task

    end if !id 

  end subroutine assemble_IBZ


  subroutine fillpickNu_range()

    integer :: nu

    !take Nf + 2*Nf on both sides - should be enough to cover bo - should be
    !enough to cover boxx
    if(.not. allocated(pickNu)) allocate(pickNu((2 * f_range + 1) * Nf))
    if(.not. allocated(pickNu_K)) allocate(pickNu_K((2 * f_range + 1) * Nf))

    do nu = 1, (2 * f_range + 1) * Nf
  
      if(nu .le. f_range * Nf) then
        pickNu(nu) = 1
        pickNu_K(nu) = Nf
      else if(f_range * Nf < nu .and. nu < (f_range + 1) * Nf + 1) then
        pickNu(nu) = nu - f_range * Nf
        pickNu_K(nu) = nu - f_range * Nf
      else 
        pickNu(nu) = Nf
        pickNu_K(nu) = Nf
      end if !nu

    end do !nu
  end subroutine fillpickNu_range


end module parquet_sus
