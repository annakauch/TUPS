module parquet_chi

  !$use omp_lib

  use global_parameter
  use math_mod
  use parquet_ini
  use parquet_util
  use parquet_formfactors

  implicit none

contains

!----------------------------------------------------
!calculate Green's function

  subroutine calc_Gkw_Grt(ite, Grt, cal_chi)

    !iteration in execution of victory
    integer, intent(in)  :: ite
    !real space & time green's function
    complex(dp), intent(out) :: Grt(Nx * Ngrain, Ny * Ngrain, Nf)
    !decide whether to calculate chi0
    logical, intent(in) :: cal_chi

    !Green's function with coarse grained arguments
    !complex(dp) :: Gkw_c(Nx*Ngrain, Ny*Ngrain, Nf)
    complex(dp), allocatable :: Gkw_c(:, :, :)
 
    !loop variables
    integer :: kx, ky, nu, iTau
    !coarse grained loop variables
    integer :: kx_c, ky_c
    !actual matsubara frequency 
    complex(dp) :: freq

    !generic array index
    integer     :: idx
    !and corresponding array indices
    integer     :: idx1, idx2

    !initial bounds for mu in binary search
    real(dp) :: mu_UpperBound, mu_LowerBound

    !guess on particle number in binary search
    real(dp)    :: t1

    real(dp)    :: tau, dummy
    complex(dp), allocatable :: dummy1D(:)

    !loop variables for calc of normal Green's function
    integer :: idx_k
    type(indxmap) :: map_k
 
    ! --- variables for calculation of chi0 ---

    !loop variables
    integer :: i, j, k, i1, j1
    !integer :: idx, idx1, idx2

    !for fitting
    real(dp) :: FD1, FD2

    !FT stuff for FT-calculation of chi
    !complex(dp) :: Chi0rt_ph(Nx*Ngrain, Ny*Ngrain, Nf), Chi0rt_pp(Nx*Ngrain, Ny*Ngrain, Nf)
    complex(dp), allocatable :: Chi0rt_ph(:, :, :)
    complex(dp), allocatable :: Chi0rt_pp(:, :, :)
    !same as above just projected on formfactors
    complex(dp), allocatable :: Chi0rt_ph_FF(:, :, :, :)
    complex(dp), allocatable :: Chi0rt_pp_FF(:, :, :, :)

    !intermediate quantities
    complex(dp), allocatable :: FFG(:, :, :)

    !loop variable over 4 projections
    integer :: form_factor
    integer :: Rx, Ry, r1x, r1y, r2x, r2y
    !actual momenta
    real(dp) :: px, py

    !loop variables for coarse graining
    integer :: idx_r1x, idx_r1y, r1x_grain, r1y_grain

    !string needed for determination of operation - probably
    character(len=30) :: Mtype

    complex(dp), allocatable :: coutdata(:)

    if(.not. allocated(dummy1D)) allocate(dummy1D(Nf))
    if(.not. allocated(coutdata)) allocate(coutdata(Nf/2))

    if (.NOT. allocated(Chi0_ph)) allocate(Chi0_ph(Nred))
    if (.NOT. allocated(Chi0_pp)) allocate(Chi0_pp(Nred))

    !Green's function of arguments momentum and frequency
    IF (.NOT. allocated(Gkw)) allocate (Gkw(Nt))

    if (.not. allocated(Gkw_c)) allocate(Gkw_c(Nx*Ngrain, Ny*Ngrain, Nf))
    if (.not. allocated(Chi0rt_ph)) allocate(Chi0rt_ph(Nx*Ngrain, Ny*Ngrain, Nf))
    if (.not. allocated(Chi0rt_pp)) allocate(Chi0rt_pp(Nx*Ngrain, Ny*Ngrain, Nf))


    !calculate hartree part of selfenergy
    do ky = 1, Nx
      do kx = 1, Ny
        idx = (Ny * (kx - 1) + ky - 1) * Nf + Nf
        Sigma_H(kx, ky) = Sigma(idx)
        !if (id == master) write(*, "('Hartree energy for', 2i4, ' is', 2f12.6 )") kx, ky, Sigma_H(kx, ky)
      end do!kx
    end do!ky    

    !adjust chemical potential via binary search
    mu_UpperBound = mu + xU
    mu_LowerBound = mu - xU
   
    t1 = 3.0d0 !which is wrong -> to cycle at least once

    !perform binary search on correct chemical potential
    do while(abs(t1 - nParticle) > 1e-6)

    mu = 0.5d0 * (mu_UpperBound + mu_LowerBound)

    !$omp parallel private(idx, idx1, idx2, freq, dummy1D, tau, dummy, ky, kx, ky_c, kx_c, nu, iTau)
    !$omp do 
      do ky = 1, Ny
        do kx = 1, Nx
          do ky_c = 1, Ngrain
            do kx_c = 1, Ngrain

              do nu = 1, Nf

                idx = (Ny * (kx - 1) + ky - 1) * Nf + nu
                
                idx1= (kx - 1) * Ngrain + (kx_c - (Ngrain + 1)/2) + 1
                if (idx1 < 1) idx1 = idx1 + Ngrain*Nx
                
                idx2= (ky - 1) * Ngrain + (ky_c - (Ngrain + 1)/2) + 1
                if (idx2 < 1) idx2 = idx2 + Ngrain*Ny
                
                !freq   = dcmplx(0.0d0, Pi/beta * (2.0d0 * (nu - Nf/2 - 1) + 1.0d0))
                freq = get_freq(nu)

                Gkw_c(idx1, idx2, nu) = 1.0d0/(freq + mu - Ek_grain(kx_c, ky_c, kx, ky) - Sigma(idx))

                !subtract bare green's function
                dummy1D(nu) = (Gkw_c(idx1, idx2, nu) - 1.0d0/(freq + mu - Ek_grain(kx_c, ky_c, kx, ky))) 

              end do !nu

            ! Fourier transform from w -> tau space with the supplemented
            ! function (without non-interacting Green's function).
              do iTau = 1, Nf
                tau = beta/(Nf - 1) * (iTau - 1)
                dummy = 0.0d0

                do nu = 1, Nf
                   !freq = get_freq(nu)
                   !dummy = dummy + dble(exp(-freq * tau) * dummy1D(nu))/beta
                   !same as the two lines above
                   dummy = dummy + dble(exp(dcmplx(0.0d0, - Pi/beta * (2.0d0 * (nu - Nf/2) - 1.0d0) * tau)) * dummy1D(nu))/beta
                end do !nu

                Grt(idx1, idx2, iTau) = dummy - &
                exp((beta - tau) * (Ek_grain(kx_c, ky_c, kx, ky) - mu))/(1.0d0 + exp(beta * (Ek_grain(kx_c, ky_c, kx, ky) - mu)))

              end do !iTau   

            end do !kx_c
          end do !ky_c
        end do !kx
      end do !ky
   !$omp end do
   !$omp end parallel
     
 
      do nu = 1, Nf
         call fftf2d(Nx * Ngrain, Ny * Ngrain, Grt(1:Nx * Ngrain, 1:Ny * Ngrain, nu), C_wave_x, C_wave_y)
      end do !nu
  

      t1 = -2.0d0 * Grt(1, 1, Nf)

      !output current estimate
      if (id == master) write(*, "(' particle number is:', f12.6, ' chemical potential is', f12.6)") t1, mu

      !in case function was called to calculate chi -> dont adjust chemnical
      !potential again
      if(cal_chi) then
        exit
      end if
 
      !now perform binary search operation 
      if (t1 > nParticle) then
         mu_UpperBound = mu
      else
         mu_LowerBound = mu
      end if

    end do !while(abs(t1 - nParticle) > 1e-4)

    !calculate the actually used Green's function with newly determined mu
    do idx_k = 1, Nt

      map_k = Index_Fermionic(idx_k)
      !get actual compelx frequency
      freq = get_freq(map_k%iw)
      Gkw(idx_k) = 1.0d0/(freq - Ek(map_k%ix, map_k%iy) + mu - Sigma(idx_k))

    end do !idx_k


    ! --- chi part from old code --- 
    ! --- not in calc_chi since I don't want it at every iteration of BSE

    if(cal_chi) then
  
 
      if(id == master) write(*, *) "calculating bare bubble 'exactly'"
 
      ! --- particle-hole bubble ---
      do k = 1, Nf
        do idx1 = 1, Nx*Ngrain
          if (idx1 == 1) then
            i1 = idx1
          else
            i1 = Nx*Ngrain- idx1 +2
          end if
          
          do idx2 = 1, Ny*Ngrain
            if (idx2 == 1) then
              j1 = 1
            else
              j1 = Ny*Ngrain- idx2 +2
            end if
            Chi0rt_ph(idx1, idx2, k) = -Grt(idx1, idx2, k)*Grt(i1, j1, Nf-k+1)     
          end do
        end do
        call fftb2d(Nx*Ngrain, Ny*Ngrain, Chi0rt_ph(1:Nx*Ngrain, 1:Ny*Ngrain, k), C_wave_x, C_wave_y)
      end do

      Mtype = 'Bosonic'
      do i = 1, Nx_IBZ
        idx1 = (i-1)*Ngrain +  1
        do j = 1, i
          idx2 = (j-1)*Ngrain +  1
          do k = 1, Nf
            dummy1D(k) = Chi0rt_ph(idx1, idx2, k)
          end do
          call FDfit(Nf, dble(dummy1D), beta/dble(Nf-1), FD1, FD2)
          call nfourier(Mtype, Nf-1, Nf/2, FD1, FD2, dble(dummy1D), coutdata)
          do k = 1, Nf/2
            idx = ((i*(i-1))/2+j-1)*Nf/2 + k
            Chi0_ph(idx) = coutdata(k)
          end do
        end do
      end do
  
      ! --- particle-particle bubble ---
      do k = 1, Nf
        do idx1 = 1, Nx*Ngrain
          do idx2 = 1, Ny*Ngrain
            Chi0rt_pp(idx1, idx2, k) = - 0.5d0 * Grt(idx1, idx2, k)**2
          end do
        end do
        call fftb2d(Nx*Ngrain, Ny*Ngrain, Chi0rt_pp(1:Nx*Ngrain, 1:Ny*Ngrain, k), C_wave_x, C_wave_y)
      end do
      
      do i = 1, Nx_IBZ
        idx1 = (i-1)*Ngrain + 1
        do j = 1, i
          idx2 = (j-1)*Ngrain + 1
          do k = 1, Nf
            dummy1D(k) = Chi0rt_pp(idx1, idx2, k)    
          end do
          
          call FDfit(Nf, dble(dummy1D), beta/dble(Nf-1), FD1, FD2)
          call nfourier(Mtype, Nf-1, Nf/2, FD1, FD2, dble(dummy1D), coutdata)
          
          do k = 1, Nf/2
            idx = ((i*(i-1))/2+j-1)*Nf/2 + k
            Chi0_pp(idx) = coutdata(k)
          end do          
        end do
      end do

      !   --- end of chi part from old code

      ! do the same for bare-bubble projected on form factors

      if (.NOT. allocated(Chi0_ph_FF)) allocate(Chi0_ph_FF(Nred, 4))
      if (.NOT. allocated(Chi0_pp_FF)) allocate(Chi0_pp_FF(Nred, 4))
      if (.not. allocated(Chi0rt_ph_FF)) allocate(Chi0rt_ph_FF(Nx*Ngrain, Ny*Ngrain, Nf, 4))
      if (.not. allocated(Chi0rt_pp_FF)) allocate(Chi0rt_pp_FF(Nx*Ngrain, Ny*Ngrain, Nf, 4))

      !start by calculating needed intermediate quantities

      if(.not. allocated(FFG)) allocate(FFG(Nx * Ngrain, Ny * Ngrain, Nf))

      do form_factor = 1, 4

      FFG = dcmplx(0.0d0, 0.0d0)

      do ky = 1, Ny
        do kx = 1, Nx
          do ky_c = 1, Ngrain
            do kx_c = 1, Ngrain

              do nu = 1, Nf

                idx = (Ny * (kx - 1) + ky - 1) * Nf + nu
                
                idx1= (kx - 1) * Ngrain + (kx_c - (Ngrain + 1)/2) + 1
                if (idx1 < 1) idx1 = idx1 + Ngrain*Nx
                
                idx2= (ky - 1) * Ngrain + (ky_c - (Ngrain + 1)/2) + 1
                if (idx2 < 1) idx2 = idx2 + Ngrain*Ny
                
                !freq   = dcmplx(0.0d0, Pi/beta * (2.0d0 * (nu - Nf/2 - 1) + 1.0d0))
                freq = get_freq(nu)

                Gkw_c(idx1, idx2, nu) = 1.0d0/(freq + mu - Ek_grain(kx_c, ky_c, kx, ky) - Sigma(idx))

                !subtract bare green's function
                dummy1D(nu) = (Gkw_c(idx1, idx2, nu) - 1.0d0/(freq + mu - Ek_grain(kx_c, ky_c, kx, ky))) 

              end do !nu

            ! Fourier transform from w -> tau space with the supplemented
            ! function (without non-interacting Green's function).
              do iTau = 1, Nf
                tau = beta/(Nf - 1) * (iTau - 1)
                dummy = 0.0d0

                do nu = 1, Nf
                   !freq = get_freq(nu)
                   !dummy = dummy + dble(exp(-freq * tau) * dummy1D(nu))/beta
                   !same as the two lines above
                   dummy = dummy + dble(exp(dcmplx(0.0d0, - Pi/beta * (2.0d0 * (nu - Nf/2) - 1.0d0) * tau)) * dummy1D(nu))/beta
                end do !nu

                FFG(idx1, idx2, iTau) = (dummy - &
                exp((beta - tau) * (Ek_grain(kx_c, ky_c, kx, ky) - mu))/(1.0d0 + exp(beta * (Ek_grain(kx_c, ky_c, kx, ky) - mu)))) * &
                FF_sdwave(kx, ky, form_factor) * FF_sdwave(kx, ky, form_factor)

              end do !iTau   

            end do !kx_c
          end do !ky_c
        end do !kx
      end do !ky

      do nu = 1, Nf
         call fftf2d(Nx * Ngrain, Ny * Ngrain, FFG(1:Nx * Ngrain, 1:Ny * Ngrain, nu), C_wave_x, C_wave_y)
      end do !nu

      !do the same calculation as above just with form factor projected
      !quantities

        ! --- particle-hole bubble ---
        do k = 1, Nf
          do idx1 = 1, Nx*Ngrain
            if (idx1 == 1) then
              i1 = idx1
            else
              i1 = Nx*Ngrain- idx1 +2
            end if
            
            do idx2 = 1, Ny*Ngrain
              if (idx2 == 1) then
                j1 = 1
              else
                j1 = Ny*Ngrain- idx2 +2
              end if
              !Chi0rt_ph(idx1, idx2, k) = -Grt(idx1, idx2, k)*Grt(i1, j1, Nf-k+1)     
              Chi0rt_ph_FF(idx1, idx2, k, form_factor) = -Grt(idx1, idx2, k) * FFG(i1, j1, Nf - k + 1)
            end do
          end do
          call fftb2d(Nx*Ngrain, Ny*Ngrain, Chi0rt_ph_FF(1:Nx*Ngrain, 1:Ny*Ngrain, k, form_factor), C_wave_x, C_wave_y)
        end do
  
        Mtype = 'Bosonic'
        do i = 1, Nx_IBZ
          idx1 = (i-1)*Ngrain +  1
          do j = 1, i
            idx2 = (j-1)*Ngrain +  1
            do k = 1, Nf
              dummy1D(k) = Chi0rt_ph_FF(idx1, idx2, k, form_factor)
            end do
            call FDfit(Nf, dble(dummy1D), beta/dble(Nf-1), FD1, FD2)
            call nfourier(Mtype, Nf-1, Nf/2, FD1, FD2, dble(dummy1D), coutdata)
            do k = 1, Nf/2
              idx = ((i*(i-1))/2+j-1)*Nf/2 + k
              Chi0_ph_FF(idx, form_factor) = coutdata(k)
            end do
          end do
        end do

        ! --- particle-particle bubble ---
        do k = 1, Nf
          do idx1 = 1, Nx*Ngrain
            do idx2 = 1, Ny*Ngrain
              !Chi0rt_pp(idx1, idx2, k) = - 0.5d0 * Grt(idx1, idx2, k)**2
              Chi0rt_pp_FF(idx1, idx2, k, form_factor) = - 0.5d0 * Grt(idx1, idx2, k) * FFG(idx1, idx2, k)
            end do
          end do
          call fftb2d(Nx*Ngrain, Ny*Ngrain, Chi0rt_pp_FF(1:Nx*Ngrain, 1:Ny*Ngrain, k, form_factor), C_wave_x, C_wave_y)
        end do
        
        do i = 1, Nx_IBZ
          idx1 = (i-1)*Ngrain + 1
          do j = 1, i
            idx2 = (j-1)*Ngrain + 1
            do k = 1, Nf
              dummy1D(k) = Chi0rt_pp_FF(idx1, idx2, k, form_factor)    
            end do
            
            call FDfit(Nf, dble(dummy1D), beta/dble(Nf-1), FD1, FD2)
            call nfourier(Mtype, Nf-1, Nf/2, FD1, FD2, dble(dummy1D), coutdata)
            
            do k = 1, Nf/2
              idx = ((i*(i-1))/2+j-1)*Nf/2 + k
              Chi0_pp_FF(idx, form_factor) = coutdata(k)
            end do          
          end do
        end do

      end do !form_factor

    end if !cal_chi

  

  end subroutine calc_Gkw_Grt


  ! --- routine to calculate chis in l-space ---
  subroutine calc_chi(map_q, chi_ph_LL, chi_pp_LL)

    !bosonic argument from out loop in BSE
    type(Indxmap), intent(in) :: map_q
    !bubble in l-space to be used in BSE
    complex(dp), intent(out) :: chi_ph_LL(Nz, Nz)
    complex(dp), intent(out) :: chi_pp_LL(Nz, Nz)

    !save k-space bubble at intermediate point in calculation
    complex(dp), allocatable :: chi_k_ph(:)
    complex(dp), allocatable :: chi_k_pp(:)

    !index for calculation of k-space bubble
    integer :: idx_k
    !arguments for calculation of k-space bubble
    type(Indxmap) :: map_k, map_kpq, map_mk, map_qmk
    !coarse graining loop varibales
    integer :: kx_grain, ky_grain 

    !l-space indices for output
    integer :: idx_l1, idx_l2
    !loop varibales for calculation of l-space bubble
    integer :: l1, l2, nu, kx, ky

    !Green's functions for calculation of k-space bubble
    complex(dp) :: Gk, Gkpq, Gqmk

    if(.not. allocated(chi_k_ph)) allocate(chi_k_ph(Nt))
    if(.not. allocated(chi_k_pp)) allocate(chi_k_pp(Nt))
    chi_k_ph = 0.0d0
    chi_k_pp = 0.0d0

    do idx_k = 1, Nt
      
      map_k = Index_Fermionic(idx_k)
      !for ph bubble
      call Index_FaddB(map_k, map_q, map_kpq)
      !for pp bubble
      call index_minusF(map_k, map_mk)
      call index_FaddB(map_mk, map_q, map_qmk)

      do kx_grain = 1, Ngrain
        do ky_grain = 1, Ngrain

          !ph bubble in k-space
          Gk = get_green_coarse(kx_grain, ky_grain, map_k)
          Gkpq = get_green_coarse(kx_grain, ky_grain, map_kpq)
          chi_k_ph(idx_k) = chi_k_ph(idx_k) + Gk * Gkpq * 1.0d0/Ngrain/Ngrain

          !pp bubble in k-space
          Gqmk = get_green_coarse(Ngrain - kx_grain + 1, Ngrain - ky_grain + 1, map_qmk)
          chi_k_pp(idx_k) = chi_k_pp(idx_k) + Gk * Gqmk * 1.0d0/Ngrain/Ngrain

        end do!kx_grain
      end do!ky_grain

    end do !idx_k

    !now transform bubble into l-space
    chi_ph_LL = 0.0d0
    chi_pp_LL = 0.0d0

    do l2 = 1, Nl
      do l1 = 1, Nl
        do nu = 1, Nf
          do ky = 1, Ny
            do kx = 1, Nx
  
              idx_l1 = (l1 - 1) * Nf + nu        
              idx_l2 = (l2 - 1) * Nf + nu        
    
              chi_ph_LL(idx_l1, idx_l2) = chi_ph_LL(idx_l1, idx_l2) + &
                                          chi_k_ph(((kx - 1) * Ny + ky - 1) * Nf + nu) * &
                                          FF(kx, ky, l1) * FF(kx, ky, l2)

              chi_pp_LL(idx_l1, idx_l2) = chi_pp_LL(idx_l1, idx_l2) + &
                                          chi_k_pp(((kx - 1) * Ny + ky - 1) * Nf + nu) * &
                                          FF(kx, ky, l1) * FF(kx, ky, l2)

            end do !kx
          end do !ky
        end do !nu
      end do !l1
    end do !l2

  end subroutine calc_chi


end module parquet_chi

