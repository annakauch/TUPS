module parquet_selfenergy

  use global_parameter
  use math_mod
  use parquet_ini
  use parquet_util
  use parquet_kernel
  use parquet_formfactors
  use parquet_equation

  !$use omp_lib

  implicit none

  private

  public :: self_energy

contains

  subroutine self_energy(ite, Grt)

    integer, intent(in)     :: ite
    !real space & time green's function
    complex(dp), intent(in) :: Grt(Nx * Ngrain, Ny * Ngrain, Nf)

    !holds intermediate result in loops
    complex(dp) sig_inter

    complex(dp) sig_U
    complex(dp) sig_ph
    complex(dp) sig_pb
    complex(dp) sig_pp

    !arguments for vertices and green's functions
    integer :: l3
    !loop indices
    integer :: idx_k, idx_q
    !loop variables for IBZ
    integer :: kx, ky, k_0
    integer :: idx_l4, idx_l3, idx_l3_mnu
    !maps
    type(Indxmap) map_k, map_q, map_kq, map_kmq, map_mk, map_qmk
    type(Indxmap_L) map_l3, map_l4, map_l3_mnu

    !for using vertex asymptotics
    integer :: l4_0 !<- loop variable for frequency in an enlarged box
    integer :: l4   !<- simple l-value

    !IBZ stuff
    !loop variable for symmetries
    integer :: is
    !symmetriezed q-map
    type(Indxmap) map_q_sym 
    !for getting the whole BZ from IBZ
    logical, dimension(8) :: symm_list

    !chi-like product of green's functions needed in l-space schwinger dyson
    !equaiton
    complex(dp), dimension(:, :), allocatable :: Chi_aux_ph
    complex(dp), dimension(:, :), allocatable :: Chi_aux_pp
    !chi-functions in enlarged frequency box
    complex(dp), dimension(:, :), allocatable :: Chi_aux_ph_L
    complex(dp), dimension(:, :), allocatable :: Chi_aux_pp_L

    !Sigma_V - sigma calculated on one node - will the be reduced
    !Sigma Reduced - memory adress to hold data from mpi call
    complex(dp), dimension(:), allocatable :: Sigma_V, Sigma_Reduced

    !split up contributions
    complex(dp), dimension(:), allocatable :: SigV_U, SigRed_U
    complex(dp), dimension(:), allocatable :: SigV_ph, SigRed_ph 
    complex(dp), dimension(:), allocatable :: SigV_pb, SigRed_pb
    complex(dp), dimension(:), allocatable :: SigV_pp, SigRed_pp

    complex(dp) :: Gkq, Gkmq, Gqmk
    !second order contribution in Schwinger dyson equation caculated via FT
    complex(dp), allocatable :: order_2nd(:)

    !to replace local part of selfenergy by DMFT
    complex(dp) :: sig_loc
    integer :: nu

    if (.NOT. allocated(Sigma_V)) allocate (Sigma_V(2 * Nred))

    if (.NOT. allocated(SigV_U)) allocate (SigV_U(2 * Nred))
    if (.NOT. allocated(SigV_ph)) allocate (SigV_ph(2 * Nred))
    if (.NOT. allocated(SigV_pb)) allocate (SigV_pb(2 * Nred))
    if (.NOT. allocated(SigV_pp)) allocate (SigV_pp(2 * Nred))


    if (.NOT. allocated(Chi_aux_ph)) allocate (Chi_aux_ph(Nz, Nb))
    if (.NOT. allocated(Chi_aux_pp)) allocate (Chi_aux_pp(Nz, Nb))
    if (.NOT. allocated(Chi_aux_ph_L)) allocate (Chi_aux_ph_L(Nl * (2 * f_range + 1) * Nf, Nb))
    if (.NOT. allocated(Chi_aux_pp_L)) allocate (Chi_aux_pp_L(Nl * (2 * f_range + 1) * Nf, Nb))
    if (.NOT. allocated(order_2nd)) allocate(order_2nd(Nx * Ny * Nf))


    !caculate auxillary chis that appear in l-space schwinger dyson
    call calc_chi_aux(chi_aux_ph, chi_aux_pp, chi_aux_ph_L, chi_aux_pp_L)

    !evaluate SD with reducible vertices which are saved in G after evaluation
    !of BS-equations
    ! -> unbiased treatment compared to using F

    Sigma = (0.0d0, 0.0d0)
    Sigma_U = (0.0d0, 0.0d0)
    Sigma_ph = (0.0d0, 0.0d0)
    Sigma_pb = (0.0d0, 0.0d0)
    Sigma_pp = (0.0d0, 0.0d0)

    Sigma_V = (0.0d0, 0.0d0)
    SigV_U = (0.0d0, 0.0d0)
    SigV_ph = (0.0d0, 0.0d0)
    SigV_pb = (0.0d0, 0.0d0)
    SigV_pp = (0.0d0, 0.0d0)

    !$omp parallel private(map_k, sig_inter, map_q, map_q_sym, map_kq, map_qmk, map_kmq, map_mk, Gkq, Gkmq, Gqmk, symm_list, map_l4, map_l3, idx_l4, idx_l3, map_l3_mnu, idx_l3_mnu, idx_k, sig_U, sig_ph, sig_pb, sig_pp)
    
    !$omp do collapse(2)
    do k_0 = 1, Nf
    do kx = 1, Nx_IBZ
    do ky = 1, kx

      map_k = indxmap(kx, ky, k_0)

      call index_minusF(map_k, map_mk) !for use later

      sig_inter = 0.0d0

      sig_U = 0.0d0
      sig_ph = 0.0d0
      sig_pb = 0.0d0
      sig_pp = 0.0d0

      DO idx_q = 1, Nb
        map_q = Index_Bosonic_IBZ(id * Nb + idx_q) !correct map on any node


        !get all symmetries that are needed to fill whole BZ
        call list_symmetries(map_q, symm_list)                  

        !loop over symmetries
        DO is = 1, Ns

          !leave out the ones that we already took
          IF(.NOT. symm_list(is)) CYCLE
          !create symmetrized q2 to index formfactors
          call symmetry_operation(is, map_q, map_q_sym)

          !calculation of green's functions
  
          call Index_FaddB(map_k, map_q_sym, map_kq)
          Gkq = get_green(map_kq)
  
          call index_FaddB(map_mk, map_q_sym, map_qmk)
          Gqmk = get_green(map_qmk)
      
          call index_minusF(map_qmk, map_kmq)
          Gkmq = get_green(map_kmq)
  
          DO l4 = 1, Nl
            !loop over enlarged frequency box - frange = 0 -> no kernels
            DO l4_0 = -f_range * Nf + 1, (f_range + 1) * Nf

              map_l4 = Indxmap_L(l4, l4_0)

              DO l3 = 1, Nl
                map_l3 = Indxmap_L(l3, map_k%iw)
                idx_l3 = List_Index_L(map_l3)

                !for inside old box just take the old code
                IF((l4_0 .LE. Nf) .AND. (l4_0 .GE. 1)) THEN

                  !calculate index for dereferencing arrays
                  idx_l4 = List_Index_L(map_l4)




                  ! ---

                  !in case of U2 correction subtract U2 term
                  if(map_l4%il == 1 .AND. map_l3%il == 1) then

                    if(use_U2) then
                      !subtract bare U in case of taken 2nd order correction
                      sig_inter = sig_inter - 2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                                  Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
                    endif !use_U2

                  end if

                  !!old biased version
                  !!ph biased
                  !sig_inter = sig_inter + 1.0d0/3.0d0 * F_d_LL(idx_l3, idx_l4, idx_q) * Gkq * &
                  !        FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
                  !sig_inter = sig_inter - 1.0d0/3.0d0 * F_m_LL(idx_l3, idx_l4, idx_q) * Gkq * &
                  !        FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)

                  !sig_inter = sig_inter - 2.0d0/3.0d0 * F_m_LL(idx_l3, idx_l4, idx_q) * Gkq * &
                  !        FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
 
                  !!pp biased
                  !sig_inter = sig_inter + 1.0d0/3.0d0 * F_s_LL(idx_l3, idx_l4, idx_q) * &
                  !        Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp(idx_l4, idx_q)
  
                  !sig_inter = sig_inter + 1.0d0 / 3.0d0 * F_t_LL(idx_l3, idx_l4, idx_q) * &
                  !        Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp(idx_l4, idx_q)

                  ! ---
 
                  !contribution from fully irreducible vertex
                  if(map_l4%il == 1 .AND. map_l3%il == 1) then
                    sig_inter = sig_inter + L_d(map_l3%iw, map_l4%iw, map_q%iw) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                            Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)

                    sig_inter = sig_inter - L_m(map_l3%iw, map_l4%iw, map_q%iw) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                            Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)

                  end if
          
                  !!ph contribution from both ph and \bar{ph}
                  sig_inter = sig_inter + G_d_LL(idx_l3, idx_l4, idx_q) * &
                          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
  
                  sig_inter = sig_inter - 1.0d0 * G_m_LL(idx_l3, idx_l4, idx_q) * &
                          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)

                  sig_inter = sig_inter - 2.0d0 * G_m_LL(idx_l3, idx_l4, idx_q) * &
                          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
  
                  !!pp contribution
  
                  sig_inter = sig_inter + G_s_LL(idx_l3, idx_l4, idx_q) * &
                          Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp(idx_l4, idx_q)
  
                  sig_inter = sig_inter + G_t_LL(idx_l3, idx_l4, idx_q) * &
                          Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp(idx_l4, idx_q)


                  !!split up of contributions

                  !if(ite == ite_max - 1) then

                  !  !contribution from fully irreducible vertex
                  !  if(map_l4%il == 1 .AND. map_l3%il == 1) then
                  !    sig_U = sig_U + L_d(map_l3%iw, map_l4%iw, map_q%iw) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                  !            Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
  
                  !    sig_U = sig_U - L_m(map_l3%iw, map_l4%iw, map_q%iw) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                  !            Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
  
                  !    if(use_U2) then
                  !      !subtract bare U in case of taken 2nd order correction
                  !      sig_U = sig_U - 2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                  !              Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
                  !    endif !use_U2
                  !  end if
            
                  !  !ph contribution from both ph and \bar{ph}
                  !  sig_ph = sig_ph + G_d_LL(idx_l3, idx_l4, idx_q) * &
                  !          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
    
                  !  sig_ph = sig_ph - G_m_LL(idx_l3, idx_l4, idx_q) * &
                  !          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
    
                  !  sig_pb = sig_pb - 2.0d0 * G_m_LL(idx_l3, idx_l4, idx_q) * &
                  !          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph(idx_l4, idx_q)
  
                  !  !pp contribution
                  !  sig_pp = sig_pp + G_s_LL(idx_l3, idx_l4, idx_q) * &
                  !          Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp(idx_l4, idx_q)
    
                  !  sig_pp = sig_pp + G_t_LL(idx_l3, idx_l4, idx_q) * &
                  !          Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp(idx_l4, idx_q)

                  !end if !ite = ite_max - 1 



  
                  !now add TR part
                  IF (map_q%iw == 1) CYCLE
  
                  !take negative nu (via copy-pasted rule)
                  !negation of l3 is taken care of by signs
                  map_l3_mnu = Indxmap_L(map_l3%il, -map_l3%iw + Nf + 1)
                  idx_l3_mnu = List_Index_L(map_l3_mnu)
 
                  ! ---
   

                  !in case of U2 correction subtract U2 term
                  if(map_l4%il == 1 .AND. map_l3%il == 1) then

                    if(use_U2) then
                      !subtract bare U in case of taking 2nd order correction
                      sig_inter = sig_inter - 2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                                  Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
                    endif !use_U2
                  endif
                

                  !!biased version
                  !!ph biased
                  !sig_inter = sig_inter + 1.0d0/3.0d0 * signs(map_l3%il) * CONJG(F_d_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !        Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
                  !sig_inter = sig_inter - 1.0d0/3.0d0 * signs(map_l3%il) * CONJG(F_m_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !        Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))

                  !sig_inter = sig_inter - 2.0d0/3.0d0 * signs(map_l3%il) * CONJG(F_m_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !        Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))

                  !!pp biased
                  !sig_inter = sig_inter + signs(map_l3%il) * &
                  !        1.0d0/3.0d0 * CONJG(F_s_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !        CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp(idx_l4, idx_q))
  
                  !sig_inter = sig_inter + signs(map_l3%il) * &
                  !        1.0d0/3.0d0 * CONJG(F_t_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !        CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp(idx_l4, idx_q))
  
                  ! ---
  
                  !contribution from fully irreducible vertex
                  if(map_l4%il == 1 .AND. map_l3%il == 1) then
                    sig_inter = sig_inter + CONJG(L_d(map_l3_mnu%iw, map_l4%iw, map_q%iw)) * &
                            signs(map_l3%il) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                            Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))

                    sig_inter = sig_inter - CONJG(L_m(map_l3_mnu%iw, map_l4%iw, map_q%iw)) * &
                            signs(map_l3%il) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                            Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))

                  end if
  
                  !!ph contributions
                  sig_inter = sig_inter + signs(map_l3%il) * &
                          CONJG(G_d_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
  
                  sig_inter = sig_inter - 1.0d0 * signs(map_l3%il) * &
                          CONJG(G_m_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
 
                  sig_inter = sig_inter - 2.0d0 * signs(map_l3%il) * &
                          CONJG(G_m_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
  
                  !pp contributions
                  sig_inter = sig_inter + signs(map_l3%il) * &
                          CONJG(G_s_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                          CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp(idx_l4, idx_q))
  
                  sig_inter = sig_inter + signs(map_l3%il) * &
                          CONJG(G_t_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                          CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp(idx_l4, idx_q))


                  !if(ite == ite_max - 1) then

                  !  !contribution from fully irreducible vertex
                  !  if(map_l4%il == 1 .AND. map_l3%il == 1) then
                  !    sig_U = sig_U + CONJG(L_d(map_l3_mnu%iw, map_l4%iw, map_q%iw)) * &
                  !            signs(map_l3%il) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                  !            Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
  
                  !    sig_U = sig_U - CONJG(L_m(map_l3_mnu%iw, map_l4%iw, map_q%iw)) * &
                  !            signs(map_l3%il) * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                  !            Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
  
                  !    if(use_U2) then
                  !      !subtract bare U in case of taking 2nd order correction
                  !      sig_U = sig_U - 2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                  !              Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
                  !    endif !use_U2
                  !  end if
    
                  !  !ph contributions
                  !  sig_ph = sig_ph + signs(map_l3%il) * &
                  !          CONJG(G_d_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
    
                  !  sig_ph = sig_ph - signs(map_l3%il) * &
                  !          CONJG(G_m_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
  
                  !  sig_pb = sig_pb - 2.0d0 * signs(map_l3%il) * &
                  !          CONJG(G_m_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph(idx_l4, idx_q))
   
    
                  !  !pp contributions
                  !  sig_pp = sig_pp + signs(map_l3%il) * &
                  !          CONJG(G_s_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !          CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp(idx_l4, idx_q))
    
                  !  sig_pp = sig_pp + signs(map_l3%il) * &
                  !          CONJG(G_t_LL(idx_l3_mnu, idx_l4, idx_q)) * &
                  !          CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp(idx_l4, idx_q))

                  !end if !ite == ite_max - 1

  
  
                !if l4_0 is not inside the old box
                !use kernels and chi_aux from larger box
                ELSE
  
                  !set index for array dereferencing - usual thinking
                  idx_l4 = (l4 - 1) * (2 * f_range + 1) * Nf + f_range * Nf + l4_0



                  if(.not. use_U2) then
                    !contribution from bare U - later do this more accurately via FT
                    if(map_l4%il == 1 .AND. map_l3%il == 1) then
                      sig_inter = sig_inter + 2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                                  Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph_L(idx_l4, idx_q)

                      if(ite == ite_max - 1) then
                        sig_U = sig_U + 2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                               Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph_L(idx_l4, idx_q)
                      end if !ite == ite_max - 1

                    end if
                  endif !use_U2
 


         
                  !ph contribution from both ph and \bar{ph}
                  !->kernel contributions get added in exactly the same way as
                  !phis - just use the maps
                  sig_inter = sig_inter + kernel( 'd', map_l3, map_l4, map_q ) * &
                          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph_L(idx_l4, idx_q)
 
                  sig_inter = sig_inter - 1.0d0 * kernel( 'm', map_l3, map_l4, map_q ) * &
                          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph_L(idx_l4, idx_q)

                  sig_inter = sig_inter - 2.0d0 * kernel( 'm', map_l3, map_l4, map_q ) * &
                          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph_L(idx_l4, idx_q)
  
                  !pp contribution
                  sig_inter = sig_inter + kernel( 's', map_l3, map_l4, map_q ) * &
                          Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp_L(idx_l4, idx_q)
  
                  sig_inter = sig_inter + kernel( 't', map_l3, map_l4, map_q ) * &
                          Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp_L(idx_l4, idx_q)


                  !if(ite == ite_max - 1) then

                  !  !ph contribution from both ph and \bar{ph}
                  !  sig_ph = sig_ph + kernel( 'd', map_l3, map_l4, map_q ) * &
                  !          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph_L(idx_l4, idx_q)
                  !  !chi_aux_L is chi in larger parameter space 
   
                  !  sig_ph = sig_ph - 1.0d0 * kernel( 'm', map_l3, map_l4, map_q ) * &
                  !          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph_L(idx_l4, idx_q)

                  !  sig_pb = sig_pb - 2.0d0 * kernel( 'm', map_l3, map_l4, map_q ) * &
                  !          Gkq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_ph_L(idx_l4, idx_q)
    
                  !  !pp contribution
                  !  sig_pp = sig_pp + kernel( 's', map_l3, map_l4, map_q ) * &
                  !          Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp_L(idx_l4, idx_q)
    
                  !  sig_pp = sig_pp + kernel( 't', map_l3, map_l4, map_q ) * &
                  !          Gqmk * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * Chi_aux_pp_L(idx_l4, idx_q)
    
                  !end if !ite = ite_max - 1
  
  
  
                  !now add TR part
                  IF (map_q%iw == 1) CYCLE
  
                  !take negative nu (via copy-pasted rule)
                  !negation of l3 is taken care of by signs
                  map_l3_mnu = Indxmap_L(map_l3%il, -map_l3%iw + Nf + 1)
                  idx_l3_mnu = List_Index_L(map_l3_mnu)
 

 
  
                  if(.not. use_U2) then
                    !contribution from bare U - later do this more accurately via FT
                    if(map_l4%il == 1 .AND. map_l3%il == 1) then
                      sig_inter = sig_inter + 2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                                  Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph_L(idx_l4, idx_q))

                      if(ite == ite_max - 1) then
                        sig_U = sig_U + 2.0d0 * xU * 1.0d0/(FF(1, 1, 1) * FF(1, 1, 1)) * &
                                    Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph_L(idx_l4, idx_q))
                      end if !ite == ite_max - 1

                    end if
                  endif !use_U2


                  !ph contributions
                  sig_inter = sig_inter + signs(map_l3%il) * &
                          CONJG(kernel( 'd', map_l3_mnu, map_l4, map_q )) * &
                          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph_L(idx_l4, idx_q))
  
                  sig_inter = sig_inter - 1.0d0 * signs(map_l3%il) * &
                          CONJG(kernel( 'm', map_l3_mnu, map_l4, map_q )) * &
                          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph_L(idx_l4, idx_q))

                  sig_inter = sig_inter - 2.0d0 * signs(map_l3%il) * &
                          CONJG(kernel( 'm', map_l3_mnu, map_l4, map_q )) * &
                          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph_L(idx_l4, idx_q))
 
  
                  !pp contributions
                  sig_inter = sig_inter + signs(map_l3%il) * &
                          CONJG(kernel( 's', map_l3_mnu, map_l4, map_q )) * &
                          CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp_L(idx_l4, idx_q))
  
                  sig_inter = sig_inter + signs(map_l3%il) * &
                          CONJG(kernel( 't', map_l3_mnu, map_l4, map_q )) * &
                          CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp_L(idx_l4, idx_q))


                  !if(ite == ite_max - 1) then

                  !  !ph contributions
                  !  sig_ph = sig_ph + signs(map_l3%il) * &
                  !          CONJG(kernel( 'd', map_l3_mnu, map_l4, map_q )) * &
                  !          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph_L(idx_l4, idx_q))
  
                  !  sig_ph = sig_ph - 1.0d0 * signs(map_l3%il) * &
                  !          CONJG(kernel( 'm', map_l3_mnu, map_l4, map_q )) * &
                  !          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph_L(idx_l4, idx_q))
 
                  !  sig_pb = sig_pb - 2.0d0 * signs(map_l3%il) * &
                  !          CONJG(kernel( 'm', map_l3_mnu, map_l4, map_q )) * &
                  !          Gkmq * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_ph_L(idx_l4, idx_q))
  
                  !  !pp contributions
                  !  sig_pp = sig_pp + signs(map_l3%il) * &
                  !          CONJG(kernel( 's', map_l3_mnu, map_l4, map_q )) * &
                  !          CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp_L(idx_l4, idx_q))
  
                  !  sig_pp = sig_pp + signs(map_l3%il) * &
                  !          CONJG(kernel( 't', map_l3_mnu, map_l4, map_q )) * &
                  !          CONJG(Gkq) * FF_inv(is, map_k%ix, map_k%iy, map_l3%il) * CONJG(Chi_aux_pp_L(idx_l4, idx_q))

                  !end if !ite == ite_max - 1
  
                !l4_0 inside box if
                END IF
              !l3
              END DO
            !l4_0 <- frequency that lives in enlarged box
            END DO
          !l4
          END DO
        !is
        END DO
      !q
      END DO

      !index corresponding to k-index -> only needed at this place
      idx_k = ( (((kx - 1) * kx))/2 + ky - 1 ) * Nf + k_0

      !add intermediate result to to be reduced sigma_V
      Sigma_V(idx_k) = Sigma_V(idx_k) - xU/beta/beta/Nc/Nc/2.0d0 * sig_inter

      !split up contributions
      SigV_U(idx_k) =  SigV_U(idx_k) -  xU/beta/beta/Nc/Nc/2.0d0 * sig_U
      SigV_ph(idx_k) = SigV_ph(idx_k) - xU/beta/beta/Nc/Nc/2.0d0 * sig_ph
      SigV_pb(idx_k) = SigV_pb(idx_k) - xU/beta/beta/Nc/Nc/2.0d0 * sig_pb
      SigV_pp(idx_k) = SigV_pp(idx_k) - xU/beta/beta/Nc/Nc/2.0d0 * sig_pp

    

    !outer fermionic loop - k-argument
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel


    !sum up sigmas from all tasks
    IF (.NOT. ALLOCATED(Sigma_Reduced)) ALLOCATE (Sigma_Reduced(Nred * 2))

    IF (.NOT. ALLOCATED(SigRed_U)) ALLOCATE (SigRed_U(Nred * 2))
    IF (.NOT. ALLOCATED(SigRed_ph)) ALLOCATE (SigRed_ph(Nred * 2))
    IF (.NOT. ALLOCATED(SigRed_pb)) ALLOCATE (SigRed_pb(Nred * 2))
    IF (.NOT. ALLOCATED(SigRed_pp)) ALLOCATE (SigRed_pp(Nred * 2))

    Sigma_Reduced = (0.0d0, 0.0d0)
    SigRed_U = (0.0d0, 0.0d0)
    SigRed_ph = (0.0d0, 0.0d0)
    SigRed_pb = (0.0d0, 0.0d0)
    SigRed_pp = (0.0d0, 0.0d0)

    !reduce only IBZ part
    call MPI_AllReduce(Sigma_V, Sigma_Reduced, Nred * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, rc)

    if(ite == ite_max - 1) then
      !split up contributions
      call MPI_AllReduce(SigV_U,  SigRed_U, Nred * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_AllReduce(SigV_ph, SigRed_ph, Nred * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_AllReduce(SigV_pb, SigRed_pb, Nred * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_AllReduce(SigV_pp, SigRed_pp, Nred * 2, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, rc)
    end if !ite == ite_max - 1

    !now symmetrize sigma - calculate rest of BZ
    call symmetrize_array(Nf, Sigma_Reduced, Sigma) 

    if(ite == ite_max - 1) then
      !split up contributions 
      call symmetrize_array(Nf, SigRed_U, Sigma_U)  
      call symmetrize_array(Nf, SigRed_ph, Sigma_ph)  
      call symmetrize_array(Nf, SigRed_pb, Sigma_pb)  
      call symmetrize_array(Nf, SigRed_pp, Sigma_pp)  
    end if !ite == ite_max - 1

    if(use_U2) then
    
      !caculate 2nd order diagram more accurately via FFT
      call selfenergy_2nd(Grt, order_2nd)
  
      !add 2nd order explicitly ...
      Sigma = Sigma + order_2nd
      ! ... also for split up selfenergy
      if(ite == ite_max - 1) Sigma_U = Sigma_U + order_2nd

    endif !use_U2

    Deallocate (Sigma_Reduced)

    !now assign to local part the DMFT selfenergy

    !if(use_sig_dmft) then

    !  do nu = 1, Nf

    !    !calculate local part of selfenergy for given frequency ...
    !    sig_loc = 0.0d0
    !    do ky = 1, Ny
    !      do kx = 1, Nx
    !        sig_loc = sig_loc + Sigma(((kx - 1) * Ny + ky - 1) * Nf + nu) / (Nx * Ny)
    !      end do !kx
    !    end do !ky

    !    ! ... and replace it by the DMFT selfenergy
    !    do ky = 1, Ny
    !      do kx = 1, Nx
    !        Sigma(((kx - 1) * Ny + ky - 1) * Nf + nu) = Sigma(((kx - 1) * Ny + ky - 1) * Nf + nu) - &
    !        sig_loc + Sigma_DMFT(nu)
    !      end do !kx
    !    end do !ky

    !  end do !nu

    !end if !use_dmft_sig

  end subroutine self_energy

  ! ---
  ! function to calculate the 2nd order contribution to the schwinger-dyson
  ! equation
  ! ---
  subroutine selfenergy_2nd(Grt, order_2nd)

    !real space & time green's function
    complex(dp), intent(in) :: Grt(Nx * Ngrain, Ny * Ngrain, Nf)
    !return value
    complex(dp), dimension(Nx * Ny * Nf), intent(out) :: order_2nd

    !loop variables
    integer :: kx, ky, kxm, kym, kx_c, ky_c, iTau, nu
    
    !index
    integer :: idx_k, idx_mk

    !holds product of 3 green's functions
    complex(dp), allocatable :: GGG(:, :, :)

    !for fitting
    real(dp) :: FD1, FD2

    !string needed for determination of operation - probably
    character(len=30) :: Mtype

    !temporary to hold fourier-transfor to w-space
    complex(dp), allocatable :: ret_temp(:)

    order_2nd = 0.0d0

    if (.NOT. allocated(GGG)) allocate(GGG(Nx * Ngrain, Ny * Ngrain, Nf))
    if (.not. allocated(ret_temp)) allocate(ret_temp(Nf/2))
    ret_temp = 0.0d0
    GGG = 0.0d0
    MType = 'Fermionic'

    do iTau = 1, Nf
      do ky = 1, Ny * Ngrain
        kym = merge(1, Ny * Ngrain - ky + 2, ky == 1)
        do kx = 1, Nx * Ngrain
          kxm = merge(1, Nx * Ngrain - kx + 2, kx == 1)
          !idea : convolution becomes product in FT-space
          GGG(kx, ky, iTau) = Grt(kx, ky, iTau)**2 * Grt(kxm, kym, Nf - iTau + 1) * xU * xU
        end do !kx
      end do !ky
    
      ! now change it to k-t space by using FFT - first momentum arguments
      call fftb2d(Nx * Ngrain, Ny * Ngrain, GGG(1:Nx * Ngrain, 1:Ny * Ngrain, iTau), C_wave_x, C_wave_y)
    end do !iTau
  
    do kx = 1, Nx
      !argument in finer grid
      kx_c = (kx - 1) * Ngrain + 1
      do ky = 1, Ny
        !argument in finer grid
        ky_c = (ky - 1) * Ngrain + 1

        !fit the time component with polynomials
        call FDfit(Nf, dble(GGG(kx_c, ky_c, 1:Nf)), beta/(Nf - 1), FD1, FD2)

        !perform trafo
        call nfourier(MType, Nf - 1, Nf/2, FD1, FD2, dble(GGG(kx_c, ky_c, 1:Nf)), ret_temp)

        !actually fill output array
        do nu = 1, Nf/2
          idx_k  = ((kx - 1) * Ny + ky - 1) * Nf + nu
          idx_mk = ((kx - 1) * Ny + ky - 1) * Nf + Nf - nu + 1
          order_2nd(idx_k)  = order_2nd(idx_k) + conjg(ret_temp(Nf/2 + 1 - nu))
          order_2nd(idx_mk) = order_2nd(idx_mk) + ret_temp(Nf/2 + 1 - nu)
        end do !nu

      end do !ky
    end do !kx

    deallocate(GGG)
    deallocate(ret_temp)

  end subroutine selfenergy_2nd



    !-------------------------------------------------------------------

  subroutine calc_chi_aux(Chi_aux_ph, Chi_aux_pp, Chi_aux_ph_L, Chi_aux_pp_L)

    !chis in small box
    complex(dp), intent(out) :: Chi_aux_ph(Nz, Nb)
    complex(dp), intent(out) :: Chi_aux_pp(Nz, Nb)
    !chis in enlarged box
    complex(dp), intent(out) :: Chi_aux_ph_L((2 * f_range + 1) * Nf * Nl, Nb)
    complex(dp), intent(out) :: Chi_aux_pp_L((2 * f_range + 1) * Nf * Nl, Nb)

    !indices in loop
    integer :: idx_q, idx_l4
    !k-momenta - have same frequency as l4
    integer :: kx, ky
    !index maps
    type(Indxmap) map_q, map_k, map_mk, map_kpq, map_qmk
    type(Indxmap_L) map_l4
    !Green's function for a particular argument
    complex(dp) Gk, Gqmk, Gkpq

    !variables for extended frequency loop
    integer :: l4_0, l4 !freq and l-arg

    !loop variables for coarse graining
    integer :: kx_grain, ky_grain

    Chi_aux_ph = 0.0d0
    Chi_aux_pp = 0.0d0
    Chi_aux_ph_L = 0.0d0
    Chi_aux_pp_L = 0.0d0

    !$omp parallel private(map_q, map_l4, idx_l4, map_k, map_kpq, map_mk, map_qmk, Gk, Gkpq, Gqmk)

    !$omp do collapse(3)
    DO idx_q = 1, Nb
      !loop over larger frequency box
      DO l4 = 1, Nl
        DO l4_0 = -f_range * Nf + 1, (f_range + 1) * Nf

          map_q = Index_Bosonic_IBZ(id * Nb + idx_q) !correct map on any node
          map_l4 = Indxmap_L(l4, l4_0)
    
          DO kx = 1, Nx
            DO ky = 1, Ny

              IF((l4_0 .LE. Nf) .AND. (l4_0 .GE. 1)) THEN
                !set index for array dereferencing
                idx_l4 = list_index_L(map_l4)

                map_k = Indxmap(kx, ky, map_l4%iw)
                call Index_FaddB(map_k, map_q, map_kpq)
  
                !caculate chi_aux_ph
  
                !Gk = get_green(map_k)
                !Gkpq = get_green(map_kpq)
      
                !Chi_aux_ph(idx_l4, idx_q) = Chi_aux_ph(idx_l4, idx_q) + &
                !                            Gk * Gkpq * FF(map_k%ix, map_k%iy, map_l4%il)  

                do kx_grain = 1, Ngrain
                  do ky_grain = 1, Ngrain
                
                    Gk = get_green_coarse(kx_grain, ky_grain, map_k)
                    Gkpq = get_green_coarse(kx_grain, ky_grain, map_kpq)
      
                    Chi_aux_ph(idx_l4, idx_q) = Chi_aux_ph(idx_l4, idx_q) + &
                                                Gk * Gkpq * &
                                                FF(map_k%ix, map_k%iy, map_l4%il) * &
                                                1.0d0/Ngrain/Ngrain
  
                  end do !kx_grain
                end do!ky_grain

                !calculate chi_aux_pp
  
                call index_minusF(map_k, map_mk) 
                call index_FaddB(map_mk, map_q, map_qmk) 
                
                !Gqmk = get_green(map_qmk)
  
                !Chi_aux_pp(idx_l4, idx_q) = Chi_aux_pp(idx_l4, idx_q) + &
                !                            Gk * Gqmk * FF(map_k%ix, map_k%iy, map_l4%il)

                do kx_grain = 1, Ngrain
                  do ky_grain = 1, Ngrain
                
                    Gk = get_green_coarse(kx_grain, ky_grain, map_k)
                    Gqmk = get_green_coarse(Ngrain - kx_grain + 1, Ngrain - ky_grain + 1, map_qmk)
      
                    Chi_aux_pp(idx_l4, idx_q) = Chi_aux_pp(idx_l4, idx_q) + &
                                                Gk * Gqmk * &
                                                FF(map_k%ix, map_k%iy, map_l4%il) * &
                                                1.0d0/(Ngrain * Ngrain)
  
                  end do !kx_grain
                end do!ky_grain

              
              END IF !l4_0 inside box

                !set index for array dereferencing
                ! addition at the end ensures that one starts at 1
                idx_l4 = (l4 - 1) * (2 * f_range + 1) * Nf + f_range * Nf + l4_0

                map_k = Indxmap(kx, ky, map_l4%iw)
                call Index_FaddB(map_k, map_q, map_kpq)
!
                !caculate chi_aux_ph

                !Gk = get_green(map_k)
                !Gkpq = get_green(map_kpq)
      
                !Chi_aux_ph_L(idx_l4, idx_q) = Chi_aux_ph_L(idx_l4, idx_q) + &
                !                              Gk * Gkpq * FF(map_k%ix, map_k%iy, map_l4%il)  

                do kx_grain = 1, Ngrain
                  do ky_grain = 1, Ngrain
                
                    Gk = get_green_coarse(kx_grain, ky_grain, map_k)
                    Gkpq = get_green_coarse(kx_grain, ky_grain, map_kpq)
      
                    Chi_aux_ph_L(idx_l4, idx_q) = Chi_aux_ph_L(idx_l4, idx_q) + &
                                                  Gk * Gkpq * &
                                                  FF(map_k%ix, map_k%iy, map_l4%il) * &
                                                  1.0d0/(Ngrain * Ngrain)
  
                  end do !kx_grain
                end do!ky_grain

                !calculate chi_aux_pp
  
                call index_minusF(map_k, map_mk) 
                call index_FaddB(map_mk, map_q, map_qmk) 
                
                !Gqmk = get_green(map_qmk)
  
                !Chi_aux_pp_L(idx_l4, idx_q) = Chi_aux_pp_L(idx_l4, idx_q) + &
                !                              Gk * Gqmk * FF(map_k%ix, map_k%iy, map_l4%il)

                do kx_grain = 1, Ngrain
                  do ky_grain = 1, Ngrain
                
                    Gk = get_green_coarse(kx_grain, ky_grain, map_k)
                    Gqmk = get_green_coarse(Ngrain - kx_grain + 1, Ngrain - ky_grain + 1, map_qmk)
      
                    Chi_aux_pp_L(idx_l4, idx_q) = Chi_aux_pp_L(idx_l4, idx_q) + &
                                                  Gk * Gqmk * &
                                                  FF(map_k%ix, map_k%iy, map_l4%il) * &
                                                  1.0d0/(Ngrain * Ngrain)
  
                  end do !kx_grain
                end do!ky_grain

            END DO!kx
          END DO!ky
        END DO !l4_0
      END DO!l4
    END DO!q
    !$omp end do

    !$omp end parallel

  end subroutine calc_chi_aux



end module parquet_selfenergy
