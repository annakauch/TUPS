module parquet_util

! .. this is a version with form factors started as a new branch of the victory code

  use global_parameter
  use parquet_ini
  use parquet_formfactors

  implicit none

contains

  !------------------------------------------------------------------------------
  subroutine Initialization

    implicit none

    !variables for fine-grain dispersion
    integer   :: ix, ix2, iy, iy2
    !hold actual momentum value in dispersion setup and smallest possible momentum
    real(dp)  :: px, py, dkx, dky
!--------------------------------------

    !for convienience at readin
    integer :: temp

    !array indices
    integer :: idx, idx_k
    !momenta in loops
    integer :: kx, ky, kx1, ky1
    !frequencies in loops
    integer :: nu
    !store data from readin
    real(dp) :: dat1, dat2, dat3
    integer   :: l

    !variables for reading data from IBZ
    integer :: idx_IBZ, idx_k1, idx_k2

    !readin data
    real(dp) dat_dr, dat_di, dat_mr, dat_mi, dat_sr, dat_si, dat_tr, dat_ti

    !just some strings for names
    integer :: input_status

    !read in of Lambdas
    !value of read in lambda - real & imag
    real(dp) :: L1_r, L1_i, L2_r, L2_i
    real(dp) :: itemp
    !frequencies to loop over
    integer :: k1_0, k2_0, q_0
    !freq. arguments in readin of Lambdas
    !integer :: f_q, f_k1, f_k2
    !loop through bosonic arguments
    !integer :: idx_q
    !map for bosonic argument
    type(Indxmap) :: map_k


    !auxillary vertices are used during the readin process
    complex(dp), allocatable :: F_d_aux(:, :, :) ! density channel
    complex(dp), allocatable :: F_m_aux(:, :, :) ! magnetic channel
    complex(dp), allocatable :: F_s_aux(:, :, :) ! singlet channel
    complex(dp), allocatable :: F_t_aux(:, :, :) ! triplet channel

    ! irreducible vertex in each channel, :
    complex(dp), allocatable :: G_d_aux(:, :, :) ! density channel
    complex(dp), allocatable :: G_m_aux(:, :, :) ! magnetic channel
    complex(dp), allocatable :: G_s_aux(:, :, :) ! singlet channel
    complex(dp), allocatable :: G_t_aux(:, :, :) ! triplet channel

    !2P irreducible vertex 

    complex(dp), allocatable :: L_d_aux(:, :, :) ! density channel
    complex(dp), allocatable :: L_m_aux(:, :, :) ! magnetic channel
    complex(dp), allocatable :: L_s_aux(:, :, :) ! singlet channel
    complex(dp), allocatable :: L_t_aux(:, :, :) ! triplet channel


    ! --- initialize variables needed for formfactor array ---
    !for instance Nl

    call Init_Formfactor_variables

    ! --- initialize the character table - atm for the square lattice

    call init_char_table(length_char_table, Ns)

    ! --- set up array to be able to convert idices to index_maps ---
    if (.NOT. allocated(Index_fermionic)) allocate (Index_fermionic(Nt))
    if (.NOT. allocated(Index_bosonic)) allocate (Index_bosonic(Nt/2))
    if (.NOT. allocated(Index_bosonic_IBZ)) allocate (Index_bosonic_IBZ(Nred))
    !additional index_map for L_indices
    if (.NOT. allocated(Index_L)) allocate (Index_L(Nz))

    do kx = 1, Nx
      do ky = 1, Ny
        do nu = 1, Nf
          idx = ((kx - 1)*Ny + ky - 1)*Nf + nu
          Index_fermionic(idx) = indxmap(kx, ky, nu)
        end do
        do nu = 1, Nf/2
          idx = ((kx - 1)*Ny + ky - 1)*Nf/2 + nu
          Index_bosonic(idx) = indxmap(kx, ky, nu)
        end do
      end do
    end do

    ! Lattice dependent
    do kx = 1, Nx_IBZ
      do ky = 1, kx
        do nu = 1, Nf/2
          idx = (((kx - 1)*kx)/2 + ky - 1)*Nf/2 + nu ! Lattice dependent
          Index_bosonic_IBZ(idx) = indxmap(kx, ky, nu)
        end do
      end do
    end do
    ! create Index for l-space
    DO l = 1, Nl
      DO nu = 1, Nf
        idx = (l - 1)*Nf + nu
        Index_L(idx) = indxmap_L(l, nu)
      END DO
    END DO

    ! --- set up the tight-binding dispersion ---

    if (.NOT. allocated(Ek)) allocate (Ek(Nx, Ny))
    dkx = Two*Pi/Nx
    dky = Two*Pi/Ny
    Ek = Zero
    do kx = 1, Nx
      px = dkx*(kx - 1)
      do ky = 1, Ny
        py = dky*(ky - 1)
        if (Nx > 1) Ek(kx, ky) = Ek(kx, ky) - 2.0d0 * cos(px)
        if (Ny > 1) Ek(kx, ky) = Ek(kx, ky) - 2.0d0 * cos(py)
      end do
    end do
    if (Nc == 1) Ek = 0.0d0

    ! --- set up the tight-binding dispersion for the enlarged case ---
    if (.NOT. allocated(Ek_grain)) allocate (Ek_grain(Ngrain, Ngrain, Nx, Ny))
    dkx = 2.0d0 * Pi/Nx
    dky = 2.0d0 * Pi/Ny
    ! --- dispersion for fine-graining
    Ek_grain = 0.0d0
    do iy = 1, Ny
      do ix = 1, Nx
        do iy2 = 1, Ngrain
          do ix2 = 1, Ngrain

            px = dkx * (ix - 1) + (ix2 - (Ngrain + 1)/2)*dkx/dble(Ngrain)
            py = dky * (iy - 1) + (iy2 - (Ngrain + 1)/2)*dky/dble(Ngrain)

            Ek_grain(ix2, iy2, ix, iy) = Ek_grain(ix2, iy2, ix, iy) - &
                                         2.0d0 * cos(px) - 2.0d0 * cos(py)

          end do
        end do
      end do
    end do
    !case of only one site (like AIM)
    if (Nc == 1) Ek_grain = 0.0d0


    ! --- initialization of Sigma ---

    !Initialize sigma
    IF (.NOT. allocated(Sigma)) allocate (Sigma(Nt))

    IF (.NOT. allocated(Sigma_U)) allocate (Sigma_U(Nt))
    IF (.NOT. allocated(Sigma_ph)) allocate (Sigma_ph(Nt))
    IF (.NOT. allocated(Sigma_pb)) allocate (Sigma_pb(Nt))
    IF (.NOT. allocated(Sigma_pp)) allocate (Sigma_pp(Nt))

    !read in Sigma if corresponding flag is set
    IF (Sigma_ini) THEN

      input_status = 0

      open (unit=2, file='data/inidata/Sigma.init', status='old')

      DO kx1 = 1, Nx
        DO ky1 = 1, Ny
          DO nu = 1, Nf

            idx = (Ny*(kx1 - 1) + ky1 - 1)*Nf + nu
            !temp temp is two momenta, t1 is frequency and t2, t3 the actual complex value
            read (2, *, IOSTAT=input_status) temp, temp, dat1, dat2, dat3
            Sigma(idx) = dcmplx(dat2, dat3)

          END DO
        END DO
      END DO

      IF (input_status > 0) write (*, *) '----------ERROR WHILE READING FILE FOR SIGMA!!!-------------'

      close (2)

    ELSE

      Sigma = (0.0d0, 0.0d0)

    END IF !Sigma_ini

    if (.NOT. allocated(Sigma_H)) allocate (Sigma_H(Nx, Ny))
    Sigma_H = 0.0d0


    ! --- read in of DMFT Sigma ---
    if(use_sig_dmft) then
    
      if (.not. allocated(Sigma_DMFT)) allocate (Sigma_DMFT(Nf))

      open (unit=3, file=fname_sig_dmft, status='old')

        do nu = 1, Nf
          read (3, *) dat1, dat2, dat3
          !Sigma_DMFT(nu) = dcmplx(dat1, dat2)
          Sigma_DMFT(nu) = dcmplx(0.0d0, dat3)
        end do

      close (3)

      !initialize selfenergy with DMFT selfenergy
      do idx_k = 1, Nt
        map_k = Index_fermionic(idx_k)
        Sigma(idx_k) = Sigma_DMFT(map_k%iw)
      end do !idx_k


    end if



    !initalization for fourier transforms

    !if (.NOT. allocated(C_wave_x)) allocate (C_wave_x(4*Nx + 15))
    !if (.NOT. allocated(C_wave_y)) allocate (C_wave_y(4*Ny + 15))

    !call Zffti(Nx, C_wave_x)
    !call Zffti(Ny, C_wave_y)

    if (.NOT. allocated(C_wave_x)) allocate(C_wave_x(4*Nx*Ngrain+15))
    if (.NOT. allocated(C_wave_y)) allocate(C_wave_y(4*Ny*Ngrain+15))   
    
    call Zffti(Nx*Ngrain, C_wave_x)
    call Zffti(Ny*Ngrain, C_wave_y)


    ! --- initialize vertices ---

    ! --- read in local gammas ---

    ! allocate arrays for the local fully irreducible vertex in each node
    ! this is in k-space (ie. constant)

    !if (.NOT. allocated(G_d_loc)) allocate (G_d_loc(Nf, Nf, Nf/2))
    !if (.NOT. allocated(G_m_loc)) allocate (G_m_loc(Nf, Nf, Nf/2))
    !if (.NOT. allocated(G_s_loc)) allocate (G_s_loc(Nf, Nf, Nf/2))
    !if (.NOT. allocated(G_t_loc)) allocate (G_t_loc(Nf, Nf, Nf/2))

    !G_d_loc = 0.0    
    !G_m_loc = 0.0    
    !G_s_loc = 0.0    
    !G_t_loc = 0.0    

    !!read in local gammas from file
    !open(unit=1, file=TRIM(fname_gamloc_ph), status='old')
    !do k1_0 = 1, Nf
    !  do k2_0 = 1, Nf
    !    do q_0 = 1, Nf/2
    !      read(1, *) itemp, itemp, itemp, L1_r, L1_i, L2_r, L2_i
    !      G_d_loc(k1_0, k2_0, q_0) = dcmplx(L1_r, L1_i)
    !      G_m_loc(k1_0, k2_0, q_0) = dcmplx(L2_r, L2_i)
    !    end do
    !  end do
    !end do
    !close(1)


    !!read in local gammas from file
    !open(unit=2, file=TRIM(fname_gamloc_pp), status='old')
    !do k1_0 = 1, Nf
    !  do k2_0 = 1, Nf
    !    do q_0 = 1, Nf/2
    !      read(2, *) itemp, itemp, itemp, L1_r, L1_i, L2_r, L2_i
    !      G_s_loc(k1_0, k2_0, q_0) = dcmplx(L1_r, L1_i)
    !      G_t_loc(k1_0, k2_0, q_0) = dcmplx(L2_r, L2_i)
    !    end do
    !  end do
    !end do
    !close(2)


    ! --- finished readin of local gammas ---


    ! allocate arrays for the local fully irreducible vertex in each node
    ! this is in k-space (ie. constant)
    if (.NOT. allocated(L_d)) allocate (L_d(Nf, Nf, Nf/2))
    if (.NOT. allocated(L_m)) allocate (L_m(Nf, Nf, Nf/2))
    if (.NOT. allocated(L_s)) allocate (L_s(Nf, Nf, Nf/2))
    if (.NOT. allocated(L_t)) allocate (L_t(Nf, Nf, Nf/2))

    ! -- initialize Lambda = 2P irreducible vertex    

    !if(DGA == .false.) then

      L_d = xU
      L_m = -xU
      L_s = 2.0d0 * xU
      L_t = 0.0d0

    !else !if calculating in DGA
    if(DGA) then

      ! - read in input files -

        open(unit=1, file=TRIM(fname_L_ph), status='old')
        do k1_0 = 1, Nf
          do k2_0 = 1, Nf
            do q_0 = 1, Nf/2
              read(1, *) itemp, itemp, itemp, L1_r, L1_i, L2_r, L2_i
              L_d(k1_0, k2_0, q_0) = dcmplx(L1_r, L1_i)
              L_m(k1_0, k2_0, q_0) = dcmplx(L2_r, L2_i)
            end do
          end do
        end do
        close(1)
        
        open(unit=1, file=TRIM(fname_L_pp), status='old')
        do k1_0 = 1, Nf
          do k2_0 = 1, Nf
            do q_0 = 1, Nf/2
              read(1, *) itemp, itemp, itemp, L1_r, L1_i, L2_r, L2_i
              L_s(k1_0, k2_0, q_0) = dcmplx(L1_r, L1_i)
              L_t(k1_0, k2_0, q_0) = dcmplx(L2_r, L2_i)
            end do
          end do
        end do
        close(1)


    end if !DGA

    !allocation of actual vertices
    IF (.NOT. allocated(F_m_LL)) allocate (F_m_LL(Nz, Nz, Nb))
    IF (.NOT. allocated(F_d_LL)) allocate (F_d_LL(Nz, Nz, Nb))
    IF (.NOT. allocated(F_s_LL)) allocate (F_s_LL(Nz, Nz, Nb))
    IF (.NOT. allocated(F_t_LL)) allocate (F_t_LL(Nz, Nz, Nb))

    IF (.NOT. allocated(G_m_LL)) allocate (G_m_LL(Nz, Nz, Nb))
    IF (.NOT. allocated(G_d_LL)) allocate (G_d_LL(Nz, Nz, Nb))
    IF (.NOT. allocated(G_s_LL)) allocate (G_s_LL(Nz, Nz, Nb))
    IF (.NOT. allocated(G_t_LL)) allocate (G_t_LL(Nz, Nz, Nb))

    !FF's needed for initailization so done later
    G_d_LL = (0.0d0, 0.0d0)
    G_m_LL = (0.0d0, 0.0d0)
    G_s_LL = (0.0d0, 0.0d0)
    G_t_LL = (0.0d0, 0.0d0)

    F_d_LL = (0.0d0, 0.0d0)
    F_m_LL = (0.0d0, 0.0d0)
    F_s_LL = (0.0d0, 0.0d0)
    F_t_LL = (0.0d0, 0.0d0)

    if (read_input) then

      !only needed at this point when one reads in vertices
      IF (.NOT. allocated(F_m)) allocate (F_m(Nt, Nt, Nb))
      IF (.NOT. allocated(F_d)) allocate (F_d(Nt, Nt, Nb))
      IF (.NOT. allocated(F_s)) allocate (F_s(Nt, Nt, Nb))
      IF (.NOT. allocated(F_t)) allocate (F_t(Nt, Nt, Nb))

      IF (.NOT. allocated(G_m)) allocate (G_m(Nt, Nt, Nb))
      IF (.NOT. allocated(G_d)) allocate (G_d(Nt, Nt, Nb))
      IF (.NOT. allocated(G_s)) allocate (G_s(Nt, Nt, Nb))
      IF (.NOT. allocated(G_t)) allocate (G_t(Nt, Nt, Nb))

      !This is only for testing; in production code we cannot read in the whole Vertex on each node!

      IF (.NOT. allocated(F_m_aux)) allocate (F_m_aux(Nt, Nt, Nt/2))
      IF (.NOT. allocated(F_d_aux)) allocate (F_d_aux(Nt, Nt, Nt/2))
      IF (.NOT. allocated(F_s_aux)) allocate (F_s_aux(Nt, Nt, Nt/2))
      IF (.NOT. allocated(F_t_aux)) allocate (F_t_aux(Nt, Nt, Nt/2))

      IF (.NOT. allocated(G_m_aux)) allocate (G_m_aux(Nt, Nt, Nt/2))
      IF (.NOT. allocated(G_d_aux)) allocate (G_d_aux(Nt, Nt, Nt/2))
      IF (.NOT. allocated(G_s_aux)) allocate (G_s_aux(Nt, Nt, Nt/2))
      IF (.NOT. allocated(G_t_aux)) allocate (G_t_aux(Nt, Nt, Nt/2))
      !----------------------------------
      !vertices are read in from IBZ

      open (unit=2, file='data/inidata/F_IBZ.init', status='old')
      open (unit=3, file='data/inidata/G_IBZ.init', status='old')

      !read F data into IBZ
      !status = 0 -> everything fine
      !status = -1 -> file ended
      !status = 1 -> some error occured

      input_status = 0

      DO WHILE (.TRUE.)
        !ordering of idx variables determined by format of file
        read (2, *, IOSTAT=input_status) idx_k2, idx_k1, idx_IBZ, &
          dat_dr, dat_di, dat_mr, dat_mi, dat_sr, dat_si, dat_tr, dat_ti

        IF (input_status .NE. 0) EXIT

        !convert IBZ-index to BZ-index
        !map_q = Index_Bosonic_IBZ(idx_IBZ)
        !idx_q_tot = List_Index_B(map_q)

        F_d_aux(idx_k1, idx_k2, idx_IBZ) = dcmplx(dat_dr, dat_di)
        F_m_aux(idx_k1, idx_k2, idx_IBZ) = dcmplx(dat_mr, dat_mi)
        F_s_aux(idx_k1, idx_k2, idx_IBZ) = dcmplx(dat_sr, dat_si)
        F_t_aux(idx_k1, idx_k2, idx_IBZ) = dcmplx(dat_tr, dat_ti)

      END DO

      IF (input_status > 0) write (*, *) '----------ERROR WHILE READING FILE !!!-------------'

      input_status = 0

      !read Gamma data into IBZ
      DO WHILE (.TRUE.)

        read (3, *, IOSTAT=input_status) idx_k2, idx_k1, idx_IBZ, &
          dat_dr, dat_di, dat_mr, dat_mi, dat_sr, dat_si, dat_tr, dat_ti

        IF (input_status .NE. 0) EXIT

        !convert IBZ-index to BZ-index
        !map_q = Index_Bosonic_IBZ(idx_IBZ)
        !idx_q_tot = List_Index_B(map_q)

        G_d_aux(idx_k1, idx_k2, idx_IBZ) = dcmplx(dat_dr, dat_di)
        G_m_aux(idx_k1, idx_k2, idx_IBZ) = dcmplx(dat_mr, dat_mi)
        G_s_aux(idx_k1, idx_k2, idx_IBZ) = dcmplx(dat_sr, dat_si)
        G_t_aux(idx_k1, idx_k2, idx_IBZ) = dcmplx(dat_tr, dat_ti)

      END DO

      IF (input_status > 0) write (*, *) '----------ERROR WHILE READING FILE !!!-------------'

      close (2)
      close (3)

      !now assign only correct bosonic part on each node to vertices
      F_d(:, :, 1:Nb) = F_d_aux(:, :, (id*Nb + 1):(id*Nb + Nb))
      F_m(:, :, 1:Nb) = F_m_aux(:, :, (id*Nb + 1):(id*Nb + Nb))
      F_s(:, :, 1:Nb) = F_s_aux(:, :, (id*Nb + 1):(id*Nb + Nb))
      F_t(:, :, 1:Nb) = F_t_aux(:, :, (id*Nb + 1):(id*Nb + Nb))

      G_d(:, :, 1:Nb) = G_d_aux(:, :, (id*Nb + 1):(id*Nb + Nb))
      G_m(:, :, 1:Nb) = G_m_aux(:, :, (id*Nb + 1):(id*Nb + Nb))
      G_s(:, :, 1:Nb) = G_s_aux(:, :, (id*Nb + 1):(id*Nb + Nb))
      G_t(:, :, 1:Nb) = G_t_aux(:, :, (id*Nb + 1):(id*Nb + Nb))

      deallocate (F_d_aux)
      deallocate (F_m_aux)
      deallocate (F_s_aux)
      deallocate (F_t_aux)

      deallocate (G_d_aux)
      deallocate (G_m_aux)
      deallocate (G_s_aux)
      deallocate (G_t_aux)

      !readin_input if statement
    end if


  end subroutine Initialization

!------------------------------------------------------------------------------

  integer pure function list_index_F(P) result(idx)

    type(Indxmap), intent(in) :: P

    idx = ((P%ix - 1)*Ny + P%iy - 1)*Nf + P%iw

  end function list_index_F

!-------------------------------------------------------------

  integer pure function list_index_L(P) result(idx)

    type(Indxmap_L), intent(in) :: P

    idx = (P%il - 1)*Nf + P%iw

  end function list_index_L

!-------------------------------------------------------------------------------

!this is for the first argument of the susceptibility
  integer pure function list_index_X(P) result(idx)

    type(Indxmap), intent(in) :: P

    idx = ((P%ix - 1)*Nl + (P%iy - 1))*Nf + P%iw

  end function list_index_X

  !------------------------------------------------------------------------------
  integer pure function list_index_IBZ(P) result(idx)

    type(Indxmap), intent(in) :: P

    ! !(only bosonic is now in the IBZ)

    idx = (((P%ix - 1) * P%ix)/2 + P%iy - 1) * Nf/2 + P%iw

  end function list_index_IBZ

  !------------------------------------------------------------------------------
  integer pure function list_index_B(P) result(idx)

    type(Indxmap), intent(in) :: P

    idx = ((P%ix - 1)*Ny + P%iy - 1)*Nf/2 + P%iw

  end function list_index_B

  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------

!transform vertices from k-space into l-space
!this is included into a seperate function as it takes some time and one doesnt always want to do it
!mainly usefull for testing purposes

  subroutine trafo_2D

    integer :: kx1, ky1, kx2, ky2, l1, l2
    integer :: w1, w2
    integer :: idxk1, idxk2
    integer :: idxl1, idxl2
    integer :: idx_q

    F_d_LL = (0.d0, 0.d0)
    F_m_LL = (0.d0, 0.d0)
    F_s_LL = (0.d0, 0.d0)
    F_t_LL = (0.d0, 0.d0)

    G_d_LL = (0.d0, 0.d0)
    G_m_LL = (0.d0, 0.d0)
    G_s_LL = (0.d0, 0.d0)
    G_t_LL = (0.d0, 0.d0)

!use so many for loops as one can keep the frequency fixed

    DO idx_q = 1, Nb
      DO l1 = 1, Nl
        DO kx1 = 1, Nx
          DO ky1 = 1, Ny
            DO w1 = 1, Nf
              DO l2 = 1, Nl
                DO kx2 = 1, Nx
                  DO ky2 = 1, Ny
                    DO w2 = 1, Nf

                      idxk1 = ((kx1 - 1)*Ny + (ky1 - 1))*Nf + w1
                      idxk2 = ((kx2 - 1)*Ny + (ky2 - 1))*Nf + w2
                      idxl1 = (l1 - 1)*Nf + w1
                      idxl2 = (l2 - 1)*Nf + w2

                      F_d_LL(idxl1, idxl2, idx_q) = F_d_LL(idxl1, idxl2, idx_q) + F_d(idxk1, idxk2, idx_q)* &
                                                    FF(kx1, ky1, l1)*FF(kx2, ky2, l2)
                      F_m_LL(idxl1, idxl2, idx_q) = F_m_LL(idxl1, idxl2, idx_q) + F_m(idxk1, idxk2, idx_q)* &
                                                    FF(kx1, ky1, l1)*FF(kx2, ky2, l2)

                      F_s_LL(idxl1, idxl2, idx_q) = F_s_LL(idxl1, idxl2, idx_q) + F_s(idxk1, idxk2, idx_q)* &
                                                    FF(kx1, ky1, l1)*FF(kx2, ky2, l2)
                      F_t_LL(idxl1, idxl2, idx_q) = F_t_LL(idxl1, idxl2, idx_q) + F_t(idxk1, idxk2, idx_q)* &
                                                    FF(kx1, ky1, l1)*FF(kx2, ky2, l2)

                      G_d_LL(idxl1, idxl2, idx_q) = G_d_LL(idxl1, idxl2, idx_q) + G_d(idxk1, idxk2, idx_q)* &
                                                    FF(kx1, ky1, l1)*FF(kx2, ky2, l2)
                      G_m_LL(idxl1, idxl2, idx_q) = G_m_LL(idxl1, idxl2, idx_q) + G_m(idxk1, idxk2, idx_q)* &
                                                    FF(kx1, ky1, l1)*FF(kx2, ky2, l2)

                      G_s_LL(idxl1, idxl2, idx_q) = G_s_LL(idxl1, idxl2, idx_q) + G_s(idxk1, idxk2, idx_q)* &
                                                    FF(kx1, ky1, l1)*FF(kx2, ky2, l2)
                      G_t_LL(idxl1, idxl2, idx_q) = G_t_LL(idxl1, idxl2, idx_q) + G_t(idxk1, idxk2, idx_q)* &
                                                    FF(kx1, ky1, l1)*FF(kx2, ky2, l2)

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  end subroutine trafo_2D

!-------------------------------------------------

!trasforms l-space vertices back to k-space
!this overwrites the original k-space vertices
!mainly for plotting results in k-space

  subroutine trafo_2D_inv

    integer :: kx1, ky1, kx2, ky2, l1, l2
    integer :: w1, w2
    integer :: idxk1, idxk2
    integer :: idxl1, idxl2
    integer :: idx_q

!initialize with zero again

    !k-space arrays needed for plotting
    IF (.NOT. allocated(F_m)) allocate (F_m(Nt, Nt, Nb))
    IF (.NOT. allocated(F_d)) allocate (F_d(Nt, Nt, Nb))
    IF (.NOT. allocated(F_s)) allocate (F_s(Nt, Nt, Nb))
    IF (.NOT. allocated(F_t)) allocate (F_t(Nt, Nt, Nb))
  
    IF (.NOT. allocated(G_m)) allocate (G_m(Nt, Nt, Nb))
    IF (.NOT. allocated(G_d)) allocate (G_d(Nt, Nt, Nb))
    IF (.NOT. allocated(G_s)) allocate (G_s(Nt, Nt, Nb))
    IF (.NOT. allocated(G_t)) allocate (G_t(Nt, Nt, Nb))


    F_d = (0.0d0, 0.0d0)
    F_m = (0.0d0, 0.0d0)
    F_s = (0.0d0, 0.0d0)
    F_t = (0.0d0, 0.0d0)

    G_d = (0.0d0, 0.0d0)
    G_m = (0.0d0, 0.0d0)
    G_s = (0.0d0, 0.0d0)
    G_t = (0.0d0, 0.0d0)

!use so many for loops as one can keep the frequency fixed
    DO idx_q = 1, Nb
      DO kx1 = 1, Nx
        DO ky1 = 1, Ny
          DO l1 = 1, Nl
            DO w1 = 1, Nf
              DO kx2 = 1, Nx
                DO ky2 = 1, Ny
                  DO l2 = 1, Nl
                    DO w2 = 1, Nf

                      idxk1 = ((kx1 - 1)*Ny + (ky1 - 1))*Nf + w1
                      idxk2 = ((kx2 - 1)*Ny + (ky2 - 1))*Nf + w2
                      idxl1 = (l1 - 1)*Nf + w1
                      idxl2 = (l2 - 1)*Nf + w2

                      F_d(idxk1, idxk2, idx_q) = F_d(idxk1, idxk2, idx_q) + F_d_LL(idxl1, idxl2, idx_q)* &
                                                 FF(kx1, ky1, l1)*FF(kx2, ky2, l2)
                      F_m(idxk1, idxk2, idx_q) = F_m(idxk1, idxk2, idx_q) + F_m_LL(idxl1, idxl2, idx_q)* &
                                                 FF(kx1, ky1, l1)*FF(kx2, ky2, l2)

                      F_s(idxk1, idxk2, idx_q) = F_s(idxk1, idxk2, idx_q) + F_s_LL(idxl1, idxl2, idx_q)* &
                                                 FF(kx1, ky1, l1)*FF(kx2, ky2, l2)
                      F_t(idxk1, idxk2, idx_q) = F_t(idxk1, idxk2, idx_q) + F_t_LL(idxl1, idxl2, idx_q)* &
                                                 FF(kx1, ky1, l1)*FF(kx2, ky2, l2)

                      G_d(idxk1, idxk2, idx_q) = G_d(idxk1, idxk2, idx_q) + G_d_LL(idxl1, idxl2, idx_q)* &
                                                 FF(kx1, ky1, l1)*FF(kx2, ky2, l2)
                      G_m(idxk1, idxk2, idx_q) = G_m(idxk1, idxk2, idx_q) + G_m_LL(idxl1, idxl2, idx_q)* &
                                                 FF(kx1, ky1, l1)*FF(kx2, ky2, l2)

                      G_s(idxk1, idxk2, idx_q) = G_s(idxk1, idxk2, idx_q) + G_s_LL(idxl1, idxl2, idx_q)* &
                                                 FF(kx1, ky1, l1)*FF(kx2, ky2, l2)
                      G_t(idxk1, idxk2, idx_q) = G_t(idxk1, idxk2, idx_q) + G_t_LL(idxl1, idxl2, idx_q)* &
                                                 FF(kx1, ky1, l1)*FF(kx2, ky2, l2)

                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  end subroutine trafo_2D_inv

!check orthogonality of basis functions
!
!  subroutine check_ortho
!
!    complex(dp), dimension(Nz, Nz) :: scal
!    integer :: i, j, kx, ky, nu1, nu2
!
!    character(len=30) :: label
!
!    scal = 0.0d0
!
!!the frequencies are just a meaningless nummy variable to be able to use the old plotfunction
!
!    DO nu1 = 1, Nf
!    DO nu2 = 1, Nf
!    DO i = 1, Nl
!      DO j = 1, Nl
!        !sum over 1. BZ
!        DO kx = 1, Nx
!          DO ky = 1, Ny
!
!            scal((i - 1)*Nf + nu1, (j - 1)*Nf + nu2) = &
!              scal((i - 1)*Nf + nu1, (j - 1)*Nf + nu2) + FF(kx, ky, i)*FF(kx, ky, j)
!
!          END DO
!        END DO
!      END DO
!    END DO
!    END DO
!    END DO
!
!    label = 'data/ortho.dat'
!    call output_L(scal, label)
!
!  end subroutine check_ortho

!--------------------------------------------------------------

  !symmetry operation in L-space - this is the inverse opeartion
  !for given symmetry and l-index <-> representation it returns the prefactor
  !for a vertex and the to be taken new l-index
  subroutine sym_op_L(is, l_in, l_sym, fac)

    !number of symmetry
    integer, intent(in)   :: is
    !original l-index
    integer, intent(in)   :: l_in
    !to be returned l-index
    integer, intent(out)  :: l_sym
    !symmetry factor
    real(dp), intent(out) :: fac

    !representation to which l_index correponds
    integer :: rep

    rep = reps(l_in)

    select case(rep)

      case(1:4)!1D representations

        l_sym = l_in
        fac = char_table(is, rep)

      case(5)
        if(l_in <= 5) then 
          select case(is)
          
            !l-index not changed
            case(1:4)
              l_sym = l_in
            !change l-index
            case(5:8)
              l_sym = l_in + 1

          end select
          
          !infact the same signs as in sym_op x-component (which is no coincidence)
          !6 and 8 are exchanged as this is the inverse operation
          select case(is)
            case (1) 
              fac = 1
            case (2) 
              fac = 1
            case (3) 
              fac = -1
            case (4) 
              fac = -1
            case (5) 
              fac = 1
            case (6) 
              fac = 1
            case (7) 
              fac = -1
            case (8) 
              fac = -1
          end select

        else

          select case(is)
            case (1) 
              fac = 1
              l_sym = l_in
            case (2) 
              fac = -1
              l_sym = l_in + 1
            case (3) 
              fac = 1
              l_sym = l_in + 1
            case (4) 
              fac = -1
              l_sym = l_in
            case (5) 
              fac = 1
              l_sym = l_in
            case (6) 
              fac = 1
              l_sym = l_in + 1
            case (7) 
              fac = -1
              l_sym = l_in
            case (8) 
              fac = -1
              l_sym = l_in + 1
          end select

        end if !l_in <= 5

      case(6)

        if(l_in .le. 5) then

          select case(is)
          
            !l-index not changed
            case(1:4)
              l_sym = l_in
            !change l-index
            case(5:8)
              l_sym = l_in - 1

          end select
    
          !infact the same signs as in sym_op y-component 
          !6 and 8 are exchanged since this is the inverse operaion
          select case(is)
            case (1) 
              fac = 1
            case (2) 
              fac = -1
            case (3) 
              fac = 1
            case (4) 
              fac = -1
            case (5) 
              fac = 1
            case (6) 
              fac = -1
            case (7) 
              fac = -1
            case (8) 
              fac = 1
          end select

        else

          select case(is)
            case (1) 
              fac = + 1
              l_sym = l_in
            case (2) 
              fac = - 1
              l_sym = l_in - 1
            case (3) 
              fac = + 1
              l_sym = l_in - 1
            case (4) 
              fac = - 1
              l_sym = l_in
            case (5) 
              fac = - 1
              l_sym = l_in
            case (6) 
              fac = - 1
              l_sym = l_in - 1
            case (7) 
              fac = + 1
              l_sym = l_in
            case (8) 
              fac = + 1
              l_sym = l_in - 1
          end select


        end if !l_in .le. 5

      end select 

  end subroutine sym_op_L

! ---------------------------------------------------------------

  !symmetry operations in the BZ
  !NOT ordered in the way platt orders them - see below
  subroutine symmetry_operation(i_sym, i_in, i_out)

    !
    !  Symmetry operation on one index map.
    !
    !  1  idetity
    !  2  (kx,ky) -> (kx,-ky)
    !  3  (kx,ky) -> (-kx,ky)
    !  4  (kx,ky) -> (-kx,-ky)
    !  5  (kx,ky) -> (ky,kx)
    !  6  (kx,ky) -> (ky,-kx) inv of 8
    !  7  (kx,ky) -> (-ky,-kx)
    !  8  (kx,ky) -> (-ky, kx) inv of 6

    integer, intent(in) :: i_sym
    type(indxmap), intent(in) :: i_in
    type(indxmap), intent(out) :: i_out

    integer :: ix, iy, wi, kx, ky

    ix = i_in%ix
    iy = i_in%iy
    wi = i_in%iw

    select case (i_sym)

    case (1)
      i_out = i_in

    case (2) ! (kx,ky) -> (kx,-ky)

      kx = ix
      ky = mod(-(iy - 1) + Ny, Ny) + 1

      i_out = indxmap(kx, ky, wi)

    case (3) ! (kx,ky) -> (-kx,ky)

      kx = mod(-(ix - 1) + Nx, Nx) + 1
      ky = iy

      i_out = indxmap(kx, ky, wi)

    case (4) ! (kx,ky) -> (-kx,-ky)

      ky = mod(-(iy - 1) + Ny, Ny) + 1
      kx = mod(-(ix - 1) + Nx, Nx) + 1

      i_out = indxmap(kx, ky, wi)

    case (5) ! (kx,ky) -> (ky,kx)

      kx = iy
      ky = ix

      i_out = indxmap(kx, ky, wi)

    case (6) !  (kx,ky) -> (ky,-kx)

      kx = iy
      ky = mod(-(ix - 1) + Nx, Nx) + 1

      i_out = indxmap(kx, ky, wi)

    case (7) ! (kx,ky) -> (-ky,-kx)

      kx = mod(-(iy - 1) + Ny, Ny) + 1
      ky = mod(-(ix - 1) + Nx, Nx) + 1

      i_out = indxmap(kx, ky, wi)

    case (8) ! (kx,ky) -> (-ky, kx)

      kx = mod(-(iy - 1) + Ny, Ny) + 1
      ky = ix

      i_out = indxmap(kx, ky, wi)

    end select

    !end if

  end subroutine symmetry_operation

!-----------------------------------------------------------

  subroutine symmetry_operation_inv(i_sym, i_in, j_in, i_out, j_out)

    !
    !  Inverse symmetry operation on two index maps.
    !
    !  1  idetity
    !  2  (kx,ky) -> (kx,-ky)
    !  3  (kx,ky) -> (-kx,ky)
    !  4  (kx,ky) -> (-kx,-ky)
    !  5  (kx,ky) -> (ky,kx)
    !  6  (kx,ky) -> (ky,-kx) inv of 8
    !  7  (kx,ky) -> (-ky,-kx)
    !  8  (kx,ky) -> (-ky, kx) inv of 6

    integer, intent(in) :: i_sym
    type(indxmap), intent(in) :: i_in, j_in
    type(indxmap), intent(out) :: i_out, j_out

    integer :: ix, iy, jx, jy, wi, wj, kx, ky

    ix = i_in%ix
    iy = i_in%iy
    wi = i_in%iw

    jx = j_in%ix
    jy = j_in%iy
    wj = j_in%iw

    select case (i_sym)

    case (1)
      i_out = i_in
      j_out = j_in

    case (2) ! (kx,ky) -> (kx,-ky)

      kx = ix
      ky = mod(-(iy - 1) + Ny, Ny) + 1
      !ky = -iy + Ny + 2
      !if (ky > Ny) ky = ky - Ny

      i_out = indxmap(kx, ky, wi)

      kx = jx
      ky = mod(-(jy - 1) + Ny, Ny) + 1
      !ky = -jy + Ny + 2
      !if (ky > Ny) ky = ky - Ny

      j_out = indxmap(kx, ky, wj)

    case (3) ! (kx,ky) -> (-kx,ky)

      kx = mod(-(ix - 1) + Nx, Nx) + 1
      !kx = -ix + Nx + 2
      !if (kx > Nx) kx = kx - Nx
      ky = iy

      i_out = indxmap(kx, ky, wi)

      kx = mod(-(jx - 1) + Nx, Nx) + 1
      !kx = -jx + Nx + 2
      !if (kx > Nx) kx = kx - Nx
      ky = jy

      j_out = indxmap(kx, ky, wj)

    case (4) ! (kx,ky) -> (-kx,-ky)

      kx = mod(-(ix - 1) + Nx, Nx) + 1
      !kx = -ix + Nx + 2
      !if (kx > Nx) kx = kx - Nx
      ky = mod(-(iy - 1) + Ny, Ny) + 1
      !ky = -iy + Ny + 2
      !if (ky > Ny) ky = ky - Ny

      i_out = indxmap(kx, ky, wi)

      kx = mod(-(jx - 1) + Nx, Nx) + 1
      !kx = -jx + Nx + 2
      !if (kx > Nx) kx = kx - Nx
      ky = mod(-(jy - 1) + Ny, Ny) + 1
      !ky = -jy + Ny + 2
      !if (ky > Ny) ky = ky - Ny

      j_out = indxmap(kx, ky, wj)

    case (5) ! (kx,ky) -> (ky,kx)

      kx = iy
      ky = ix

      i_out = indxmap(kx, ky, wi)

      kx = jy
      ky = jx

      j_out = indxmap(kx, ky, wj)

    case (6) ! inverse of 6 is 8 (kx,ky) -> (-ky, kx)

      kx = mod(-(iy - 1) + Ny, Ny) + 1
      !kx = -iy + Ny + 2
      !if (kx > Nx) kx = kx - Nx
      ky = ix

      i_out = indxmap(kx, ky, wi)

      kx = mod(-(jy - 1) + Ny, Ny) + 1
      !kx = -jy + Ny + 2
      !if (kx > Nx) kx = kx - Nx
      ky = jx

      j_out = indxmap(kx, ky, wj)

    case (7) ! (kx,ky) -> (-ky,-kx)

      kx = mod(-(iy - 1) + Ny, Ny) + 1
      !kx = -iy + Ny + 2
      !if (kx > Nx) kx = kx - Nx
      ky = mod(-(ix - 1) + Nx, Nx) + 1
      !ky = -ix + Nx + 2
      !if (ky > Ny) ky = ky - Ny

      i_out = indxmap(kx, ky, wi)

      kx = mod(-(jy - 1) + Ny, Ny) + 1
      !kx = -jy + Ny + 2
      !if (kx > Nx) kx = kx - Nx
      ky = mod(-(jx - 1) + Nx, Nx) + 1
      !ky = -jx + Nx + 2
      !if (ky > Ny) ky = ky - Ny

      j_out = indxmap(kx, ky, wj)

    case (8) ! inverse of 8 is 6: (kx,ky) -> (ky,-kx)

      kx = iy
      ky = mod(-(ix - 1) + Nx, Nx) + 1
      !ky = -ix + Nx + 2
      !if (ky > Ny) ky = ky - Ny

      i_out = indxmap(kx, ky, wi)

      kx = jy
      ky = mod(-(jx - 1) + Nx, Nx) + 1
      !ky = -jx + Nx + 2
      !if (ky > Ny) ky = ky - Ny

      j_out = indxmap(kx, ky, wj)

    end select

    !end if

  end subroutine symmetry_operation_inv

!---------------------------------------------

  subroutine list_symmetries(map, symm_list)

    !---------Purpose-------------
    !  This subroutine lists all operations that need to be performed to create
    !  star k for a given k. We need it to restore full BZ from IBZ
    !  The operations are numbered as follows (now square lattice only):
    !
    ! --- ONLY WORKS FOR EVEN LATTICES! ---
    !
    !  1  idetity
    !  2  (kx,ky) -> (kx,-ky)
    !  3  (kx,ky) -> (-kx,ky)
    !  4  (kx,ky) -> (-kx,-ky)
    !  5  (kx,ky) -> (ky,kx)
    !  6  (kx,ky) -> (ky,-kx)
    !  7  (kx,ky) -> (-ky,-kx)
    !  8  (kx,ky) -> (-ky, kx)

    type(indxmap), intent(in) :: map
    logical, intent(out) :: symm_list(Ns)

    integer ::  kx, ky

    kx = map%ix
    ky = map%iy

    symm_list = .False.

    if (((kx == 1) .and. (ky == 1)) .or. ((kx == Nx_IBZ) .and. (ky == Nx_IBZ))) then ! (0,0),(Pi,Pi)

      symm_list(1) = .True. ! identity

    else if ((kx == Nx_IBZ) .and. (ky == 1)) then !(Pi,0) (2x)

      symm_list(1) = .True. ! identity
      symm_list(5) = .True. ! 5  (kx,ky) -> (ky,kx)

    else if ((kx == Nx_IBZ) .and. ((ky > 1) .and. (ky < Nx_IBZ))) then !(Pi,j<Pi) (4x)

      symm_list(1) = .True. ! identity
      symm_list(2) = .True. !  2  (kx,ky) -> (kx,-ky)
      symm_list(5) = .True. ! 5  (kx,ky) -> (ky,kx)
      symm_list(8) = .True. !  8  (kx,ky) -> (-ky, kx)

    else if ((kx > 1) .and. (ky < Nx_IBZ)) then

      if (ky == 1) then !(0<i<Pi,0) (4x)

        symm_list(1) = .True. ! identity
        symm_list(3) = .True. !  3  (kx,ky) -> (-kx,ky)
        symm_list(5) = .True. !  5  (kx,ky) -> (ky,kx)
        symm_list(6) = .True. !  6  (kx,ky) -> (ky, -kx)

      else if (ky == kx) then ! (diagonal) (4x)

        symm_list(1) = .True. ! identity
        symm_list(2) = .True. !  2  (kx,ky) -> (kx,-ky)
        symm_list(3) = .True. !  3  (kx,ky) -> (-kx,ky)
        symm_list(4) = .True. !  4  (kx,ky) -> (-kx,-ky)

      else ! all other points (8x)

        symm_list = .True.

      end if
    end if

  end subroutine list_symmetries

!-------------------------------------------------------

  
  !assembel an array in whole BZ from one in IBZ
  !this is for only one index
  !num_freqs can be Nf or Nf/2 for fermionic or bosonic vertex
  subroutine symmetrize_array(num_freqs, arr_in, arr_out)

    !number of frequencies -> Nf for fermionic Nf/2 for bosonic
    integer, intent(in) :: num_freqs  
    !input - lives in IBZ - number of momenta gives as below
    complex(dp), dimension( (Nx_IBZ * (Nx_IBZ + 1))/2 * num_freqs), intent(in) :: arr_in
    !output lives in whole BZ
    complex(dp), dimension(Nx * Ny * num_freqs), intent(out) :: arr_out
 
    !momenta and frequencies in IBZ
    integer :: qx, qy, nu

    !maps to hold arguments and symmetrized arguments
    type(indxmap) :: map_q, map_q_sym
    !correspondung indices to dereference arrays
    integer :: idx_q, idx_q_sym  

    !loop variable for symmetries
    integer :: is
 
    !list for symmetrization of whole BZ
    logical, dimension(8) :: symm_list

    !array is assumed to hold no information at input
    arr_out = 0.0d0

      !loop over IBZ
      do qx = 1, Nx_IBZ
        do qy = 1,qx
          do nu = 1, num_freqs
           
            !create IBZ arguments
            map_q = indxmap(qx, qy, nu)
            idx_q = (((qx - 1) * qx)/2 + qy - 1) * num_freqs + nu
            
            call list_symmetries(map_q, symm_list)

            !loop over all symmetries
            do is = 1, Ns
              !to exclude double counting of high symmetry points 
              if (.not. symm_list(is)) cycle
                
              !generate symmetrized arguments
              call symmetry_operation(is, map_q, map_q_sym) 
              !calculates index in normal BZ
              idx_q_sym = ((map_q_sym%ix - 1) * Ny + map_q_sym%iy - 1) * num_freqs + map_q_sym%iw

              !calculate arr in whole IBZ
              arr_out(idx_q_sym) = arr_out(idx_q_sym) + arr_in(idx_q)
                
            end do !is -> symmetries

          end do !nu
        end do !qy
      end do !qx


  end subroutine symmetrize_array



  ! ---
  ! --- caculate green's functions ---
  pure function get_green(map_k) result(G)

    type(Indxmap), intent(in) :: map_k
    complex(dp) :: G

    complex(dp) :: nu

    !actual matsubara frequency
    nu = get_freq(map_k%iw)

    !proper inside/outside the box treatment
    if(map_k%iw < 1 .or. map_k%iw > Nf) then
 
      G = 1.0d0/(nu - Ek(map_k%ix, map_k%iy) + mu - Sigma_H(map_k%ix, map_k%iy))

    else

      G = Gkw( ( (map_k%ix - 1) * Ny + map_k%iy - 1 ) * Nf + map_k%iw ) 

    end if

   end function get_green

  ! ----
  !green_s functions with coarse graining
  pure function get_green_coarse(kx_g, ky_g, map_k) result(G)

    type(Indxmap), intent(in) :: map_k
    integer, intent(in) :: kx_g, ky_g
    complex(dp) :: G

    complex(dp) :: nu

    !actual matsubara frequency
    nu = get_freq(map_k%iw)

    !proper inside/outside the box treatment
    if(map_k%iw < 1 .or. map_k%iw > Nf) then
 
      G = 1.0d0/(nu - Ek_grain(kx_g, ky_g, map_k%ix, map_k%iy) + mu - &
          Sigma_H(map_k%ix, map_k%iy))

    else

      G = 1.0d0/(nu - Ek_grain(kx_g, ky_g, map_k%ix, map_k%iy) + mu - &
                 Sigma( ( (map_k%ix - 1) * Ny + map_k%iy - 1 ) * Nf + map_k%iw ) )

    end if

   end function get_green_coarse

  ! ------
  
  !function to get the actual matsubara frequency from an index
  pure function get_freq(nu) result(ret_nu)
    
    integer, intent(in) :: nu 
    complex(dp) :: ret_nu

    ret_nu = dcmplx(0.0d0, PI/beta * (2.0d0 * (nu - Nf/2 - 1) + 1.0d0))

  end function get_freq


! --- index operations ---

  !add a fermionic and a bosonic index map
  subroutine index_FaddB(map1, map2, map_out)

    type(Indxmap), intent(in)     :: map1, map2
    type(Indxmap), intent(out)    :: map_out

    ! ... local vars ...
    integer :: kx, ky, w

    kx = mod(map1%ix - 1 + map2%ix - 1, Nx) + 1

    ky = mod(map1%iy - 1 + map2%iy - 1, Ny) + 1

    w = map1%iw + map2%iw - 1

    map_out = indxmap(kx, ky, w)

  end subroutine index_FaddB

  !-------------------------------------------------------

  !negate a fermionic map
  subroutine index_minusF(map1, map_out)
    
    type(Indxmap), intent(in)     :: map1
    type(Indxmap), intent(out)    :: map_out

    integer :: kx, ky, w

    kx = mod( - (map1%ix - 1) + Nx, Nx) + 1

    ky = mod( - (map1%iy - 1) + Ny, Ny) + 1

    w = -map1%iw + Nf + 1

    map_out = indxmap(kx, ky, w)

  end subroutine index_minusF


end module parquet_util

