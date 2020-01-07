module parquet_ini

  use mpi_mod

  implicit none


  !variables that specify the model

  !to be read in parameters
  real(dp) :: beta! = 5.d0
  real(dp) :: xU  ! = 4.d0
  integer :: Nx != 6! linear dimension
  integer :: Nf != 32! num of Matsubara frequencies
  integer :: neighbours != 1 !neighbours in BZ that can be specified
  integer :: ite_max !maximum iterations of parquet loop

  integer :: Ny !just set to Nx later
  real(dp) :: mu 
  real(dp) :: nParticle
  integer, parameter :: Nl_sus = 1 ! number of FF's for which susceptibitlity should be calculated
  integer :: Ngrain ! number of additonal k points in the coarse-grained k-sums
  integer, parameter :: Ns = 8 ! number of point group symmetries for a given lattice
  real(dp), parameter :: tprime = 0.0d0 !next nearest neighbour hopping
  integer, parameter :: Nf_0 = 65 !freqs for which bare chi should be plotted
  !boolean to decide wether to read to be taken formfactors from a file
  logical :: readinpick 

  !variables specific to the numerics

  real(dp), parameter :: f_sigma = real(1e-8, 8) !accuracy on sigma
  real(dp), parameter :: f_vert = real(1e-8, 8) !accuracy on sigma
  real(dp), parameter :: f_kernel = real(1e-8, 8) !accuracy on sigma
  integer :: f_range !specifies how many more frequencies one wants to sum in SDE
  real(dp):: f_damping !damping in BSE

  real(dp) :: relative_err = 1.0d0 !actual relative error on sigma
  logical :: Sigma_ini !boolean to specify whether to readin sigma
  logical :: use_sig_dmft !boolean to specify wether to use as local selfenergy DMFT selfenergy
  logical :: read_input !boolean to specify whether to readin vertices
  logical :: use_U2 !boolean to specify whether to use U2 correction
  logical :: save_F !boolean to specify whether to save F's
  logical :: oneshot !boolean to specify whether to calculate without selfenergy

  !variables derived from these

  integer :: Nx_IBZ 
  integer :: Nred
  integer :: NIBZ
  integer :: Nc
  integer :: Nl ! linear dimension formfactors
  integer :: Nt ! linear dimension of combined momentum and frequency
  integer :: Nz ! linear dimension of combined formfactor and frequency
  integer :: Nb! the number of bosonic variables on each node

  !things for DGA calculations
  logical :: DGA !determine whether to start from DMFT
  character(len=100) :: fname_L !name of input file for Lambda
  character(len=100) :: fname_L_ph, fname_L_pp

  ! -- for ladder calculations
  character(len=100) :: fname_gamloc_ph
  character(len=100) :: fname_gamloc_pp
  character(len=100) :: fname_sig_dmft
    

  !arrays

  real(dp), allocatable :: Ek(:, :) ! tight-binding dispersion
  real(dp), allocatable :: Ek_grain(:, :, :, :) ! tight-binding dispersion for a finer k grid

  complex(dp), allocatable :: Gkw(:) !Green's function in kw-space
  complex(dp), allocatable :: Grt(:, :, :) !Green's function in rt-space

  complex(dp), allocatable :: Sigma(:) ! Self-energy

  !split up contributions to selfenergy
  complex(dp), allocatable :: Sigma_U(:) ! Self-energy
  complex(dp), allocatable :: Sigma_ph(:) ! Self-energy
  complex(dp), allocatable :: Sigma_pb(:) ! Self-energy
  complex(dp), allocatable :: Sigma_pp(:) ! Self-energy

  complex(dp), allocatable :: Sigma_H(:, :) ! Hatree energy
  complex(dp), allocatable :: Sigma_DMFT(:) ! DMFT Self-energy

  complex(dp), allocatable :: SigmaOld(:) ! Self-energy from previous iteration

  ! complete vertex in each channel
  complex(dp), allocatable :: F_d(:, :, :) ! density channel
  complex(dp), allocatable :: F_m(:, :, :) ! magnetic channel
  complex(dp), allocatable :: F_s(:, :, :) ! singlet channel
  complex(dp), allocatable :: F_t(:, :, :) ! triplet channel

  ! irreducible vertex in each channel
  complex(dp), allocatable :: G_d(:, :, :) ! density channel
  complex(dp), allocatable :: G_m(:, :, :) ! magnetic channel
  complex(dp), allocatable :: G_s(:, :, :) ! singlet channel
  complex(dp), allocatable :: G_t(:, :, :) ! triplet channel

  ! -- ladder stuff

  ! local irreducible vertex in each channel
  complex(dp), allocatable :: G_d_loc(:, :, :) ! density channel
  complex(dp), allocatable :: G_m_loc(:, :, :) ! magnetic channel
  complex(dp), allocatable :: G_s_loc(:, :, :) ! singlet channel
  complex(dp), allocatable :: G_t_loc(:, :, :) ! triplet channel

  !local phis
  complex(dp), allocatable :: phi_m_loc(:, :, :)
  complex(dp), allocatable :: phi_d_loc(:, :, :)

  !phi_pp + Lambda for use in SDE
  complex(dp), allocatable :: Lphi(:, :, :)

  ! --

  !completey irreducible vertex
  complex(dp), allocatable :: L_d(:, :, :) ! density channel
  complex(dp), allocatable :: L_m(:, :, :) ! magnetic channel
  complex(dp), allocatable :: L_s(:, :, :) ! singlet channel
  complex(dp), allocatable :: L_t(:, :, :) ! triplet channel

  !array for storing temporary vertices in parquet
  complex(dp), allocatable :: mat(:, :, :)

  !L-space vertices
  ! complete vertex in each channel
  complex(dp), allocatable :: F_d_LL(:, :, :) ! density channel
  complex(dp), allocatable :: F_m_LL(:, :, :) ! magnetic channel
  complex(dp), allocatable :: F_s_LL(:, :, :) ! singlet channel
  complex(dp), allocatable :: F_t_LL(:, :, :) ! triplet channel

  ! irreducible vertex in each channel
  complex(dp), allocatable :: G_d_LL(:, :, :) ! density channel
  complex(dp), allocatable :: G_m_LL(:, :, :) ! magnetic channel
  complex(dp), allocatable :: G_s_LL(:, :, :) ! singlet channel
  complex(dp), allocatable :: G_t_LL(:, :, :) ! triplet channel

  !kernel functions for vertex asymptotics in l-space
  complex(dp), allocatable :: K2_d1(:, :, :)
  complex(dp), allocatable :: K2_d2(:, :, :)

  complex(dp), allocatable :: K2_m1(:, :, :)
  complex(dp), allocatable :: K2_m2(:, :, :)

  complex(dp), allocatable :: K2_s1(:, :, :)
  complex(dp), allocatable :: K2_s2(:, :, :)

  complex(dp), allocatable :: K2_t1(:, :, :)
  complex(dp), allocatable :: K2_t2(:, :, :)

  !save from previous iteration to check convergence
  complex(dp), allocatable :: K2_d1_prev(:, :, :)
  complex(dp), allocatable :: K2_m1_prev(:, :, :)
  complex(dp), allocatable :: K2_s1_prev(:, :, :)
  complex(dp), allocatable :: K2_t1_prev(:, :, :)

  real(dp), allocatable :: C_wave_x(:) ! work array for FFT in x
  real(dp), allocatable :: C_wave_y(:) ! work array for FFT in y

  !susceptibility for d and m channel - in whole BZ
  complex(dp), allocatable :: sus_d(:, :)
  complex(dp), allocatable :: sus_m(:, :)
  complex(dp), allocatable :: sus_s(:, :)
  complex(dp), allocatable :: sus_t(:, :)
  !susceptibility in d and m channel - only in IBZ
  complex(dp), allocatable :: sus_d_IBZ(:, :)
  complex(dp), allocatable :: sus_m_IBZ(:, :)
  complex(dp), allocatable :: sus_s_IBZ(:, :)
  complex(dp), allocatable :: sus_t_IBZ(:, :)

  ! same susceptibilities as above just including channel coupling
  complex(dp), allocatable :: sus_d_opt(:, :)
  complex(dp), allocatable :: sus_m_opt(:, :)
  complex(dp), allocatable :: sus_s_opt(:, :)
  complex(dp), allocatable :: sus_t_opt(:, :)
  complex(dp), allocatable :: sus_d_IBZ_opt(:, :)
  complex(dp), allocatable :: sus_m_IBZ_opt(:, :)
  complex(dp), allocatable :: sus_s_IBZ_opt(:, :)
  complex(dp), allocatable :: sus_t_IBZ_opt(:, :)

  complex(dp), allocatable :: sus_pp_L_IBZ(:, :)
  complex(dp), allocatable :: sus_pp_ph_IBZ(:, :)
  complex(dp), allocatable :: sus_pp_phb_IBZ(:, :)
  complex(dp), allocatable :: sus_pp_pp_IBZ(:, :)
  
  complex(dp), allocatable :: sus_pp_L(:, :)
  complex(dp), allocatable :: sus_pp_ph(:, :)
  complex(dp), allocatable :: sus_pp_phb(:, :)
  complex(dp), allocatable :: sus_pp_pp(:, :)

  !double occupancie
  complex(dp), allocatable :: double_occ_1P
  complex(dp), allocatable :: double_occ_2P

  !the same but with exact fourier transform
  complex(dp), allocatable :: chi0_ph(:) ! bubble diagram in p-h channel
  complex(dp), allocatable :: chi0_pp(:) ! bubble diagram in p-p channel
  complex(dp), allocatable :: chi0_ph_BZ(:) ! bubble in whole BZ
  complex(dp), allocatable :: chi0_pp_BZ(:) ! bubble in whole BZ

  !same as above just projected on different formfactors
  complex(dp), allocatable :: chi0_ph_FF(:, :)
  complex(dp), allocatable :: chi0_pp_FF(:, :)
  complex(dp), allocatable :: chi0_ph_BZ_FF(:, :)
  complex(dp), allocatable :: chi0_pp_BZ_FF(:, :)

  !local chi0 as needed for ladder DGA
  !saved as diagonal matrix in order to use GEMM 
  complex(dp), allocatable :: chi0_loc(:, :) 

!new data types - index-maps
!for k space arguments
  type :: indxmap
    integer :: ix
    integer :: iy
    integer :: iw
  end type indxmap


!for l-arguments
  type :: indxmap_L
    integer :: il
    integer :: iw
  end type indxmap_L

  ! --- variables for formfactors functionalities
  !array with the actual formfactors
  real(dp), dimension(:, :, :), allocatable :: FF
  !formfactors with inverse symmetry operation acted upon them
  real(dp), dimension(:, :, :, :), allocatable :: FF_inv
  !product of 4 formfactors needed in parquet
  real(dp), dimension(:, :, :, :, :, :, :), allocatable :: FFFF_safe
  !character table for the square lattice
  integer, dimension(:, :), allocatable :: char_table
  !this is the length of the character table
  integer, parameter :: length_char_table = 5
  !s and d-wave formfactors for susceptibility
  real(dp), dimension(:, :, :), allocatable :: FF_sdwave

  !list of formfactors to be read
  logical, dimension(:), allocatable :: pick
  
  !list of index maps to convert an integer to an index map
  type(indxmap), allocatable :: Index_fermionic(:), Index_bosonic(:), Index_bosonic_IBZ(:)
  !same in l-space
  type(indxmap_L), allocatable :: Index_L(:)

  !bond list in initialization of formfactors 
  integer, dimension(:, :), allocatable :: thelist
  !sign under trafo of vertices with respect to mirroring on origin
  integer, dimension(:), allocatable    :: signs
  !save which l-index corrsponds to which representation
  integer, dimension(:), allocatable    :: reps

  !to check convergence in vertices
  complex(dp) :: F_d_max
  complex(dp) :: F_m_max 
  complex(dp) :: F_s_max
  complex(dp) :: F_t_max

  !stuff for saving of Phi_R
  integer :: wperTask
  !number of lattice points that one saves Phis for
  integer :: NR
  !save Phi's as function of lattice argument
  complex(dp), allocatable :: PhiRd(:, :, :), PhiRm(:, :, :), PhiRs(:, :, :), PhiRt(:, :, :)
  !saved exponential function for FT 
  complex(dp), allocatable :: expSqR(:, :, :)
  complex(dp), allocatable :: expqSR(:, :, :)

  !global name to distiguish calculations
  character(len=100) :: kind_str

  !for plotting
  character(len=120) :: file_main

  !name of parameter file
  character(len=100) :: parameter_file
  !name of picklist file to pick specific formfactors
  character(len=100) :: pick_file

  ! --- decide on which version of SDE and PAE to take ---
  logical :: old
  logical :: SDE_old
  
contains

  subroutine set_parameters

    call get_command_argument(1, parameter_file)

    if(id == master) write(*, *) 'reading parameters from : ', TRIM(parameter_file)

    !parameter_file = 'prms_victory_Nl6'

    !now read in parameters
    open (unit=1, file=trim(parameter_file), status='old')

    read(1, *) beta
    read(1, *) xU
    read(1, *) nParticle
    mu = 0.0d0
    read(1, *) Nx
    Ny = Nx !only calculate square lattice atm
    read(1, *) Nf
    read(1, *) neighbours
    read(1, *) readinpick
    read(1, *) pick_file
    read(1, *) Ngrain
    read(1, *) f_range
    read(1, *) f_damping
    read(1, *) ite_max
    read(1, *) DGA
    read(1, *) fname_L
    read(1, *) use_sig_dmft
    read(1, *) fname_sig_dmft
    read(1, *) use_U2
    read(1, *) save_F
    read(1, *) oneshot
    read(1, *) old
    read(1, *) SDE_old

    close(1)

    Nx_IBZ = (Nx + 2)/2! linear dimension of the IBZ in x direction
    NIBZ = (Nx_IBZ * (Nx_IBZ + 1))/2 ! #momenta in IBZ
    Nred = (Nx_IBZ * (Nx_IBZ + 1))/2 * Nf/2! linear dimension of combined momentum and frequency
    Nc = Nx * Ny! total number of sites

    Sigma_ini = .false.
    read_input = .false.

    !number of bosonic arguments on each task
    !-> only determined at runtime due to need of ntasks
    Nb = Nred/ntasks  !<- for IBZ

    !number of used formfactors
    Nl = 0 !set later
    !number of formfactor and frequency arguments - fermionic
    Nz = 0 !set later
    !number of momentum and frequency arguments - fermionic
    Nt = Nx * Ny * Nf

    !set how many frequencies of Phi_R are saved on one Task
    wperTask = int((Nf - 1)/(2 * ntasks) + 1)
    !number of lattice points that one saves phis for
    NR = 1


    !fname_L = "data/inidata/Lambda_DMFT_b5_u4/Lambda_out"

    fname_gamloc_ph = "data/inidata/Gam_loc_B4_U4.init"
    fname_gamloc_pp = "data/inidata/Gam_loc_B4_U4.init"

    if(DGA) then

    if(id == master) write(*, *) "readind in : ", fname_L

    end if

    !set actual filenames from this
    fname_L_ph = TRIM(fname_L) // "_ph.init"
    fname_L_pp = TRIM(fname_L) // "_pp.init"

!=1 -> 1 FF
!=2 -> 5 FF
!=3 -> 9 FF
!=4 -> 13 FF       *
!=5 -> 21 FF     * *   <- Neighbours in IBZ
!=6 -> 25 FF   * * *


  end subroutine

  subroutine Memory_Release

    !  if (allocated(Delta)) Deallocate(Delta)

    ! clear the dispersion
    if (allocated(Ek)) Deallocate (Ek)
    !  if (allocated(Ek_grain)) Deallocate(Ek_grain)

    ! deallocate the Hatree energy
    !  if (allocated(Sigma_H)) Deallocate(Sigma_H)

    ! deallocate the index array
    if (allocated(Index_bosonic)) Deallocate (Index_bosonic)
    if (allocated(Index_bosonic_IBZ)) Deallocate (Index_bosonic_IBZ)
    if (allocated(Index_fermionic)) Deallocate (Index_fermionic)

    ! deallocate arries for FFT
    !  if (allocated(C_wave_x)) Deallocate(C_wave_x)
    !  if (allocated(C_wave_y)) Deallocate(C_wave_y)

    ! deallocate the local fully irreducible vertex
    if (allocated(L_d)) Deallocate (L_d)
    if (allocated(L_m)) Deallocate (L_m)
    if (allocated(L_s)) Deallocate (L_s)
    if (allocated(L_t)) Deallocate (L_t)

    ! deallocate the complete veretx
    if (allocated(F_d)) Deallocate (F_d)
    if (allocated(F_m)) Deallocate (F_m)
    if (allocated(F_s)) Deallocate (F_s)
    if (allocated(F_t)) Deallocate (F_t)

    ! deallocate the irreducible veretx in each channel
    if (allocated(G_d)) Deallocate (G_d)
    if (allocated(G_m)) Deallocate (G_m)
    if (allocated(G_s)) Deallocate (G_s)
    if (allocated(G_t)) Deallocate (G_t)

    ! deallocate the single-particle Green's function
    !  if (allocated(Gkw))     Deallocate(Gkw)
    if (allocated(Sigma)) Deallocate (Sigma)
    !  if (allocated(Sigma_DMFT))   Deallocate(Sigma_DMFT)

    !  if (allocated(Sigma_compare))   Deallocate(Sigma_compare)
    if (allocated(Chi0_ph)) Deallocate (Chi0_ph)
    if (allocated(Chi0_pp)) Deallocate (Chi0_pp)

!    if (allocated(K1_d))   Deallocate(K1_d)
!    if (allocated(K1_m))   Deallocate(K1_m)
!    if (allocated(K1_s))   Deallocate(K1_s)
!    if (allocated(K1_t))   Deallocate(K1_t)
!    if (allocated(K2_d1))  Deallocate(K2_d1)
!    if (allocated(K2_m1))  Deallocate(k2_m1)
!    if (allocated(K2_s1))  Deallocate(k2_s1)
!    if (allocated(K2_t1))  Deallocate(k2_t1)
!    if (allocated(K2_d2))  Deallocate(K2_d2)
!    if (allocated(K2_m2))  Deallocate(k2_m2)
!    if (allocated(K2_s2))  Deallocate(k2_s2)
!    if (allocated(K2_t2))  Deallocate(k2_t2)

    ! deallocate temparary array
    if (allocated(mat)) Deallocate (mat)

  end subroutine Memory_Release

end module parquet_ini
