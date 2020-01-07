module parquet_plot

  use parquet_ini
  use parquet_formfactors 
  use parquet_util
  use hdf5_wrapper
  use hdf5

  implicit none

contains

!calculate susceptibility & symmetrize it 
subroutine symmetrize_sus(ite, converged)
    integer, intent(in) :: ite
    logical, intent(in) :: converged

  !loop index = 1, 2 for chosing s-wave or d-wave component
  integer :: form_factor

  if(id == master) then

    if(.not. allocated(sus_d)) allocate(sus_d(Nx * Ny * Nf/2, 4))
    if(.not. allocated(sus_m)) allocate(sus_m(Nx * Ny * Nf/2, 4))
    if(.not. allocated(sus_s)) allocate(sus_s(Nx * Ny * Nf/2, 4))
    if(.not. allocated(sus_t)) allocate(sus_t(Nx * Ny * Nf/2, 4))


    if(.not. allocated(chi0_ph_BZ)) allocate(chi0_ph_BZ(Nx * Ny * Nf/2))
    if(.not. allocated(chi0_pp_BZ)) allocate(chi0_pp_BZ(Nx * Ny * Nf/2))

    if(.not. allocated(chi0_ph_BZ_FF)) allocate(chi0_ph_BZ_FF(Nx * Ny * Nf/2, 4))
    if(.not. allocated(chi0_pp_BZ_FF)) allocate(chi0_pp_BZ_FF(Nx * Ny * Nf/2, 4))
  
    sus_d = 0.0d0
    sus_m = 0.0d0
    sus_s = 0.0d0
    sus_t = 0.0d0
  
    do form_factor = 1, 4
      !calculate sus_ph in whole BZ
      call symmetrize_array(Nf/2, sus_d_IBZ(:, form_factor), sus_d(:, form_factor))
      call symmetrize_array(Nf/2, sus_m_IBZ(:, form_factor), sus_m(:, form_factor))
      call symmetrize_array(Nf/2, sus_s_IBZ(:, form_factor), sus_s(:, form_factor))
      call symmetrize_array(Nf/2, sus_t_IBZ(:, form_factor), sus_t(:, form_factor))

      call symmetrize_array(Nf/2, chi0_ph_FF(:, form_factor), chi0_ph_BZ_FF(:, form_factor))
      call symmetrize_array(Nf/2, chi0_pp_FF(:, form_factor), chi0_pp_BZ_FF(:, form_factor))
    end do !form_factor
    
    !calculate chi0 in whole BZ
    call symmetrize_array(Nf/2, chi0_ph(:), chi0_ph_BZ(:))
    call symmetrize_array(Nf/2, chi0_pp(:), chi0_pp_BZ(:))

 
  end if !id == master

end subroutine symmetrize_sus

!-------------------------------------------

! stuff for hdf5 output

! ----------------

  !sets up the one main output file that contains all data and parameters
  !this function creates the file and saves the parameters from the calculation
  subroutine setup_hdf5()


    integer(hid_t) :: file_ident
    integer :: bool

    ! create file
    call hdf5_create_file(file_main)
    ! open file with file identifier ifile 
    call hdf5_open_file(file_main, file_ident)

    ! save general parameters
    call hdf5_write_data(file_ident, 'params/beta',         beta)    
    call hdf5_write_data(file_ident, 'params/xU',           xU)
    call hdf5_write_data(file_ident, 'params/Nx',           Nx)
    call hdf5_write_data(file_ident, 'params/Nf',           Nf)
    call hdf5_write_data(file_ident, 'params/Nl',           Nl)
    bool = merge(1, 0, DGA)
    call hdf5_write_data(file_ident, 'params/DGA',          bool)
    call hdf5_write_data(file_ident, 'params/Ngrain',       Ngrain)
    call hdf5_write_data(file_ident, 'params/f_range',      f_range)
    call hdf5_write_data(file_ident, 'params/nParticle',    nParticle)


    call hdf5_close_file(file_ident)


  end subroutine setup_hdf5


  !subroutine to write 'sigma-like' arrays into main file
  subroutine write_sigma(ite, converged)

    !current iteration
    integer, intent(in) :: ite
    logical, intent(in) :: converged
    character(len=20) :: ite_str, sig_str

    integer(hid_t) :: file_ident

    if(id == master) then

      ! open file with file identifier ifile 
      call hdf5_open_file(file_main, file_ident)

      !create one dataset with label specific for current iteration
      write (ite_str, '(I5.5)') ite
      sig_str = 'Sigma/sig_' // trim(ite_str)

      call hdf5_write_data(file_ident, sig_str, reshape(Sigma, (/ Nf, Nx * Ny /)))

      !also output splitup sigma in case of being at last iteration
      if(ite == ite_max .or. converged) then

        sig_str = 'Sigma/Final/sig'
        call hdf5_write_data(file_ident, sig_str, reshape(Sigma, (/ Nf, Nx * Ny /)))
 
        !sig_str = 'Sigma/Final/sig_U'
        !call hdf5_write_data(file_ident, sig_str, reshape(Sigma_U, (/ Nf, Nx * Ny /)))

        !sig_str = 'Sigma/Final/sigph'
        !call hdf5_write_data(file_ident, sig_str, reshape(Sigma_ph, (/ Nf, Nx * Ny /)))

        !sig_str = 'Sigma/Final/sigpb'
        !call hdf5_write_data(file_ident, sig_str, reshape(Sigma_pb, (/ Nf, Nx * Ny /)))

        !sig_str = 'Sigma/Final/sigpp'
        !call hdf5_write_data(file_ident, sig_str, reshape(Sigma_pp, (/ Nf, Nx * Ny /)))

      end if !ite == ite_max


      call hdf5_close_file(file_ident)

    end if! id == master

  end subroutine write_sigma


  !routine to output susceptibilities
  subroutine write_chis(ite, converged)

    !current iteration
    integer, intent(in) :: ite
    logical, intent(in) :: converged
    character(len=32) :: ite_str, chi_str

    integer(hid_t) :: file_ident

    if(id == master) then

      ! open file with file identifier ifile 
      call hdf5_open_file(file_main, file_ident)
  
      !create one dataset with label specific for current iteration
      write (ite_str, '(I5.5)') ite
  
      !Sus_d
      chi_str = 'Chi/chi_d_' // trim(ite_str)
      call hdf5_write_data(file_ident, chi_str, reshape(sus_d, (/ Nf/2, Nx * Ny, 4 /)))
  
      !Sus_m
      chi_str = 'Chi/chi_m_' // trim(ite_str)
      call hdf5_write_data(file_ident, chi_str, reshape(-sus_m, (/ Nf/2, Nx * Ny, 4 /)))
  
      !Sus_s
      chi_str = 'Chi/chi_s_' // trim(ite_str)
      call hdf5_write_data(file_ident, chi_str, reshape(sus_s, (/ Nf/2, Nx * Ny, 4 /)))

      !Sus_t
      chi_str = 'Chi/chi_t_' // trim(ite_str)
      call hdf5_write_data(file_ident, chi_str, reshape(sus_t, (/ Nf/2, Nx * Ny, 4 /)))
  
      !chi0
      chi_str = 'Chi/chi_0_' // trim(ite_str)
      call hdf5_write_data(file_ident, chi_str, reshape(-chi0_ph_BZ, (/ Nf/2, Nx * Ny /)))

      chi_str = 'Chi/chi_0_pp_' // trim(ite_str)
      call hdf5_write_data(file_ident, chi_str, reshape(-chi0_pp_BZ, (/ Nf/2, Nx * Ny /)))
 
      chi_str = 'Chi/chi_0_FF_' // trim(ite_str)
      call hdf5_write_data(file_ident, chi_str, reshape(-chi0_ph_BZ_FF, (/ Nf/2, Nx * Ny, 4 /)))

      chi_str = 'Chi/chi_0_pp_FF_' // trim(ite_str)
      call hdf5_write_data(file_ident, chi_str, reshape(-chi0_pp_BZ_FF, (/ Nf/2, Nx * Ny, 4 /)))
 
      if(ite == ite_max - 1 .or. converged) then
  
        !Sus_d
        chi_str = 'Chi/Final/chi_d'
        call hdf5_write_data(file_ident, chi_str, reshape(sus_d, (/ Nf/2, Nx * Ny, 4 /)))
    
        !Sus_m
        chi_str = 'Chi/Final/chi_m'
        call hdf5_write_data(file_ident, chi_str, reshape(-sus_m, (/ Nf/2, Nx * Ny, 4 /)))
    
        !Sus_s
        chi_str = 'Chi/Final/chi_s'
        call hdf5_write_data(file_ident, chi_str, reshape(sus_s, (/ Nf/2, Nx * Ny, 4 /)))

        !Sus_t
        chi_str = 'Chi/Final/chi_t'
        call hdf5_write_data(file_ident, chi_str, reshape(sus_t, (/ Nf/2, Nx * Ny, 4 /)))


        !chi0
        chi_str = 'Chi/Final/chi_0'
        call hdf5_write_data(file_ident, chi_str, reshape(-chi0_ph_BZ, (/ Nf/2, Nx * Ny /)))
 
        chi_str = 'Chi/Final/chi_0_pp'
        call hdf5_write_data(file_ident, chi_str, reshape(-chi0_pp_BZ, (/ Nf/2, Nx * Ny /)))
 
        chi_str = 'Chi/Final/chi_0_FF'
        call hdf5_write_data(file_ident, chi_str, reshape(-chi0_ph_BZ_FF, (/ Nf/2, Nx * Ny, 4 /)))

        chi_str = 'Chi/Final/chi_0_pp_FF'
        call hdf5_write_data(file_ident, chi_str, reshape(-chi0_pp_BZ_FF, (/ Nf/2, Nx * Ny, 4 /)))
 
      end if !ite == ite_max
  
      call hdf5_close_file(file_ident)

    end if! id == master

  end subroutine write_chis


  
  !routine to output full vertex
  subroutine write_F(ite)

    !current iteration
    integer, intent(in) :: ite
    character(len=16) :: F_str, id_str
    character(len = 120) :: file_F

    integer(hid_t) :: file_ident


    !create one dataset with label specific for current iteration
    write (id_str, '(I3.3)') id
    file_F = 'data/F_' // trim(kind_str) // '_' // trim(id_str) //'.hdf5'

    if(id == master) write(*, *) 'master writes F to ', file_F 
    if(id == ntasks - 1) write(*, *) 'last task writes F to ', file_F 

    !if(id == master) then

    ! create file
    call hdf5_create_file(file_F)
    call hdf5_open_file(file_F, file_ident)

    F_str = 'F/F_d'
    call hdf5_write_data(file_ident, F_str, reshape(F_d_LL, (/ Nf, Nl, Nf, Nl, Nb/)))

    F_str = 'F/F_m'
    call hdf5_write_data(file_ident, F_str, reshape(F_m_LL, (/ Nf, Nl, Nf, Nl, Nb/)))

    F_str = 'F/F_s'
    call hdf5_write_data(file_ident, F_str, reshape(F_s_LL, (/ Nf, Nl, Nf, Nl, Nb/)))

    F_str = 'F/F_t'
    call hdf5_write_data(file_ident, F_str, reshape(F_t_LL, (/ Nf, Nl, Nf, Nl, Nb/)))

    call hdf5_close_file(file_ident)

    !end if

  end subroutine write_F


  !routine to output full vertex
  subroutine write_G(ite)

    !current iteration
    integer, intent(in) :: ite
    character(len=16) :: G_str, id_str
    character(len = 120) :: file_G

    integer(hid_t) :: file_ident


    !create one dataset with label specific for current iteration
    write (id_str, '(I3.3)') id
    file_G = 'data/G_' // trim(kind_str) // '_' // trim(id_str) //'.hdf5'

    if(id == master) write(*, *) 'master writes G to ', file_G 
    if(id == ntasks - 1) write(*, *) 'last task writes G to ', file_G 

    !if(id == master) then

    ! create file
    call hdf5_create_file(file_G)
    call hdf5_open_file(file_G, file_ident)

    G_str = 'G/G_d'
    call hdf5_write_data(file_ident, G_str, reshape(G_d_LL, (/ Nf, Nl, Nf, Nl, Nb/)))

    G_str = 'G/G_m'
    call hdf5_write_data(file_ident, G_str, reshape(G_m_LL, (/ Nf, Nl, Nf, Nl, Nb/)))

    G_str = 'G/G_s'
    call hdf5_write_data(file_ident, G_str, reshape(G_s_LL, (/ Nf, Nl, Nf, Nl, Nb/)))

    G_str = 'G/G_t'
    call hdf5_write_data(file_ident, G_str, reshape(G_t_LL, (/ Nf, Nl, Nf, Nl, Nb/)))

    call hdf5_close_file(file_ident)

    !end if

  end subroutine write_G

end module parquet_plot
