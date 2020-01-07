!module to create arrays containing formfactors
module parquet_formfactors

  use parquet_ini

  implicit none

contains

!create the formfactor array
  subroutine FF_array

    type(indxmap), dimension(8, Nl) :: bonds
    type(indxmap) :: save_first
    integer :: l, i, kx, ky
!to first store the complex values
    complex(dp), dimension(Nx, Ny, Nl) :: FF_comp

    !norm of formfactors - needed in FFarray
    real(dp), dimension(:), allocatable :: norms

    real(dp) px, py

    IF (.NOT. allocated(FF)) allocate (FF(Nx, Ny, Nl))
    IF (.NOT. allocated(norms)) allocate (norms(Nl))

    FF = 0.0d0
    FF_comp = (0.0d0, 0.0d0)

!fill first element in bonds - on this symmetry operations will then act
    call initialize_bonds(Nl, bonds)

!fill all other elements including prefactor from character table
    DO l = 1, Nl

      save_first = bonds(1, l) !save since last argument is being overwritten

      DO i = 1, 8
        call my_symmetry_operation(i, save_first, bonds(i, l))
        !write character into iw argument
        bonds(i, l)%iw = char_table(i, save_first%iw)
      END DO

    END DO

    DO l = 1, Nl

      DO kx = 1, Nx
        DO ky = 1, Ny
          DO i = 1, 8
            px = 2.0d0/Nx*PI*(kx - 1)
            py = 2.0d0/Ny*PI*(ky - 1)
            FF_comp(kx, ky, l) = FF_comp(kx, ky, l) + &
                                 bonds(i, l)%iw * exp(dcmplx(0.0d0, bonds(i, l)%ix * px + bonds(i, l)%iy * py))

          END DO!i = symmetries
        END DO!kx
      END DO!ky

    END DO!l

!determine norm numerically - since this is the norm that is of interest
    norms = 0.0d0

!determinde norm
    DO l = 1, Nl
      DO kx = 1, Nx
        DO ky = 1, Ny
          norms(l) = norms(l) + abs(FF_comp(kx, ky, l))*abs(FF_comp(kx, ky, l))

        END DO
      END DO
    END DO

!contributions should be either purely real or purely imaginary
    DO l = 1, Nl
      DO kx = 1, Nx
        DO ky = 1, Ny

          FF(kx, ky, l) = 1.0d0/SQRT(norms(l))*real(FF_comp(kx, ky, l)) + &
                          1.0d0/SQRT(norms(l))*aimag(FF_comp(kx, ky, l))

        END DO!ky
      END DO!kx
    END DO!l

  !now also initialize FF-array with symmetries    
  call init_FF_inv

  end subroutine

!---------------------------------------------------------------

!initialize bonds by setting the first element
  subroutine initialize_bonds(Nl, bonds)
    integer, intent(in) :: Nl
    type(indxmap), dimension(8, Nl), intent(out) :: bonds

    !sym_pos - position in bonds
    !lat_pos - position in lattice
    !log_pos - position in pick
    integer :: b_pos, lat_pos, log_pos
    integer :: i

    !characters for (n,0) and edge bonds
    integer, dimension(4) :: chars_41
    !characters for (n,n) bonds
    integer, dimension(4) :: chars_42
    !characters for general bonds in interior
    integer, dimension(8) :: chars_8

    chars_41 = (/1, 4, 5, 5/)
    chars_42 = (/1, 3, 5, 5/)
    chars_8 = (/1, 2, 3, 4, 5, 5, 5, 5/)

    b_pos = 1
    lat_pos = 1
    log_pos = 1

    !allocate sign array and set to 1 at first
    IF (.NOT. allocated(signs)) allocate (signs(Nl))
    IF (.NOT. allocated(reps)) allocate (reps(Nl))
    signs = 1

    DO WHILE (.TRUE.)
      !origin
      IF (thelist(1, lat_pos) == 0 .AND. thelist(2, lat_pos) == 0) THEN

        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = thelist(1, lat_pos)
          bonds(1, b_pos)%iy = thelist(2, lat_pos)
          bonds(1, b_pos)%iw = 1

          reps(b_pos) = 1

          b_pos = b_pos + 1

          IF (b_pos > Nl) EXIT

        END IF

        log_pos = log_pos + 1
        lat_pos = lat_pos + 1

        CYCLE

        !corner - only different in even case
      ELSE IF (abs(thelist(1, lat_pos) - (Nx/2.0d0)) < 0.01 .AND. abs(thelist(2, lat_pos) - (Ny/2.0d0)) < 0.01) THEN

        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = thelist(1, lat_pos)
          bonds(1, b_pos)%iy = thelist(2, lat_pos)
          bonds(1, b_pos)%iw = 1

          reps(b_pos) = 1

          b_pos = b_pos + 1
          IF (b_pos > Nl) EXIT !should always be true in this case

        END IF
        log_pos = log_pos + 1
        lat_pos = lat_pos + 1

        CYCLE

        !edge with y = 0 - only different in even case
      ELSE IF (abs(thelist(1, lat_pos) - (Nx/2.0)) < 0.01 .AND. thelist(2, lat_pos) == 0) THEN

        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = thelist(1, lat_pos)
          bonds(1, b_pos)%iy = thelist(2, lat_pos)
          bonds(1, b_pos)%iw = 1 !only first ...

          reps(b_pos) = 1

          b_pos = b_pos + 1
          IF (b_pos > Nl) EXIT

        END IF

        log_pos = log_pos + 1

        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = thelist(1, lat_pos)
          bonds(1, b_pos)%iy = thelist(2, lat_pos)
          bonds(1, b_pos)%iw = 3 !... and third entry give non zero result - in even case

          reps(b_pos) = 3

          b_pos = b_pos + 1
          IF (b_pos > Nl) EXIT

        END IF

        log_pos = log_pos + 1
        lat_pos = lat_pos + 1

        CYCLE

        !diagonal
      ELSE IF (thelist(1, lat_pos) == thelist(2, lat_pos)) THEN

        DO i = 1, 3
          IF (pick(log_pos)) THEN
            bonds(1, b_pos)%ix = thelist(1, lat_pos)
            bonds(1, b_pos)%iy = thelist(2, lat_pos)
            bonds(1, b_pos)%iw = chars_41(i)

            reps(b_pos) = chars_41(i)

            IF (i == 3) signs(b_pos) = -1

            b_pos = b_pos + 1
            IF (b_pos > Nl) EXIT !only exits inner loop

          END IF

          log_pos = log_pos + 1

        END DO
        IF (b_pos > Nl) EXIT !one also has to exit outer while loop in this case

        !start with different bond for 2D representation
        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = -thelist(2, lat_pos)
          bonds(1, b_pos)%iy = thelist(1, lat_pos)
          bonds(1, b_pos)%iw = chars_41(4)

          reps(b_pos) = 6!the second bond for 2D rep

          signs(b_pos) = -1

          b_pos = b_pos + 1
          IF (b_pos > Nl) EXIT

        END IF
        log_pos = log_pos + 1
        lat_pos = lat_pos + 1

        CYCLE

        !vertical bonds (n,0) or rest of edge in even case
      ELSE IF ((thelist(2, lat_pos) == 0) .OR. abs(thelist(1, lat_pos) - (Nx/2.0)) < 0.01) THEN

        DO i = 1, 3
          IF (pick(log_pos)) THEN
            bonds(1, b_pos)%ix = thelist(1, lat_pos)
            bonds(1, b_pos)%iy = thelist(2, lat_pos)
            bonds(1, b_pos)%iw = chars_42(i)

            reps(b_pos) = chars_42(i)

            IF (i == 3) signs(b_pos) = -1

            b_pos = b_pos + 1
            IF (b_pos > Nl) EXIT

          END IF

          log_pos = log_pos + 1

        END DO
        IF (b_pos > Nl) EXIT !one also has to exit outer while loop in this case

        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = -thelist(2, lat_pos)
          bonds(1, b_pos)%iy = thelist(1, lat_pos)
          bonds(1, b_pos)%iw = chars_42(i)

          reps(b_pos) = 6

          signs(b_pos) = -1

          b_pos = b_pos + 1
          IF (b_pos > Nl) EXIT

        END IF
        log_pos = log_pos + 1
        lat_pos = lat_pos + 1

        CYCLE

        !8 equivalent neighbours
      ELSE

        DO i = 1, 5
          IF (pick(log_pos)) THEN
            bonds(1, b_pos)%ix = thelist(1, lat_pos)
            bonds(1, b_pos)%iy = thelist(2, lat_pos)
            bonds(1, b_pos)%iw = chars_8(i)

            reps(b_pos) = chars_8(i)

            IF (i == 5) signs(b_pos) = -1

            b_pos = b_pos + 1
            IF (b_pos > Nl) EXIT

          END IF

          log_pos = log_pos + 1

        END DO
        IF (b_pos > Nl) EXIT !one also has to exit outer while loop in this case

        !now choose different starting bonds for 2D representation

        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = thelist(1, lat_pos)
          bonds(1, b_pos)%iy = -thelist(2, lat_pos)
          bonds(1, b_pos)%iw = chars_8(6)

          reps(b_pos) = 6

          signs(b_pos) = -1

          b_pos = b_pos + 1
          IF (b_pos > Nl) EXIT

        END IF
        log_pos = log_pos + 1

        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = thelist(2, lat_pos)
          bonds(1, b_pos)%iy = thelist(1, lat_pos)
          bonds(1, b_pos)%iw = chars_8(7)

          reps(b_pos) = 5

          signs(b_pos) = -1

          b_pos = b_pos + 1
          IF (b_pos > Nl) EXIT

        END IF
        log_pos = log_pos + 1

        IF (pick(log_pos)) THEN
          bonds(1, b_pos)%ix = -thelist(2, lat_pos)
          bonds(1, b_pos)%iy = thelist(1, lat_pos)
          bonds(1, b_pos)%iw = chars_8(8)

          reps(b_pos) = 6

          signs(b_pos) = -1

          b_pos = b_pos + 1
          IF (b_pos > Nl) EXIT

        END IF
        log_pos = log_pos + 1
        lat_pos = lat_pos + 1

        CYCLE

      END IF

    END DO

  end subroutine

!-------------------------------------------------

!not from victory code
!since checking modulu is cumbersome
!but the one from victory should work when subtracting 1's in the correct places
  subroutine my_symmetry_operation(i_sym, i_in, i_out)

    !  Symmetry operation on one index map.
    !
    !  1  idetity
    !  2  (kx,ky) -> (kx,-ky) - sig_v
    !  3  (kx,ky) -> (-kx,ky) - sig_v
    !  4  (kx,ky) -> (-kx,-ky) - C2
    !  5  (kx,ky) -> (ky,kx) - sig_d
    !  6  (kx,ky) -> (ky,-kx) inv of 8 - C4
    !  7  (kx,ky) -> (-ky,-kx) - sig_d
    !  8  (kx,ky) -> (-ky, kx) inv of 6 - C4

    integer, intent(in) :: i_sym
    type(indxmap), intent(in) :: i_in
    type(indxmap), intent(out) :: i_out

    integer :: ix, iy, info, kx, ky

    ix = i_in%ix
    iy = i_in%iy
    info = i_in%iw

    select case (i_sym)

    case (1)
      i_out = i_in

    case (2) ! (kx,ky) -> (kx,-ky)

      kx = ix
      ky = -iy

      i_out = indxmap(kx, ky, info)

    case (3) ! (kx,ky) -> (-kx,ky)

      kx = -ix
      ky = iy

      i_out = indxmap(kx, ky, info)

    case (4) ! (kx,ky) -> (-kx,-ky)

      ky = -iy
      kx = -ix

      i_out = indxmap(kx, ky, info)

    case (5) ! (kx,ky) -> (ky,kx)

      kx = iy
      ky = ix

      i_out = indxmap(kx, ky, info)

    case (6) !  (kx,ky) -> (ky,-kx)

      kx = iy
      ky = -ix

      i_out = indxmap(kx, ky, info)

    case (7) ! (kx,ky) -> (-ky,-kx)

      kx = -iy
      ky = -ix

      i_out = indxmap(kx, ky, info)

    case (8) ! (kx,ky) -> (-ky, kx)

      kx = -iy
      ky = ix

      i_out = indxmap(kx, ky, info)

    end select

  end subroutine my_symmetry_operation

!--------------------------------------------

!here variables important for the initialization of the Formfactor array are set
  recursive subroutine Init_Formfactor_variables

    integer :: i, stat
    logical :: bool

!determine length of pick

    IF (.NOT. allocated(pick)) allocate (pick(Nx*Ny))

    Nl = 0

!first initialize list of matrix elements
    call create_list

    IF (readinpick) THEN

      if(id == master) write(*, *) "reading picklist from: ", pick_file

      open (unit=2, file=pick_file, IOSTAT=stat, status='old')
      !if file doesnt exist for instance one just takes the specified number of neighbours
      IF (stat > 0) THEN
        write (*, *) 'Cannot open picklist file for reading - continuing with specified number of neighbours'
        readinpick = .false.
        close (2)

        call Init_Formfactor_variables

        !rest of function not needed in this case
        return

      ELSE

        DO i = 1, Nx*Ny
          read (2, *) bool
          pick(i) = bool
          !determine Nl on the fly
          IF (bool .EQV. .TRUE.) Nl = Nl + 1
        END DO

      END IF

      close (2)

    ELSE
!create logical list for the specified number of neighbours
      DO i = 1, neighbours
        !origin
        IF (thelist(1, i) == 0 .AND. thelist(2, i) == 0) THEN
          Nl = Nl + 1
          !corner - only different in even case
        ELSE IF (abs(thelist(1, i) - (Nx/2.0)) < 0.01 .AND. abs(thelist(2, i) - (Ny/2.0)) < 0.01) THEN
          Nl = Nl + 1
          !right with y = 0 - only different in even case
        ELSE IF (abs(thelist(1, i) - (Nx/2.0)) < 0.01 .AND. thelist(2, i) == 0) THEN
          Nl = Nl + 2
          !edge - only different in even case
        ELSE IF (abs(thelist(1, i) - (Nx/2.0)) < 0.01) THEN
          Nl = Nl + 4
          !four equivalent neighbours
        ELSE IF (thelist(1, i) == thelist(2, i) .OR. thelist(2, i) == 0) THEN
          Nl = Nl + 4
          !8 equivalent neighbours
        ELSE
          NL = NL + 8
        END IF
      END DO

!total number of formfactor functions
      IF (.NOT. allocated(pick)) allocate (pick(Nx*Ny))

      !now fill in first Nl elements of
      DO i = 1, Nx*Ny
        IF (i <= Nl) THEN
          pick(i) = .true.
        ELSE
          pick(i) = .false.
        END IF
      END DO

      !if of init_picklist
    END IF

!set Nz after Nl has been determined

    Nz = Nl*Nf

  end subroutine

!------------------------------------------


!this subroutine initializes the array FF_inv
!it takes the old FF-array and fills in the functions with all
!INVERSE symmetry opeartions acted upon
  subroutine init_FF_inv

    integer :: idx_l, idx_kx, idx_ky, s
    type(indxmap) :: map_k, map_k_sym
   
    IF(.NOT. allocated(FF_inv)) allocate(FF_inv(Ns, Nx, Ny, Nl)) 

    DO idx_l = 1, Nl
      DO idx_ky  = 1, Ny
        DO idx_kx  = 1, Nx
          DO s = 1, Ns

            !create an indxmap just with frequency zero
            map_k = Index_Fermionic(((idx_kx-1) * Ny + (idx_ky-1)) * Nf + 1)

            call my_symmetry_operation_inv(s, map_k, map_k, map_k_sym, map_k_sym)
            !take accoding value from existing array
            FF_inv(s, idx_kx, idx_ky, idx_l) = FF(map_k_sym%ix, map_k_sym%iy, idx_l) 

          END DO
        END DO
      END DO
    END DO

  end subroutine init_FF_inv

!create character table - for square lattice
  subroutine init_char_table(length_char_table, Ns)
    integer, intent(in) :: length_char_table
    integer, intent(in) :: Ns

    IF (.NOT. allocated(char_table)) allocate (char_table(Ns, length_char_table))

    !order of symmetry operations as in victory code
    !NOT as in platt review
    char_table = &
      reshape((/1, 1, 1, 1, 1, 1, 1, 1, &
                1, -1, -1, 1, -1, 1, -1, 1, &
                1, 1, 1, 1, -1, -1, -1, -1, &
                1, -1, -1, 1, 1, -1, 1, -1, &
                2, 0, 0, -2, 0, 0, 0, 0/), &
              (/Ns, length_char_table/))

  end subroutine

!------------------------------------------

!creates a list of actual lattice positions in IBZ
!this is important for the formfactor initialization
  subroutine create_list

    integer i, j
    integer :: pos, many

    many = ((Nx/2 + 1)*((Nx/2 + 1) + 1))/2 ! watch paranthesis to avoid integer rounding

    IF (.NOT. allocated(thelist)) allocate (thelist(3, many))

    pos = 1

    DO i = 1, Nx/2 + 1
      DO j = 1, i
        thelist(1, pos) = i - 1
        thelist(2, pos) = j - 1
        thelist(3, pos) = (i - 1)*(i - 1) + (j - 1)*(j - 1)
        pos = pos + 1
      END DO
    END DO

    call sort(many, thelist)

    !if (id == master) then

    !  DO i = 1, many
    !    write (*, *) thelist(1, i), " ", thelist(2, i), " ", thelist(3, i)
    !  END DO

    !end if

  end subroutine

!-------------------------------------------------------------

!create a sorted lattice cite list
  subroutine sort(many, thelist)
    integer, intent(in) :: many
    integer, dimension(3, many), intent(inout) :: thelist

    integer :: i, j
    integer, dimension(3) :: temp

    DO i = 1, many
      DO j = i, 2, -1
        IF (thelist(3, j) < thelist(3, j - 1)) THEN
          temp = thelist(:, j - 1)
          thelist(:, j - 1) = thelist(:, j)
          thelist(:, j) = temp(:)
        END IF
      END DO
    END DO

  end subroutine


  !this is absolutely the same is in parquet util
  !I need to fix the circular dependency at some point ... dont know how atm
  subroutine my_symmetry_operation_inv(i_sym, i_in, j_in, i_out, j_out)

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
      ky = -iy + Ny + 2
      if (ky > Ny) ky = ky - Ny

      i_out = indxmap(kx, ky, wi)

      kx = jx
      ky = -jy + Ny + 2
      if (ky > Ny) ky = ky - Ny

      j_out = indxmap(kx, ky, wj)

    case (3) ! (kx,ky) -> (-kx,ky)

      kx = -ix + Nx + 2
      if (kx > Nx) kx = kx - Nx
      ky = iy

      i_out = indxmap(kx, ky, wi)

      kx = -jx + Nx + 2
      if (kx > Nx) kx = kx - Nx
      ky = jy

      j_out = indxmap(kx, ky, wj)

    case (4) ! (kx,ky) -> (-kx,-ky)

      kx = -ix + Nx + 2
      if (kx > Nx) kx = kx - Nx
      ky = -iy + Ny + 2
      if (ky > Ny) ky = ky - Ny

      i_out = indxmap(kx, ky, wi)

      kx = -jx + Nx + 2
      if (kx > Nx) kx = kx - Nx
      ky = -jy + Ny + 2
      if (ky > Ny) ky = ky - Ny

      j_out = indxmap(kx, ky, wj)

    case (5) ! (kx,ky) -> (ky,kx)

      kx = iy
      ky = ix

      i_out = indxmap(kx, ky, wi)

      kx = jy
      ky = jx

      j_out = indxmap(kx, ky, wj)

    case (6) ! inverse of 6 is 8 (kx,ky) -> (-ky, kx)

      kx = -iy + Ny + 2
      if (kx > Nx) kx = kx - Nx
      ky = ix

      i_out = indxmap(kx, ky, wi)

      kx = -jy + Ny + 2
      if (kx > Nx) kx = kx - Nx
      ky = jx

      j_out = indxmap(kx, ky, wj)

    case (7) ! (kx,ky) -> (-ky,-kx)

      kx = -iy + Ny + 2
      if (kx > Nx) kx = kx - Nx
      ky = -ix + Nx + 2
      if (ky > Ny) ky = ky - Ny

      i_out = indxmap(kx, ky, wi)

      kx = -jy + Ny + 2
      if (kx > Nx) kx = kx - Nx
      ky = -jx + Nx + 2
      if (ky > Ny) ky = ky - Ny

      j_out = indxmap(kx, ky, wj)

    case (8) ! inverse of 8 is 6: (kx,ky) -> (ky,-kx)

      kx = iy
      ky = -ix + Nx + 2
      if (ky > Ny) ky = ky - Ny

      i_out = indxmap(kx, ky, wi)

      kx = jy
      ky = -jx + Nx + 2
      if (ky > Ny) ky = ky - Ny

      j_out = indxmap(kx, ky, wj)

    end select

    !end if

  end subroutine my_symmetry_operation_inv

!fills seperate array with s and d-wave formfactors
subroutine fill_FF_sdwave()

  !loop variables
  integer kx, ky, l
  real(dp) px, py

  if(.not. allocated(FF_sdwave)) allocate(FF_sdwave(Nx, Ny, 4))
  
  FF_sdwave(:, :, 1) = 1.0d0

  do ky = 1, Nx
    do kx = 1, Ny

      px = 2.0d0/Nx*PI*(kx - 1)
      py = 2.0d0/Ny*PI*(ky - 1)

      FF_sdwave(kx, ky, 2) = cos(px) - cos(py)
      FF_sdwave(kx, ky, 3) = sin(px)
      FF_sdwave(kx, ky, 4) = sin(px) * sin(py)

    end do !kx
  end do !ky


end subroutine fill_FF_sdwave

end module parquet_formfactors


