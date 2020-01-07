module parquet_equation

  use mpi_mod
  use parquet_ini
  use parquet_util
  use parquet_kernel


  implicit none

contains

!solver for parquet equations in l-space
  subroutine solve_parquet_equation(ite)
    integer, intent(in) :: ite

    integer :: idx_l1, idx_l2, idx_l3, idx_l4, idx_q1, idx_q2, idx_q2_mom
    type(indxmap_L) :: map_l3, map_l4
    !for initialization of F with Lambda
    integer :: idx_q, q_0
    !bosonic index map ... _sym for map in whole BZ
    type(indxmap)   :: map_q1, map_q2
    !for FFFF_safe call
    integer :: q_ind, q1x_prev, q1y_prev

    integer :: w1, w2, w3, w4, l1, l2, l3, l4

    integer :: inode, ichannel
    !integer :: id, ntasks

    real(dp) :: f_temp, aux

    !for timing initialization of FFFF_safe
    integer(4) :: count_start, count_end
    integer(4) :: count_rate, count_max

    !write(*, *) "starting parquet with task : ", id 

    !initialize FFFF_safe
    call system_clock(count_start, count_rate, count_max)

    if(ite == 0) call FFFF_mem

    call system_clock(count_end, count_rate, count_max)
    if (id == master .and. ite == 0) write (*, "(a, f12.6)") ' FFFF init took: ', &
                                              dble(count_end - count_start)/dble(count_rate)


    IF (.NOT. allocated(mat)) allocate (mat(Nz, Nz, Nb))

    ! --- reducible vertex + fully irreducible vertex
    F_d_LL = G_d_LL
    F_m_LL = G_m_LL
    F_s_LL = G_s_LL
    F_t_LL = G_t_LL

    !only local arguments have a contribution from local vertex
    !but one has to inlcude an appropriate normalization
    aux = abs( 1.0d0/FF(1, 1, 1) * 1.0d0/FF(1, 1, 1) )

    do idx_q = 1, Nb

      q_0 = index_bosonic_IBZ(id * Nb + idx_q)%iw

      F_d_LL(1:Nf, 1:Nf, idx_q) = F_d_LL(1:Nf, 1:Nf, idx_q) + L_d(1:Nf, 1:Nf, q_0) * aux 
      F_m_LL(1:Nf, 1:Nf, idx_q) = F_m_LL(1:Nf, 1:Nf, idx_q) + L_m(1:Nf, 1:Nf, q_0) * aux
      F_s_LL(1:Nf, 1:Nf, idx_q) = F_s_LL(1:Nf, 1:Nf, idx_q) + L_s(1:Nf, 1:Nf, q_0) * aux
      F_t_LL(1:Nf, 1:Nf, idx_q) = F_t_LL(1:Nf, 1:Nf, idx_q) + L_t(1:Nf, 1:Nf, q_0) * aux

    end do

    ! Important note:
    ! ===============
    ! The communication between nodes needed here is provided by BCAST
    ! Depending on the MPI version, the message size in BCAST is limited to ~2GB (2^31 bytes);
    ! In case you need to pass larger messages, it is necessary to split them in smaller chunks.
    ! Alternatively, if your Nb >1 you may use more cores for paralelization.
    ! Newest versions of MPI may support larger messages;
    ! in this case, please remove the following line of code.

    if (Nz*Nz*Nb > 134217728) write (*, *) 'The MPI message size exceeds the 2^31 byte limit.'

    DO ichannel = 1, 4

      DO inode = 0, ntasks - 1

        SELECT CASE (ichannel)
        CASE (1)
          mat = G_d_LL ! (1) density channel
        CASE (2)
          mat = G_m_LL ! (2) magnetic channel
        CASE (3)
          mat = G_s_LL ! (3) singlet channel
        CASE (4)
          mat = G_t_LL ! (4) triplet channel
        END SELECT

        !write(*, *) "before cast"

        call MPI_BCAST(mat, Nz * Nz * Nb, MPI_DOUBLE_COMPLEX, inode, MPI_COMM_WORLD, rc)

        !write(*, *) "after cast"

        SELECT CASE (ichannel)
        CASE (1, 2) ! density and magnetic channel

          !  rotation 1: Phi(k, k+q; k'-k) -> Phi(k1, k1'; q1)

          DO idx_q1 = 1, Nb

            !map from 'node-independent' index
            map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)

            ! --- some stuff to loop momenta independent of frequencies

            if(idx_q1 == 1) then

              !variables to loop momenta independent of frequencies
              q1x_prev = 0
              q1y_prev = 0
              q_ind = 0

            end if

            !only advance in calculation if momentum changes
            if ( (map_q1%ix .NE. q1x_prev) .OR. (map_q1%iy .NE. q1y_prev) ) then
      
              q1x_prev = map_q1%ix
              q1y_prev = map_q1%iy
      
              q_ind = q_ind + 1
            
            end if

            ! --- end of momentum loop stuff

            DO idx_q2 = 1, Nb

              map_q2 = Index_Bosonic_IBZ(inode * Nb + idx_q2)

              !get momentum part only of index
              idx_q2_mom = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1


                do l4 = 1, Nl
                  do l2 = 1, Nl
                    do l1 = 1, Nl
                      do l3 = 1, Nl
                        
                        !calculate product of four formfactors
                        !f_temp = FFFF(map_q1%ix, map_q1%iy, map_q2%ix, map_q2%iy, &
                        !              l1, l2, l3, l4, 1)

                        !load product of four formactors from saved value
                        f_temp = FFFF_safe(l3, l1, l2, l4, idx_q2_mom, q_ind, 1)

                        !fermionic frequency 
                        do w1 = 1, Nf

                          idx_l1 = (l1 - 1) * Nf + w1
                          w2 = w1 + map_q2%iw - 1 !nu_2 = nu_1 + w2
  
                          IF (w2 >= 1 .AND. w2 <= Nf) THEN !check if inside frequency box
                            idx_l2 = (l2 - 1) * Nf + w2
                            w4 = w1 + map_q1%iw - 1 !nu_4 = nu_1 + w1
  
                            IF (w4 >= 1 .AND. w4 <= Nf) THEN !check if inside frequency box
  
                              idx_l4 = (l4 - 1) * Nf + w4
                              idx_l3 = (l3 - 1) * Nf + w1 !List_Index_L(map_l3)
  
                              SELECT CASE (ichannel)
  
                              CASE (1) !density
                                F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-0.5d0)*mat(idx_l3, idx_l4, idx_q2) * &
                                                                 f_temp
                                F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                                 f_temp
  
                              CASE (2) !magnetic
                                F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-1.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                                 f_temp
                                F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                                 f_temp
  
                              END SELECT
  
                            !kernels
                            ELSE
 
                              map_l3 = Indxmap_L(l3, w1)
                              map_l4 = Indxmap_L(l4, w4)
  
                              SELECT CASE (ichannel)
  
                              CASE (1) !density
                                F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-0.5d0)*kernel('d', map_l3, map_l4, map_q2)* &
                                                                 f_temp
                                F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-0.5d0)*kernel('d', map_l3, map_l4, map_q2)* &
                                                                 f_temp
  
                              CASE (2) !magnetic
                                F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-1.5d0)*kernel('m', map_l3, map_l4, map_q2)* &
                                                                 f_temp
                                F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (0.5d0)*kernel('m', map_l3, map_l4, map_q2)* &
                                                                 f_temp
 
                              END SELECT
  
                            END IF !w4 inside box
                          END IF !w2 inside box
  
                        end do !w1

                      END DO !l4
                    end do !l3
                  end do !l2
                end do !l1
                !!!$omp end parallel do
            END DO !idx_q2
          END DO !idx_q1
        ! END DO !idx_l1


          !Now time - reversal part

          DO idx_q1 = 1, Nb

            map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)    !make sure this is correct on any node

            ! --- some stuff to loop momenta independent of frequencies

            if(idx_q1 == 1) then

              !variables to loop momenta independent of frequencies
              q1x_prev = 0
              q1y_prev = 0
              q_ind = 0

            end if

            !only advance in calculation if momentum changes
            if ( (map_q1%ix .NE. q1x_prev) .OR. (map_q1%iy .NE. q1y_prev) ) then
      
              q1x_prev = map_q1%ix
              q1y_prev = map_q1%iy
      
              q_ind = q_ind + 1
            
            end if
            ! ---

            DO idx_q2 = 1, Nb
              map_q2 = Index_Bosonic_IBZ(inode * Nb + idx_q2)
              if (map_q2%iw > 1) then

              !get momentum part only of index
              idx_q2_mom = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1

                !!!$omp parallel do collapse(2), &
                !!!$omp& private(l1, l2, l3, l4, idx_l1, idx_l2, idx_l3, idx_l4, map_l3, map_l4, f_temp, w1, w2, w3, w4)
                do l4 = 1, Nl
                  do l2 = 1, Nl
                    do l1 = 1, Nl
                    !do idx_l1 = 1, Nz
                    !  w1 = mod(idx_l1 - 1, Nl) + 1
                    !  l1 = int((idx_l1 - 1)/Nf) + 1
                      do l3 = 1, Nl
  
                        !calculate product of four formfactors
                        !f_temp = FFFF(map_q1%ix, map_q1%iy, map_q2%ix, map_q2%iy, &
                        !              l1, l2, l3, l4, 2)
 
                        !load product of four formactors from saved value
                        f_temp = FFFF_safe(l3, l1, l2, l4, idx_q2_mom, q_ind, 2)

                        !fermionic frequency
                        do w1 = 1, Nf

                          !map_l1 = Indxmap_L(l1, w1)
                          !idx_l1 = List_index_L(map_l1)

                          idx_l1 = (l1 - 1)*Nf + w1

                          !map_l2 = Indxmap_L(l2, w1 - map_q2%iw + 1) !nu_2 = nu_1 + w2
                          w2 = w1 - map_q2%iw + 1 !nu_2 = nu_1 - w2
                          IF (w2 >= 1 .AND. w2 <= Nf) THEN
                            idx_l2 = (l2 - 1)*Nf + w2 !List_Index_L(map_l2)
                            w4 = -(w1 + map_q1%iw - 1) + Nf + 1 !nu_4 = -nu_1 - w1
                            !map_l4 = Indxmap_L(l4, -(w1 + map_q1%iw -1) +Nf+1) !nu_4 = -nu_1 - w1
                            w3 = -w1 + Nf + 1
                            !idx_l3 = List_Index_L(map_l3)
  
                            IF (w4 >= 1 .AND. w4 <= Nf) THEN
  
                              idx_l4 = (l4 - 1)*Nf + w4 !List_Index_L(map_l4)
                              !map_l3 = Indxmap_L(l3, -w1 + Nf +1) !nu_3 = -nu_1
                              idx_l3 = (l3 - 1)*Nf + w3
  
                              SELECT CASE (ichannel)
  
                              CASE (1) !density
                                F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 signs(l3)*signs(l4)* & !for negative l-arguments
                                                                 (-0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                                 f_temp
                                F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 signs(l3)*signs(l4)* & !for negative l-arguments
                                                                 (-0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                                 f_temp
  
                              CASE (2) !magnetic
                                F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 signs(l3)*signs(l4)* & !for negative l-arguments
                                                                 (-1.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                                 f_temp
                                F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 signs(l3)*signs(l4)* & !for negative l-arguments
                                                                 (0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                                 f_temp
  
                              END SELECT
  
                              !put kernel function treatment here
                            ELSE
  
                              map_l3 = Indxmap_L(l3, w3)
                              map_l4 = Indxmap_L(l4, w4)
  
                              SELECT CASE (ichannel)
  
                              CASE (1) !density
                                F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 signs(l3)*signs(l4)* & !for negative l-arguments
                                                                 (-0.5d0)*CONJG(kernel('d', map_l3, map_l4, map_q2))* &
                                                                 f_temp
                                F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 signs(l3)*signs(l4)* & !for negative l-arguments
                                                                 (-0.5d0)*CONJG(kernel('d', map_l3, map_l4, map_q2))* &
                                                                 f_temp
  
                              CASE (2) !magnetic
                                F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 signs(l3)*signs(l4)* & !for negative l-arguments
                                                                 (-1.5d0)*CONJG(kernel('m', map_l3, map_l4, map_q2))* &
                                                                 f_temp
                                F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 signs(l3)*signs(l4)* & !for negative l-arguments
                                                                 (0.5d0)*CONJG(kernel('m', map_l3, map_l4, map_q2))* &
                                                                 f_temp
  
                              END SELECT
  
                            END IF
  
                          END IF
  
                        end do !w1

                      END DO !l4
                    end do !l3
                  end do !l2
                end do !l1
                !!!$omp end parallel do
              end if !time reversal
            END DO !idx_q2
          END DO !idx_q1
          ! END DO !idx_l1

!-----------------------------------------------------------------

          !1st contributions to singlet and triplet
          !rotation 3: Phi(k, k'; q-k-k') -> Phi(k1, k1'; q1)

          !these variables are arguments of LHS of equation

          ! DO idx_l1 = 1, Nz
          !   map_l1 = Index_L(idx_l1)
          DO idx_q1 = 1, Nb

            map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)    !make sure this is correct on any node

            ! --- some stuff to loop momenta independent of frequencies

            if(idx_q1 == 1) then

              !variables to loop momenta independent of frequencies
              q1x_prev = 0
              q1y_prev = 0
              q_ind = 0

            end if

            !only advance in calculation if momentum changes
            if ( (map_q1%ix .NE. q1x_prev) .OR. (map_q1%iy .NE. q1y_prev) ) then
      
              q1x_prev = map_q1%ix
              q1y_prev = map_q1%iy
      
              q_ind = q_ind + 1
            
            end if
            ! ---

            DO idx_q2 = 1, Nb
              map_q2 = Index_Bosonic_IBZ(inode * Nb + idx_q2)

              !get momentum part only of index
              idx_q2_mom = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1

                !!!$omp parallel do collapse(2), &
                !!!$omp& private(l1, l2, l3, l4, idx_l1, idx_l2, idx_l3, idx_l4, map_l3, map_l4, f_temp, w1, w2, w4)
                do l4 = 1, Nl
                  do l2 = 1, Nl
                    do l1 = 1, Nl
                    !do idx_l1 = 1, Nz
                    !  w1 = mod(idx_l1 - 1, Nl) + 1
                    !  l1 = int((idx_l1 - 1)/Nf) + 1
                      do l3 = 1, Nl

                        !calculate product of four formfactors
                        !f_temp = FFFF(map_q1%ix, map_q1%iy, map_q2%ix, map_q2%iy, &
                        !              l1, l2, l3, l4, 3)
 
                        !load product of four formactors from saved value
                        f_temp = FFFF_safe(l3, l1, l2, l4, idx_q2_mom, q_ind, 3)
 

                        !fermionic frequency
                        do w1 = 1, Nf

                          !map_l1 = Indxmap_L(l1, w1)

                          idx_l1 = (l1 - 1)*Nf + w1 !List_index_L(map_l1)

                          w2 = -(w1 + map_q2%iw - 1) + Nf + map_q1%iw !nu_2 = w1 - w2 - nu1
                          !map_l2 = Indxmap_L(l2, Add_MFBB( map_l1%iw, map_q1%iw,- map_q2%iw)) !nu_2 = w1 - w2 - nu1
                          IF (w2 >= 1 .AND. w2 <= Nf) THEN !check if inside frequency box
                            idx_l2 = (l2 - 1)*Nf + w2 !List_Index_L(map_l2)
  
                            !map_l4 = Indxmap_L(l4, Add_MFBB( map_l1%iw, map_q1%iw,- map_q2%iw)) !nu_4 = nu_2 = w1 - w2 - nu1
                            !IF(map_l4%iw >= 1 .AND. map_l4%iw <= Nf) THEN !check if inside frequency box
                            !idx_l4 = List_Index_L(map_l4)
                            idx_l4 = (l4 - 1)*Nf + w2 !nu_4 = nu_2 = w1 - w2 - nu1
                            !map_l3 = Indxmap_L(l3, w1) !nu_3 = nu_1
                            !idx_l3 = List_Index_L(map_l3)
                            idx_l3 = (l3 - 1)*Nf + w1 !nu_3 = nu_1
  
                            SELECT CASE (ichannel)
  
                            CASE (1) !density
                              F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                               (0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                               f_temp
                              F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                               (0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                               f_temp
  
                            CASE (2) !magnetic
                              F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                               (-1.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                               f_temp
                              F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                               (0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                               f_temp
  
                            END SELECT
  
                          END IF
  
                        end do !w1

                      END DO !l4
                    end do !l3
                  end do !l2
                end do !l1
                !!!$omp end parallel do
            END DO !idx_q2
          END DO !idx_q1
          ! END DO !idx_l1

         !Now time - reversal part

          DO idx_q1 = 1, Nb

            map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)    !make sure this is correct on any node

            ! --- some stuff to loop momenta independent of frequencies

            if(idx_q1 == 1) then

              !variables to loop momenta independent of frequencies
              q1x_prev = 0
              q1y_prev = 0
              q_ind = 0

            end if

            !only advance in calculation if momentum changes
            if ( (map_q1%ix .NE. q1x_prev) .OR. (map_q1%iy .NE. q1y_prev) ) then
      
              q1x_prev = map_q1%ix
              q1y_prev = map_q1%iy
      
              q_ind = q_ind + 1
            
            end if
            ! ---

            DO idx_q2 = 1, Nb
              map_q2 = Index_Bosonic_IBZ(inode * Nb + idx_q2)
              if (map_q2%iw > 1) then

              !get momentum part only of index
              idx_q2_mom = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1

                !!!$omp parallel do collapse(2), &
                !!!$omp& private(l1, l2, l3, l4, idx_l1, idx_l2, idx_l3, idx_l4, map_l3, map_l4, f_temp, w1, w2, w3, w4)
                do l4 = 1, Nl
                  do l2 = 1, Nl
                    do l1 = 1, Nl
                    !do idx_l1 = 1, Nz
                    !  w1 = mod(idx_l1 - 1, Nl) + 1
                    !  l1 = int((idx_l1 - 1)/Nf) + 1
                      do l3 = 1, Nl

                        !calculate product of four formfactors
                        !f_temp = FFFF(map_q1%ix, map_q1%iy, map_q2%ix, map_q2%iy, &
                        !              l1, l2, l3, l4, 4) 

                        !load product of four formactors from saved value
                        f_temp = FFFF_safe(l3, l1, l2, l4, idx_q2_mom, q_ind, 4)
 

                        !fermionic frequency
                        do w1 = 1, Nf

                          !map_l1 = Indxmap_L(l1, w1)

                          idx_l1 = (l1 - 1)*Nf + w1 !List_index_L(map_l1)

                          w2 = -w1 + Nf + 1 + map_q1%iw + map_q2%iw - 2 !nu_2 = w1 + w2 - nu1
                          !map_l2 = Indxmap_L(l2, Add_MFBB( w1, map_q1%iw, map_q2%iw)) !nu_2 = w1 + w2 - nu1
                          IF (w2 >= 1 .AND. w2 <= Nf) THEN
                            idx_l2 = (l2 - 1)*Nf + w2 !List_Index_L(map_l2)
  
                            !map_l4 = Indxmap_L(l4, Add_FBB(map_l1%iw, - map_q1%iw, - map_q2%iw))!nu_4 = nu_1 - w1 - w2
                            !IF(map_l4%iw >= 1 .AND. map_l4%iw <= Nf) THEN
                            w4 = w1 - map_q1%iw - map_q2%iw + 2 !nu_4 = -nu_2 = nu_1 - w1 - w2
                            !idx_l4 = List_Index_L(map_l4)
                            idx_l4 = (l4 - 1)*Nf + w4
                            !map_l3 = Indxmap_L(l3, Add_MFBB( map_l1%iw, 1, 1)) !nu3 = -nu1
                            w3 = -w1 + Nf + 1 !nu3 = -nu1
                            idx_l3 = (l3 - 1)*Nf + w3 !List_Index_L(map_l3) !nu_3 = -nu_1
  
                            SELECT CASE (ichannel)
  
                            CASE (1) !density
                              F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                               signs(l3)*signs(l4)* & !for negative l-arguments
                                                               (0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                               f_temp
                              F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                               signs(l3)*signs(l4)* & !for negative l-arguments
                                                               (0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                               f_temp
  
                            CASE (2) !magnetic
                              F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                               signs(l3)*signs(l4)* & !for negative l-arguments
                                                               (-1.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                               f_temp
                              F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                               signs(l3)*signs(l4)* & !for negative l-arguments
                                                               (0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                               f_temp
  
                            END SELECT
  
                          END IF
  
                        end do !w1

                      END DO !l4
                    end do !l3
                  end do !l2
                end do !l1
                !!!$omp end parallel do
              end if ! time reversal
            END DO !idx_q2
          END DO !idx_q1
          ! END DO !idx_l1

!---------------------------------------------------------------
!2nd contributions to singlet and triplet from density
!rotation 2: Phi(k, q-k'; k'-k) -> Phi(k1, k1'; q1)

          ! DO idx_l1 = 1, Nz
          !   map_l1 = Index_L(idx_l1)
          DO idx_q1 = 1, Nb

            map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)    !make sure this is correct on any node

            ! --- some stuff to loop momenta independent of frequencies

            if(idx_q1 == 1) then

              !variables to loop momenta independent of frequencies
              q1x_prev = 0
              q1y_prev = 0
              q_ind = 0

            end if

            !only advance in calculation if momentum changes
            if ( (map_q1%ix .NE. q1x_prev) .OR. (map_q1%iy .NE. q1y_prev) ) then
      
              q1x_prev = map_q1%ix
              q1y_prev = map_q1%iy
      
              q_ind = q_ind + 1
            
            end if
            ! ---

            DO idx_q2 = 1, Nb
              map_q2 = Index_Bosonic_IBZ(inode * Nb + idx_q2)

              !get momentum part only of index
              idx_q2_mom = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1

                !!!$omp parallel do collapse(2), &
                !!!$omp& private(l1, l2, l3, l4, idx_l1, idx_l2, idx_l3, idx_l4, map_l3, map_l4, f_temp, w1, w2, w4)
                do l4 = 1, Nl
                  do l2 = 1, Nl
                    do l1 = 1, Nl
                    !do idx_l1 = 1, Nz
                    !  w1 = mod(idx_l1 - 1, Nl) + 1
                    !  l1 = int((idx_l1 - 1)/Nf) + 1
                      do l3 = 1, Nl

                        !calculate product of four formfactors
                        !f_temp = FFFF(map_q1%ix, map_q1%iy, map_q2%ix, map_q2%iy, &
                        !              l1, l2, l3, l4, 5) 

                        !load product of four formactors from saved value
                        f_temp = FFFF_safe(l3, l1, l2, l4, idx_q2_mom, q_ind, 5)

 
                        !fermionic frequency
                        do w1 = 1, Nf
                          idx_l1 = (l1 - 1)*Nf + w1
                          !map_l2 = Indxmap_L(l2, Add_FBB(map_l1%iw, map_q2%iw, 1)) !nu_2 = nu_1 + w2
                          w2 = w1 + map_q2%iw - 1 !nu_2 = nu_1 + w2
  
                          IF (w2 >= 1 .AND. w2 <= Nf) THEN !check if inside frequency box
                            !idx_l2 = List_Index_L(map_l2)
                            idx_l2 = (l2 - 1)*Nf + w2
                            !map_l4 = Indxmap_L(l4, Add_MFBB( map_l1%iw, map_q1%iw, - map_q2%iw)) !nu_4 = w1 - w2 - nu1
                            w4 = -(w1 + map_q2%iw - 1) + Nf + map_q1%iw !nu_4 = w1 - w2 - nu1
  
                            IF (w4 >= 1 .AND. w4 <= Nf) THEN !check if inside frequency box
  
                              idx_l4 = (l4 - 1)*Nf + w4 !List_Index_L(map_l4)
                              !map_l3 = Indxmap_L(l3, map_l1%iw) !nu_3 = nu_1
                              idx_l3 = (l3 - 1)*Nf + w1 !List_Index_L(map_l3) !nu_3 = nu_1
  
                              SELECT CASE (ichannel)
  
                              CASE (1) !density
                                F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                                 f_temp
                                F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                                 f_temp
  
                              CASE (2) !magnetic
                                F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-1.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                                 f_temp
                                F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                                 f_temp
  
                              END SELECT
  
                              !later put kernels here
                            ELSE
  
                              map_l3 = Indxmap_L(l3, w1)
                              map_l4 = Indxmap_L(l4, w4)
  
                              SELECT CASE (ichannel)
  
                              CASE (1) !density
                                F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (0.5d0)*kernel('d', map_l3, map_l4, map_q2)* &
                                                                 f_temp
                                F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-0.5d0)*kernel('d', map_l3, map_l4, map_q2)* &
                                                                 f_temp
  
                              CASE (2) !magnetic
                                F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-1.5d0)*kernel('m', map_l3, map_l4, map_q2)* &
                                                                 f_temp
                                F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                                 (-0.5d0)*kernel('m', map_l3, map_l4, map_q2)* &
                                                                 f_temp
  
                              END SELECT
  
                            END IF
                          END IF
  
                        end do !w1
                      END DO !l4
                    end do !l3
                  end do !l2
                end do !l1
                !!!$omp end parallel do
            END DO !idx_q2
          END DO !idx_q1
          ! END DO !idx_l1

          !Now time - reversal part

          DO idx_q1 = 1, Nb

            map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)    !make sure this is correct on any node

            ! --- some stuff to loop momenta independent of frequencies

            if(idx_q1 == 1) then

              !variables to loop momenta independent of frequencies
              q1x_prev = 0
              q1y_prev = 0
              q_ind = 0

            end if

            !only advance in calculation if momentum changes
            if ( (map_q1%ix .NE. q1x_prev) .OR. (map_q1%iy .NE. q1y_prev) ) then
      
              q1x_prev = map_q1%ix
              q1y_prev = map_q1%iy
      
              q_ind = q_ind + 1
            
            end if
            ! ---

            DO idx_q2 = 1, Nb
              map_q2 = Index_Bosonic_IBZ(inode * Nb + idx_q2)
              if (map_q2%iw > 1) then

              !get momentum part only of index
              idx_q2_mom = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1

                  !!!$omp parallel do collapse(2), &
                  !!!$omp& private(l1, l2, l3, l4, idx_l1, idx_l2, idx_l3, idx_l4, map_l3, map_l4, f_temp, w1, w2, w3, w4)
                  do l4 = 1, Nl
                    do l2 = 1, Nl
                      do l1 = 1, Nl
                      !do idx_l1 = 1, Nz
                      !  w1 = mod(idx_l1 - 1, Nl) + 1
                      !  l1 = int((idx_l1 - 1)/Nf) + 1
                        do l3 = 1, Nl
 
                        !calculate product of four formfactors
                        !f_temp = FFFF(map_q1%ix, map_q1%iy, map_q2%ix, map_q2%iy, &
                        !              l1, l2, l3, l4, 6) 
 
                        !load product of four formactors from saved value
                        f_temp = FFFF_safe(l3, l1, l2, l4, idx_q2_mom, q_ind, 6)
 
  
                          !fermionic frequency
                          do w1 = 1, Nf
                            idx_l1 = (l1 - 1)*Nf + w1
                            w2 = w1 - map_q2%iw + 1 !nu_2 = nu1 - w2
                            !map_l2 = Indxmap_L(l2, Add_FBB(map_l1%iw, - map_q2%iw, 1)) !nu_2 = nu1 - w2
                            IF (w2 >= 1 .AND. w2 <= Nf) THEN
                              idx_l2 = (l2 - 1)*Nf + w2 !List_Index_L(map_l2)
                              !map_l4 = Indxmap_L(l4, Add_FBB(map_l1%iw, - map_q1%iw, - map_q2%iw))!nu_4 = nu_1 - w1 - w2
                              w4 = w1 - map_q1%iw - map_q2%iw + 2 !nu_4 = nu_1 - w1 - w2
                              w3 = -w1 + Nf + 1 !nu3 = -nu1
  
                              IF (w4 >= 1 .AND. w4 <= Nf) THEN
  
                                idx_l4 = (l4 - 1)*Nf + w4 !List_Index_L(map_l4)
                                !map_l3 = Indxmap_L(l3, Add_MFBB( map_l1%iw, 1, 1)) !nu3 = -nu1
                                !w3 = -w1 + Nf + 1 !nu3 = -nu1
                                idx_l3 = (l3 - 1)*Nf + w3 !List_Index_L(map_l3) !nu_3 = -nu_1
  
                                SELECT CASE (ichannel)
  
                                CASE (1) !density
                                  F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                                   signs(l3)*signs(l4)* & !for negative l-arguments
                                                                   (0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                                   f_temp
                                  F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                                   signs(l3)*signs(l4)* & !for negative l-arguments
                                                                   (-0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                                   f_temp
  
                                CASE (2) !magnetic
                                  F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                                   signs(l3)*signs(l4)* & !for negative l-arguments
                                                                   (-1.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                                   f_temp
                                  F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                                   signs(l3)*signs(l4)* & !for negative l-arguments
                                                                   (-0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                                   f_temp
  
                                END SELECT
  
                                !put kernel function treatment here
                              ELSE
  
                                map_l3 = Indxmap_L(l3, w3)
                                map_l4 = Indxmap_L(l4, w4)
  
                                SELECT CASE (ichannel)
  
                                CASE (1) !density
                                  F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                                   signs(l3)*signs(l4)* & !for negative l-arguments
                                                                   (0.5d0)*CONJG(kernel('d', map_l3, map_l4, map_q2))* &
                                                                   f_temp
                                  F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                                   signs(l3)*signs(l4)* & !for negative l-arguments
                                                                   (-0.5d0)*CONJG(kernel('d', map_l3, map_l4, map_q2))* &
                                                                   f_temp
  
                                CASE (2) !magnetic
                                  F_s_LL(idx_l1, idx_l2, idx_q1) = F_s_LL(idx_l1, idx_l2, idx_q1) + &
                                                                   signs(l3)*signs(l4)* & !for negative l-arguments
                                                                   (-1.5d0)*CONJG(kernel('m', map_l3, map_l4, map_q2))* &
                                                                   f_temp
                                  F_t_LL(idx_l1, idx_l2, idx_q1) = F_t_LL(idx_l1, idx_l2, idx_q1) + &
                                                                   signs(l3)*signs(l4)* & !for negative l-arguments
                                                                   (-0.5d0)*CONJG(kernel('m', map_l3, map_l4, map_q2))* &
                                                                   f_temp
  
                                END SELECT
  
                              END IF
                            END IF
  
                          end do !w1
                        END DO !l4
                      end do !l3
                    end do !l2
                  end do !l1
                  !!!$omp end parallel do
                end if ! time reversal
            END DO !idx_q2
          END DO !idx_q1
          ! END DO !idx_l1

          !contributions from singlet & triplet
        CASE (3, 4)

          !rotation 2: Phi(k, q-k'; k'-k) -> Phi(k1, k1'; q1)

          !these variables are arguments of LHS of equation
          ! DO idx_l1 = 1, Nz
          !  map_l1 = Index_L(idx_l1)
          DO idx_q1 = 1, Nb

            map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)    !make sure this is correct on any node

            ! --- some stuff to loop momenta independent of frequencies

            if(idx_q1 == 1) then

              !variables to loop momenta independent of frequencies
              q1x_prev = 0
              q1y_prev = 0
              q_ind = 0

            end if

            !only advance in calculation if momentum changes
            if ( (map_q1%ix .NE. q1x_prev) .OR. (map_q1%iy .NE. q1y_prev) ) then
      
              q1x_prev = map_q1%ix
              q1y_prev = map_q1%iy
      
              q_ind = q_ind + 1
            
            end if
            ! ---

            DO idx_q2 = 1, Nb
              map_q2 = Index_Bosonic_IBZ(inode * Nb + idx_q2)

              !get momentum part only of index
              idx_q2_mom = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1

                !!!$omp parallel do collapse(2), &
                !!!$omp& private(l1, l2, l3, l4, idx_l1, idx_l2, idx_l3, idx_l4, map_l3, map_l4, f_temp, w1, w2, w4)
                do l4 = 1, Nl
                  do l2 = 1, Nl
                    do l1 = 1, Nl
                    !do idx_l1 = 1, Nz
                    !  w1 = mod(idx_l1 - 1, Nl) + 1
                    !  l1 = int((idx_l1 - 1)/Nf) + 1
                      do l3 = 1, Nl

                        !calculate product of four formfactors
                        !f_temp = FFFF(map_q1%ix, map_q1%iy, map_q2%ix, map_q2%iy, &
                        !              l1, l2, l3, l4, 7) 

                        !load product of four formactors from saved value
                        f_temp = FFFF_safe(l3, l1, l2, l4, idx_q2_mom, q_ind, 7)

 
                        !fermionic frequency
                        do w1 = 1, Nf
                          !map_l1 = Indxmap_L(l1, w1)
                          idx_l1 = (l1 - 1)*Nf + w1 !List_index_L(map_l1)
                          w2 = -(w1 + map_q1%iw - 1) + Nf + map_q2%iw !nu_2 = w2 - w1 - nu1
                          !map_l2 = Indxmap_L(l2, w_aux) !nu_2 = w2 - w1 - nu1
                          IF (w2 >= 1 .AND. w2 <= Nf) THEN !check if inside frequency box
                            idx_l2 = (l2 - 1)*Nf + w2 ! List_Index_L(map_l2)
  
                            !map_l4 = Indxmap_L(l4, w2) !nu_4 = w2 - w1 - nu1
                            !IF(w_aux >= 1 .AND. w_aux <= Nf) THEN !check if inside frequency box
                            idx_l4 = (l4 - 1)*Nf + w2 !nu4 = nu2 !List_Index_L(map_l4)
                            !map_l3 = Indxmap_L(l3, w1) !nu_3 = nu_1
                            idx_l3 = (l3 - 1)*Nf + w1 !!nu_3 = nu_1 !List_Index_L(map_l3)
  
                            SELECT CASE (ichannel)
  
                            CASE (3) !singlet
                              F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                               0.5d0*mat(idx_l3, idx_l4, idx_q2)* &
                                                               f_temp
                              F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                               (-0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                               f_temp
  
                            CASE (4) !triplet
                              F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                               (1.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                               f_temp
                              F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                               (0.5d0)*mat(idx_l3, idx_l4, idx_q2)* &
                                                               f_temp
  
                            END SELECT
  
                            !END IF
                          END IF
  
                        end do !w1
                      END DO !l4
                    end do !l3
                  end do !l2
                end do !l1
                !!!$omp end parallel do
              !end do !s
            END DO !idx_q2
          END DO !idx_q1
          ! END DO !idx_l1
!
          !Now time - reversal part

          !these variables are arguments of LHS of equation

          !DO idx_l1 = 1, Nz
          ! map_l1 = Index_L(idx_l1)
          DO idx_q1 = 1, Nb

            map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)    !make sure this is correct on any node

            ! --- some stuff to loop momenta independent of frequencies

            if(idx_q1 == 1) then

              !variables to loop momenta independent of frequencies
              q1x_prev = 0
              q1y_prev = 0
              q_ind = 0

            end if

            !only advance in calculation if momentum changes
            if ( (map_q1%ix .NE. q1x_prev) .OR. (map_q1%iy .NE. q1y_prev) ) then
      
              q1x_prev = map_q1%ix
              q1y_prev = map_q1%iy
      
              q_ind = q_ind + 1
            
            end if
            ! ---

            DO idx_q2 = 1, Nb
              map_q2 = Index_Bosonic_IBZ(inode * Nb + idx_q2)
              if (map_q2%iw > 1) then

              !get momentum part only of index
              idx_q2_mom = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1

                !!!$omp parallel do collapse(2), &
                !!!$omp& private(l1, l2, l3, l4, idx_l1, idx_l2, idx_l3, idx_l4, map_l3, map_l4, f_temp, w1, w2, w3, w4)
                do l4 = 1, Nl
                  do l2 = 1, Nl
                    do l1 = 1, Nl
                    !do idx_l1 = 1, Nz
                    !  w1 = mod(idx_l1 - 1, Nl) + 1
                    !  l1 = int((idx_l1 - 1)/Nf) + 1
                      do l3 = 1, Nl

                        !calculate product of four formfactors
                        !f_temp = FFFF(map_q1%ix, map_q1%iy, map_q2%ix, map_q2%iy, &
                        !              l1, l2, l3, l4, 8) 

                        !load product of four formactors from saved value
                        f_temp = FFFF_safe(l3, l1, l2, l4, idx_q2_mom, q_ind, 8)

 
                        !fermionic frequency
                        do w1 = 1, Nf
                          !map_l1 = Indxmap_L(l1, w1)
                          !idx_l1 = List_index_L(map_l1)
                          idx_l1 = (l1 - 1)*Nf + w1
                          w4 = w1 + map_q1%iw + map_q2%iw - 2 !nu4 = nu1 + w1 + w2
                          !map_l2 = Indxmap_L(l2, -w_aux + Nf + 1) !nu_2 = -nu4 = -nu1 - w1 - w2
                          w2 = -w4 + Nf + 1 !nu_2 = -nu4 = -nu1 - w1 - w2
                          IF (w2 >= 1 .AND. w2 <= Nf) THEN
                            !idx_l2 = List_Index_L(map_l2)
                            idx_l2 = (l2 - 1)*Nf + w2
                            !map_l4 = Indxmap_L(l4, w_aux)!nu4 = nu1 + w1 + w2
                            !IF(map_l4%iw >= 1 .AND. map_l4%iw <= Nf) THEN
                            !idx_l4 = List_Index_L(map_l4)
                            idx_l4 = (l4 - 1)*Nf + w4
                            !map_l3 = Indxmap_L(l3, -w1 + Nf +1) !nu3 = -nu1
                            w3 = -w1 + Nf + 1 !nu3 = -nu1
                            idx_l3 = (l3 - 1)*Nf + w3 !List_Index_L(map_l3) !nu_3 = -nu_1
  
                            SELECT CASE (ichannel)
  
                            CASE (3) !singlet
                              F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                               signs(l3)*signs(l4)* & !for negative l-arguments
                                                               (0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                               f_temp
                              F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                               signs(l3)*signs(l4)* & !for negative l-arguments
                                                               (-0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                               f_temp
  
                            CASE (4) !triplet
                              F_d_LL(idx_l1, idx_l2, idx_q1) = F_d_LL(idx_l1, idx_l2, idx_q1) + &
                                                               signs(l3)*signs(l4)* & !for negative l-arguments
                                                               (1.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                               f_temp
                              F_m_LL(idx_l1, idx_l2, idx_q1) = F_m_LL(idx_l1, idx_l2, idx_q1) + &
                                                               signs(l3)*signs(l4)* & !for negative l-arguments
                                                               (0.5d0)*CONJG(mat(idx_l3, idx_l4, idx_q2))* &
                                                               f_temp
  
                            END SELECT
  
                          END IF
  
                        end do !w1
                      END DO !l4
                    end do !l3
                  end do !l2
                end do !l1
                !!!$omp end parallel do
              end if !time reversal
            END DO !idx_q2
          END DO !idx_q1
          ! END DO !idx_l1

        !select of channel contributions
        END SELECT

!end of inode loop
      END DO
!end of ichannel loop
    END DO

  end subroutine solve_parquet_equation

! ------------------------

  !function returns product of four formfactors
  real(dp) function FFFF(q1x, q1y, q2x, q2y, l1, l2, l3, l4, cse) result(res) 
  
    !x and y component of bosonic momenta
    integer, intent(in) :: q1x, q1y, q2x, q2y
    !FF indices
    integer, intent(in) :: l1, l2, l3, l4 
    !place at which function is called
    integer, intent(in) :: cse

    !symmetry for IBZ treatment
    integer :: is
    !to save and symmetrize q2 arguemnts
    type(indxmap) map_q2, map_q2_sym
    !for getting the whole BZ from IBZ
    logical, dimension(8) :: symm_list

    !inner loop variables
    integer :: kx, ky

    res = 0.0d0

    !map_q2_sym = indxmap(q2x, q2y, 0)

    map_q2 = indxmap(q2x, q2y, 0)

    !get all symmetries that are needed to fill whole BZ
    call list_symmetries(map_q2, symm_list)                  

    !loop over symmetries
    do is = 1, Ns

      !leave out the ones that we already took
      IF(.NOT. symm_list(is)) CYCLE
      !create symmetrized q2 to index formfactors
      call symmetry_operation(is, map_q2, map_q2_sym)

      select case(cse)

        case(1)
  
          do kx =1, Nx
            do ky = 1,Ny
  
              res = res + FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                       FF(mod(kx + map_q2_sym%ix - 2, Nx) + 1, &
                       mod(ky + map_q2_sym%iy -2, Ny) + 1, l2) * &
                       FF_inv(is, mod(kx + q1x - 2, Nx) + 1, &
                       mod(ky + q1y -2, Ny) + 1, l4)
            END DO
          END DO
  
        case(2)
 
          do kx =1, Nx
            do ky = 1,Ny
          
              res = res + &
                       FF(kx, ky, l1) * &
                       FF(mod(kx - map_q2_sym%ix + Nx, Nx) + 1, &
                       mod(ky - map_q2_sym%iy + Ny, Ny) + 1, l2) * &
                       FF_inv(is, kx, ky, l3) * &
                       FF_inv(is, mod(kx + q1x -2, Nx) + 1, &
                       mod(ky + q1y - 2, Ny) + 1, l4)
            END DO
          END DO

        case(3)

          do kx =1, Nx
            do ky = 1,Ny

              res = res + FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                       FF(mod(q1x - map_q2_sym%ix - kx + 2 * Nx + 1, Nx) + 1, &
                       mod(q1y - map_q2_sym%iy - ky + 2 * Ny + 1, Ny) + 1, l2) * &
                       FF_inv(is, mod(q1x - map_q2_sym%ix - kx + 2 * Nx + 1, Nx) + 1, &
                       mod(q1y - map_q2_sym%iy - ky + 2 * Ny + 1, Ny) + 1, l4)
            END DO
          END DO
  
        case(4)

          do kx =1, Nx
            do ky = 1,Ny

               res = res + FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                        FF(mod( q1x + map_q2_sym%ix - kx + Nx -1, Nx) + 1, &
                        mod( q1y + map_q2_sym%iy - ky + Ny - 1, Ny) + 1, l2) * &
                        FF_inv(is, mod( q1x + map_q2_sym%ix - kx + Nx -1, Nx) + 1, &
                        mod( q1y + map_q2_sym%iy - ky +  Ny -1, Ny) + 1, l4)
            END DO
          END DO
  
        case(5)

          do kx =1, Nx
            do ky = 1,Ny

              res = res + FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                       FF(mod(map_q2_sym%ix + kx - 2, Nx) + 1, &
                       mod(map_q2_sym%iy + ky - 2, Ny) + 1, l2) * &
                       FF_inv(is, mod(q1x - map_q2_sym%ix - kx + 2 * Nx + 1, Nx) + 1, &
                       mod(q1y - map_q2_sym%iy - ky + 2 * Ny + 1, Ny) + 1, l4)
            END DO
          END DO
  
        case(6)

          do kx =1, Nx
            do ky = 1,Ny

              res = res + FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                       FF(mod( - map_q2_sym%ix + kx + Nx + 0, Nx) + 1, &
                       mod( - map_q2_sym%iy + ky + Ny + 0, Ny) + 1, l2) * &
                       FF_inv(is, mod(q1x + map_q2_sym%ix - kx + Nx - 1, Nx) + 1, &
                       mod(q1y + map_q2_sym%iy - ky + Ny - 1, Ny) + 1, l4)
            END DO
          END DO
      
        case(7)

          do kx =1, Nx
            do ky = 1,Ny

              res = res + FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                       FF(mod(map_q2_sym%ix - q1x - kx + 1 + 2 * Nx, Nx) + 1, &
                       mod(map_q2_sym%iy - q1y - ky + 1 + 2 * Ny, Ny) + 1, l2) * &
                       FF_inv(is, mod(map_q2_sym%ix - q1x - kx + 1 + 2 * Nx, Nx) + 1, &
                       mod(map_q2_sym%iy - q1y - ky + 1 + 2 * Ny, Ny) + 1, l4)
            END DO
          END DO
  
        case(8)

          do kx =1, Nx
            do ky = 1,Ny

              res = res + FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                       FF(mod(- map_q2_sym%ix - q1x - kx + 3 * Nx + 3, Nx) + 1, &
                       mod(- map_q2_sym%iy - q1y - ky + 3 * Ny + 3, Ny) + 1, l2) * &
                       FF_inv(is, mod(- map_q2_sym%ix - q1x - kx + 3 * Nx + 3, Nx) + 1, & 
                       mod(- map_q2_sym%iy - q1y - ky + 3 * Ny + 3, Ny) + 1, l4)
            END DO
          END DO

      end select !cse

    end do !is

  end function FFFF

! --------


  !initializes the array FFFF_safe which holds the product of four formfactors
  subroutine FFFF_mem
 
  !$use omp_lib
   
    !x and y component of bosonic momenta
    integer :: q_ind, q1x_prev, q1y_prev, idx_q1, idx_q2
    integer :: q2x, q2y
    !map q1 for determining correct arguments on each node
    type(indxmap) :: map_q1
    !FF indices
    integer :: l1, l2, l3, l4 
    !symmetry for IBZ treatment
    integer :: is

    integer :: kx, ky

    !to save and symmetrize q2 arguemnts
    type(indxmap) map_q2, map_q2_sym
    !for getting the whole BZ from IBZ
    logical, dimension(8) :: symm_list


    q1x_prev = 0
    q1y_prev = 0
    q_ind = 0

    !argument for q1 saves less than whole IBZ but a little more than necessary
    if(.NOT. allocated(FFFF_safe)) then

      !Nb + Nf - 1 since if Nb = Nf this should be 2 
      !and if Nb = Nf + 1 this should be 3 ... and so on ...
      allocate(FFFF_safe( Nl, Nl, Nl, Nl, Nred * 2/Nf , int((Nb + Nf/2 - 1)/(Nf/2)) + 1, 8 ))

    end if

    FFFF_safe = 0.0d0

    do idx_q1 = 1, Nb
      !to get correct arguments on any node
      map_q1 = Index_Bosonic_IBZ(id * Nb + idx_q1)
      !only advance in calculation if momentum changes
      if (map_q1%ix == q1x_prev .AND. map_q1%iy == q1y_prev) CYCLE

      q1x_prev = map_q1%ix
      q1y_prev = map_q1%iy

      q_ind = q_ind + 1 


!      !$omp parallel do !!loop variables are by default private
      !do q2x = 1, Nx
      !  do q2y = 1, Ny
      do q2x = 1, Nx_IBZ !loop q2 only over IBZ
        do q2y = 1, q2x
          do is = 1, Ns

          map_q2 = indxmap(q2x, q2y, 1)

          idx_q2 = (List_Index_IBZ(map_q2) - 1) / (Nf/2) + 1

          !get all symmetries that are needed to fill whole BZ
          call list_symmetries(map_q2, symm_list)                  

          !loop over symmetries
          !do is = 1, Ns

          !leave out the ones that we already took
          IF(.NOT. symm_list(is)) CYCLE
          !create symmetrized q2 to index formfactors
          call symmetry_operation(is, map_q2, map_q2_sym)


            do l4 = 1, Nl
              do l2 = 1, Nl
                do l1 = 1, Nl
                  do l3 = 1, Nl

                    do kx =1, Nx
                      do ky = 1, Ny


    FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 1) = &
    FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 1) + &
    FF(kx, ky, l1) * &
    FF_inv(is, kx, ky, l3) * &
    FF(mod(kx + map_q2_sym%ix - 2, Nx) + 1, mod(ky + map_q2_sym%iy - 2, Ny) + 1, l2) * &
    FF_inv(is, mod(kx + map_q1%ix - 2, Nx) + 1, mod(ky + map_q1%iy - 2, Ny) + 1, l4)


                        !FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 1) = &
                        !FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 1) + &
                        !FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                        !FF(mod(kx + map_q2_sym%ix - 2, Nx) + 1, &
                        !mod(ky + map_q2_sym%iy -2, Ny) + 1, l2) * &
                        !FF_inv(is, mod(kx + map_q1%ix - 2, Nx) + 1, &
                        !mod(ky + map_q1%iy -2, Ny) + 1, l4)
                   
 
                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 2) = &
                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 2) + &
                        FF(kx, ky, l1) * &
                        FF(mod(kx - map_q2_sym%ix + Nx, Nx) + 1, &
                        mod(ky - map_q2_sym%iy + Ny, Ny) + 1, l2) * &
                        FF_inv(is, kx, ky, l3) * &
                        FF_inv(is, mod(kx + map_q1%ix -2, Nx) + 1, &
                        mod(ky + map_q1%iy - 2, Ny) + 1, l4)

                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 3) = &
                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 3) + &
                        FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                        FF(mod(map_q1%ix - map_q2_sym%ix - kx + 2 * Nx + 1, Nx) + 1, &
                        mod(map_q1%iy - map_q2_sym%iy - ky + 2 * Ny + 1, Ny) + 1, l2) * &
                        FF_inv(is, mod(map_q1%ix - map_q2_sym%ix - kx + 2 * Nx + 1, Nx) + 1, &
                        mod(map_q1%iy - map_q2_sym%iy - ky + 2 * Ny + 1, Ny) + 1, l4)

                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 4) = &
                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 4) + &
                        FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                        FF(mod( map_q1%ix + map_q2_sym%ix - kx + Nx -1, Nx) + 1, &
                        mod( map_q1%iy + map_q2_sym%iy - ky + Ny - 1, Ny) + 1, l2) * &
                        FF_inv(is, mod( map_q1%ix + map_q2_sym%ix - kx + Nx -1, Nx) + 1, &
                        mod( map_q1%iy + map_q2_sym%iy - ky +  Ny -1, Ny) + 1, l4)

                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 5) = &
                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 5) + &
                        FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                        FF(mod(map_q2_sym%ix + kx - 2, Nx) + 1, &
                        mod(map_q2_sym%iy + ky - 2, Ny) + 1, l2) * &
                        FF_inv(is, mod(map_q1%ix - map_q2_sym%ix - kx + 2 * Nx + 1, Nx) + 1, &
                        mod(map_q1%iy - map_q2_sym%iy - ky + 2 * Ny + 1, Ny) + 1, l4)

                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 6) = &
                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 6) + &
                        FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                        FF(mod( - map_q2_sym%ix + kx + Nx + 0, Nx) + 1, &
                        mod( - map_q2_sym%iy + ky + Ny + 0, Ny) + 1, l2) * &
                        FF_inv(is, mod(map_q1%ix + map_q2_sym%ix - kx + Nx - 1, Nx) + 1, &
                        mod(map_q1%iy + map_q2_sym%iy - ky + Ny - 1, Ny) + 1, l4)

                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 7) = &
                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 7) + &
                        FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                        FF(mod(map_q2_sym%ix - map_q1%ix - kx + 1 + 2 * Nx, Nx) + 1, &
                        mod(map_q2_sym%iy - map_q1%iy - ky + 1 + 2 * Ny, Ny) + 1, l2) * &
                        FF_inv(is, mod(map_q2_sym%ix - map_q1%ix - kx + 1 + 2 * Nx, Nx) + 1, &
                        mod(map_q2_sym%iy - map_q1%iy - ky + 1 + 2 * Ny, Ny) + 1, l4)

                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 8) = &
                        FFFF_safe(l3, l1, l2, l4, idx_q2, q_ind, 8) + &
                        FF(kx, ky, l1) * FF_inv(is, kx, ky, l3) * &
                        FF(mod(- map_q2_sym%ix - map_q1%ix - kx + 3 * Nx + 3, Nx) + 1, &
                        mod(- map_q2_sym%iy - map_q1%iy - ky + 3 * Ny + 3, Ny) + 1, l2) * &
                        FF_inv(is, mod(- map_q2_sym%ix - map_q1%ix - kx + 3 * Nx + 3, Nx) + 1, & 
                        mod(- map_q2_sym%iy - map_q1%iy - ky + 3 * Ny + 3, Ny) + 1, l4)


                      END DO !kx
                    END DO !ky

                  end do !l3
                end do !l1
              end do !l2
            end do !l4
          end do !is
        end do !q2y
      end do !q2x
!      !$omp end parallel do
    end do !idx_q1


  end subroutine FFFF_mem


end module parquet_equation

