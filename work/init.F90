module init_module


    use amrex_fort_module, only : amrex_spacedim, amrex_real
    
      implicit none
    
      public
    
    contains
    
    subroutine init_phi(lo, hi, mf, mf_lo, mf_hi, prob_lo, prob_hi, dx, itf_position) &
        bind(C, name="init_phi")
        integer,            intent(in   ) :: lo(3),hi(3),mf_lo(3),mf_hi(3)
        real(amrex_real),   intent(inout) :: mf(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3))
        real(amrex_real),   intent(in)    :: prob_lo(3), prob_hi(3), itf_position
        real(amrex_real),   intent(in   ) :: dx(3) 


        ! local 
        integer :: i, j, k
        real    :: x, y, z, y_height, tep_position
        real    :: width
        real(amrex_real) noise(mf_lo(1): mf_hi(1), mf_lo(2): mf_hi(2), mf_lo(3): mf_hi(3))
        
        
        ! compute flux locally 
        ! call random_seed()
        call random_number(noise)
        width = 0.2
        y_height = prob_hi(2) - prob_lo(2)
        mf = -1
        do k = mf_lo(3), mf_hi(3)
            z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
            do j = mf_lo(2), mf_hi(2)
                y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
                do i = mf_lo(1), mf_hi(1)
                    x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

                    mf(i, j, k) = 0.5 * ( tanh((itf_position - x + 1.5 * noise(i, j, k)) * 1.5) + 1)
                    ! mf(i, j, k) = 0.5 * ( tanh((itf_position - x + 1.5 * noise(i, j, k) + 5.0 * sin(2 * (y + z) * 3.14159 / 20.0) ) * 1.5) + 1)

                    ! if (x .le. 10) then 
                    !     mf(i, j, k) = 1.0
                    ! else if (x .le. 25) then
                    !     if (abs(y - 50) - 5 .lt. 5 / 7.5 / 7.5 * abs(x - 17.5) * abs(x - 17.5)) then
                    !         mf(i, j, k) = 1
                    !     else
                    !         mf(i, j, k) = 0
                    !     end if
                    ! else if (x .le. 65) then
                    !     if (y .gt. 35 .and. y .lt. 65) then
                    !         mf(i, j, k) = 1
                    !     else
                    !         mf(i, j, k) = 0
                    !     end if
                    ! else
                    !     mf(i, j, k) = 0
                    ! end if

                end do  ! end i
            end do  ! end j
        end do  ! end k
        
    end subroutine init_phi
      

    subroutine init_bcoef(lo, hi, mf, mf_lo, mf_hi, phi,phi_lo, phi_hi, prob_lo, prob_hi, dx, cond_liq, cond_sld) &
        bind(C, name="init_bcoef")
        integer,            intent(in   ) :: lo(3),hi(3),mf_lo(3),mf_hi(3), phi_lo(3),phi_hi(3)
        real(amrex_real),   intent(inout) :: mf(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3))
        real(amrex_real),   intent(in   ) :: phi(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3))

        real(amrex_real),   intent(in)    :: prob_lo(3), prob_hi(3)
        real(amrex_real),   intent(in   ) :: dx(3), cond_liq, cond_sld


        ! local 
        integer i, j, k
        real(amrex_real)  x, y, z
        real(amrex_real)  interpolation ! h(x) = x^3 (6 x^2 - 15x + 10)

        do k = mf_lo(3), mf_hi(3)
            z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
            do j = mf_lo(2), mf_hi(2)
                y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
                do i = mf_lo(1), mf_hi(1)
                    x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

                    interpolation = (phi(i, j, k) ** 3) * (6 * phi(i, j, k) ** 2 - 15 * phi(i, j, k) + 10)
                    mf(i, j, k) = interpolation * cond_sld + (1 - interpolation) * cond_liq

                end do  ! end i
            end do  ! end j
        end do  ! end k

    end subroutine init_bcoef

    subroutine init_potential(lo, hi, mf, mf_lo, mf_hi, phi, phi_lo, phi_hi, prob_lo, prob_hi, dx, voltage, itf_position) &
        bind(C, name="init_potential")
        integer,            intent(in   ) :: lo(3),hi(3),mf_lo(3),mf_hi(3),phi_lo(3),phi_hi(3)
        real(amrex_real),   intent(inout) :: mf(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3))
        real(amrex_real),   intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))

        real(amrex_real),   intent(in)    :: prob_lo(3), prob_hi(3)
        real(amrex_real),   intent(in   ) :: dx(3), voltage, itf_position

        ! local 
        integer :: i, j, k
        real    :: x, y, z

        do k = lo(3), hi(3)
            z = prob_lo(3) + (dble(k)+0.5d0) * dx(3)
            do j = lo(2), hi(2)
                y = prob_lo(2) + (dble(j)+0.5d0) * dx(2)
                do i = lo(1), hi(1)
                    x = prob_lo(1) + (dble(i)+0.5d0) * dx(1)

                end do  ! end i
            end do  ! end j
        end do  ! end k

    end subroutine init_potential

    subroutine update_solute(lo, hi, mf, mf_lo, mf_hi, phi, phi_lo, phi_hi, mu, mu_lo, mu_hi, prob_lo, prob_hi, dx, &
        epsilon_liq, epsilon_sld &
        )&
        bind(C, name="update_solute") 

        integer,            intent(in   ) :: lo(3),hi(3),mf_lo(3),mf_hi(3),phi_lo(3),phi_hi(3),mu_lo(3),mu_hi(3)
        real(amrex_real),   intent(inout) :: mf(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3))
        real(amrex_real),   intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
        real(amrex_real),   intent(in   ) :: mu(mu_lo(1):mu_hi(1),mu_lo(2):mu_hi(2),mu_lo(3):mu_hi(3))

        real(amrex_real),   intent(in   ) :: prob_lo(3), prob_hi(3)
        real(amrex_real),   intent(in   ) :: dx(3), epsilon_liq, epsilon_sld

        ! local
        integer i, j, k
        ! distribution: exp((mu - epsilon)), here mu = mu(real) / RT, epsilon = e(real) / RT
        ! interpolation: h(x) = x^3 (6 x^2 - 15 x + 10) 
        real exp_liq, exp_sld, exp_mu, interpolation
        exp_liq = exp(epsilon_liq)
        exp_sld = exp(epsilon_sld)

        do k = mf_lo(3), mf_hi(3)
            do j = mf_lo(2), mf_hi(2)
                do i = mf_lo(1), mf_hi(1)
                        ! exp_mu = exp(mu(i, j, k))
                        interpolation = phi(i, j, k) ** 3 * (6 * phi(i, j, k) ** 2 - 15 * phi(i, j, k) + 10)
                        mf(i, j, k) = 1 - phi(i, j, k)
                end do  ! end i
            end do  ! end j
        end do  ! end k


    end subroutine update_solute


end module init_module
    