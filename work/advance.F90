module advance_mod

    ! since this is a .F90 file (instead of .f90) we run this through a C++ preprocessor
    ! for e.g., #if (amrex_spacedim == 1) statements.

    use amrex_fort_module, only : amrex_real, amrex_spacedim
    
    implicit none
    
    public advance_phase_field
    
    contains
    ! cal -L_sigma(g:x - kappa laplacian phi) -(BV)
    ! BV =  L_eta h:x (exp((1-alpha) nF/RT potential) - c_Li/c_0 exp(alpha nF/RT potential))
    subroutine advance_phase_field(lo, hi, &
        dom_lo, dom_hi, &
        phi, phi_lo, phi_hi, &
        phi_dt, phi_dt_lo, phi_dt_hi, &
        solute, solute_lo, solute_hi, & 
        potential, potential_lo, potential_hi, &
        result,  result_lo, result_hi, &
        output,  output_lo, output_hi, &
        dx, dt, bc, &
        itf_mobi,&
        itf_thickness, nFRT, alpha, voltage)&
        bind(C, name="advance_phase_field")

        use tool_mod, only: compute_flux

        integer lo(3), hi(3), dom_hi(3), dom_lo(3)
        integer phi_hi(3), phi_lo(3)
        integer phi_dt_hi(3), phi_dt_lo(3)
        integer solute_hi(3), solute_lo(3)
        integer potential_hi(3), potential_lo(3)
        integer result_hi(3), result_lo(3)
        integer output_hi(3), output_lo(3)


        real(amrex_real), intent(in   )  :: phi(phi_lo(1): phi_hi(1), phi_lo(2): phi_hi(2), phi_lo(3): phi_hi(3))
        real(amrex_real), intent(inout)  :: phi_dt(phi_dt_lo(1): phi_dt_hi(1), phi_dt_lo(2): phi_dt_hi(2), phi_dt_lo(3): phi_dt_hi(3))
        real(amrex_real), intent(in   )  :: potential(potential_lo(1): potential_hi(1), potential_lo(2): potential_hi(2), potential_lo(3): potential_hi(3))
        real(amrex_real), intent(in   )  :: solute(solute_lo(1): solute_hi(1), solute_lo(2): solute_hi(2), solute_lo(3): solute_hi(3))
        real(amrex_real), intent(inout)  :: result(result_lo(1): result_hi(1), result_lo(2): result_hi(2), result_lo(3): result_hi(3))
        real(amrex_real), intent(inout)  :: output(output_lo(1): output_hi(1), output_lo(2): output_hi(2), output_lo(3): output_hi(3))
        real(amrex_real), intent(in   )  :: dx(3), dt
        real(amrex_real), intent(in   )  :: itf_mobi, itf_thickness, nFRT, alpha, voltage
        integer, intent(in)              :: bc(amrex_spacedim,2,1) ! (dim,lohi,ncomp)

        real(amrex_real) flux_phi(lo(1): hi(1)+1, lo(2): hi(2)+1, lo(3): hi(3)+1, amrex_spacedim)
        integer i, j, k, dim
        
        real(amrex_real) d_dwell            ! derivate double well function g:phi(phi) = [ W phi^2 (1 -phi)^2 ]'
        real(amrex_real) d_interpolation    ! derivate interpolation function h:phi(phi) = [ phi^3 (6 phi^2 - 15 phi + 10) ]'
        real(amrex_real) phi_laplacian      ! phi_laplacian = div grad phi
        real(amrex_real) diffusion          ! diffusion = -itf_mobi ( d_dwell - kappa phi_laplacian )

        real(amrex_real) exp1               ! exp((1 - alpha) nF/RT potential_hi)
        real(amrex_real) exp2               ! exp((alpha nF/RT potential_hi)
        real(amrex_real) bv, delta_phi                 ! Bultervolmer = rec_const d_interpolation (exp1 - c_Li/ c_0 exp2)
        real(amrex_real) middle
        real(amrex_real) noise(phi_lo(1): phi_hi(1), phi_lo(2): phi_hi(2), phi_lo(3): phi_hi(3))
        
        
        ! compute flux locally 
        call random_seed()
        call random_number(noise)
        ! print *, noise(lo(1), lo(2): lo(2) + 1, lo(3))
        call compute_flux(lo, hi, dom_hi, dom_lo, phi, phi_lo, phi_hi, flux_phi, lo, hi+1, dx, bc)
        
        do k=lo(3),hi(3)
            do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                    d_dwell         = 24 * phi(i, j, k) * (1 - phi(i, j, k)) * (1 - 2 * phi(i, j, k)) 
                    d_interpolation = 30 * phi(i, j, k) ** 2 * (1 - phi(i, j, k)) ** 2

                    phi_laplacian = flux_phi(i+1, j, k , 1) - flux_phi(i, j, k , 1)
#if (AMREX_SPACEDIM >= 2)
                    phi_laplacian = phi_laplacian + flux_phi(i, j+1, k, 2) - flux_phi(i, j, k , 2)
#endif

#if (AMREX_SPACEDIM == 3)   
                    phi_laplacian = phi_laplacian + flux_phi(i, j, k+1, 3) - flux_phi(i, j, k , 3)
#endif
                    phi_laplacian = 1.5 * phi_laplacian / dx(1)

                    diffusion = - itf_mobi * (d_dwell - itf_thickness * itf_thickness * phi_laplacian)

                    exp1 = exp((1 - alpha) * ( nFRT  * potential(i, j, k) ))
                    exp2 = exp(   - alpha  * ( nFRT  * potential(i, j, k) ))
                    bv   = - d_interpolation * (exp1 -  solute(i, j, k) * exp2)


                    ! phi_dt(i, j, k) = (diffusion + bv) ! * (1 + noise(i, j, k) * 0.005d0)
                    phi_dt(i, j, k) = (diffusion + bv) ! * (1 + noise(i, j, k) * 0.005d0)
                    ! phi_dt(i, j, k) = bv
                    ! bv contains minus-hypen
                    result(i, j, k) = phi(i, j, k) + phi_dt(i, j, k) * dt
                    if (result(i, j, k) .gt. 1.d0) then
                        result(i, j, k) = 1.d0
                    else if (result(i, j, k) .lt. 0.d0) then
                        result(i, j, k) = 0.d0
                    end if
                    ! result(i, j, k) = phi(i, j, k) + phi_laplacian * dt

                end do ! i
            end do ! j
        end do ! k

    end subroutine advance_phase_field

    ! cal \chi mu:t = [del dot D c_Li (grad mu + nFRT grad potenttial)] - [h:t (c^s C^s_m / C^l_m - c^l)] 
    subroutine advance_solute(lo, hi, &
        dom_lo, dom_hi, &
        phi, phi_lo, phi_hi, &
        phi_dt, phi_dt_lo, phi_dt_hi, &
        solute, solute_lo, solute_hi, & 
        potential, potential_lo, potential_hi, &
        result,  result_lo, result_hi, &
        output,  output_lo, output_hi, &
        dx, dt, bc,&
        c_0, diff_sld, diff_liq, nFRT &
        )&
        bind(C, name="advance_solute")

        use tool_mod, only: compute_flux

        integer lo(3), hi(3), dom_hi(3), dom_lo(3)
        integer phi_hi(3), phi_lo(3)
        integer phi_dt_hi(3), phi_dt_lo(3)
        integer solute_hi(3), solute_lo(3)
        integer potential_hi(3), potential_lo(3) 
        integer result_hi(3), result_lo(3)
        integer output_hi(3), output_lo(3)

        real(amrex_real), intent(in   )  :: phi(phi_lo(1): phi_hi(1), phi_lo(2): phi_hi(2), phi_lo(3): phi_hi(3))
        real(amrex_real), intent(in   )  :: phi_dt(phi_dt_lo(1): phi_dt_hi(1), phi_dt_lo(2): phi_dt_hi(2), phi_dt_lo(3): phi_dt_hi(3))
        real(amrex_real), intent(in   )  :: potential(potential_lo(1): potential_hi(1), potential_lo(2): potential_hi(2), potential_lo(3): potential_hi(3))
        real(amrex_real), intent(in   )  :: solute(solute_lo(1): solute_hi(1), solute_lo(2): solute_hi(2), solute_lo(3): solute_hi(3))
        real(amrex_real), intent(inout)  :: result(result_lo(1): result_hi(1), result_lo(2): result_hi(2), result_lo(3): result_hi(3))
        real(amrex_real), intent(inout)  :: output(output_lo(1): output_hi(1), output_lo(2): output_hi(2), output_lo(3): output_hi(3))

        real(amrex_real), intent(in   )  :: dx(3), dt, nFRT, c_0, diff_sld, diff_liq
        integer, intent(in)              :: bc(amrex_spacedim,2,1) ! (dim,lohi,ncomp)
        
        ! local
        real(amrex_real)  flux_potential(lo(1): hi(1)+1, lo(2): hi(2)+1, lo(3): hi(3)+1, amrex_spacedim)
        real(amrex_real)  flux_solute   (lo(1): hi(1)+1, lo(2): hi(2)+1, lo(3): hi(3)+1, amrex_spacedim)
        real(amrex_real)  diff          (lo(1): hi(1)+1, lo(2): hi(2)+1, lo(3): hi(3)+1, amrex_spacedim)
        real(amrex_real)  solute_face   (lo(1): hi(1)+1, lo(2): hi(2)+1, lo(3): hi(3)+1, amrex_spacedim)

        integer i, j, k
        real(amrex_real)  interpolation     ! h(x) = x^3 (6 x^2 - 15x + 10)
        real(amrex_real)  d_interpolation   ! h:x(x) = 30  x^2 (1 - x)^2
        real(amrex_real)  flux_total(3, 2)  ! 1st index: x y z; 2nd index: lo hi
        real(amrex_real)  laplacian         ! div [diff solute (grad mu + nFRT grad phi)], mu = mu(real) / RT
        real(amrex_real)  migration         ! h:phi phi:t * (c^s sdt - c^l), sdt:site_density_ratio
        real(amrex_real)  temp_phi, temp_interpolation

        call compute_flux(lo, hi, dom_lo, dom_hi, solute, solute_lo, solute_hi, flux_solute, lo, hi+1, dx, bc)
        call compute_flux(lo, hi, dom_lo, dom_hi, potential, potential_lo, potential_hi, flux_potential, lo, hi+1, dx, bc)

        ! interpolation in x
        do k=lo(3),hi(3)
            do j=lo(2),hi(2)
                do i=lo(1),hi(1) + 1
                    temp_phi            = (phi(i, j, k) + phi(i-1, j, k)) / 2
                    temp_interpolation  = temp_phi ** 3 * (6 * temp_phi ** 2 - 15 * temp_phi + 10)
                    diff(i, j, k, 1)    = diff_sld * temp_interpolation + diff_liq * (1 - temp_interpolation)
                    solute_face(i, j, k, 1) = (solute(i, j, k) + solute(i-1, j, k)) / 2
                end do 
            end do
        end do

#if (AMREX_SPACEDIM >= 2)
        do k=lo(3),hi(3)
            do j=lo(2),hi(2) + 1
                do i=lo(1),hi(1)
                    temp_phi            = (phi(i, j, k) + phi(i, j-1, k)) / 2
                    temp_interpolation  = temp_phi ** 3 * (6 * temp_phi ** 2 - 15 * temp_phi + 10)
                    diff(i, j, k, 2)    = diff_sld * temp_interpolation + diff_liq * (1 - temp_interpolation)
                    solute_face(i, j, k, 2) = (solute(i, j, k) + solute(i, j-1, k)) / 2
                end do 
            end do
        end do
#endif

#if (AMREX_SPACEDIM == 3)
        ! interpolation in z
        do k=lo(3),hi(3) + 1
            do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                    temp_phi           = (phi(i, j, k) + phi(i, j, k-1)) / 2
                    temp_interpolation = temp_phi ** 3 * (6 * temp_phi ** 2 - 15 * temp_phi + 10)
                    diff(i, j, k, 3)   = diff_sld * temp_interpolation + diff_liq * (1 - temp_interpolation)
                    solute_face(i, j, k, 3) = (solute(i, j, k) + solute(i, j, k-1)) / 2
                end do 
            end do
        end do
#endif

        ! flux(i - 1, j, k, 1) stores wrong numbers (not initilized)!
        flux_total = 0
        do k=lo(3),hi(3)
            do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                    flux_total(1, 1) = diff(i  , j, k, 1) * ( flux_solute(i  , j, k, 1)  +  nFRT * flux_potential(i  , j, k, 1) * solute_face(i  , j, k, 1) )
                    flux_total(1, 2) = diff(i+1, j, k, 1) * ( flux_solute(i+1, j, k, 1)  +  nFRT * flux_potential(i+1, j, k, 1) * solute_face(i+1, j, k, 1) )
                    laplacian = flux_total(1, 2) - flux_total(1, 1)
#if (AMREX_SPACEDIM >= 2)
                    flux_total(2, 1) = diff(i, j  , k, 2) * ( flux_solute(i, j  , k, 2)  +  nFRT * flux_potential(i, j  , k, 2) * solute_face(i, j  , k, 2) )
                    flux_total(2, 2) = diff(i, j+1, k, 2) * ( flux_solute(i, j+1, k, 2)  +  nFRT * flux_potential(i, j+1, k, 2) * solute_face(i, j+1, k, 2) )
                    laplacian = laplacian + flux_total(2, 2) - flux_total(2, 1)
#endif
                
#if (AMREX_SPACEDIM == 3)   
                    flux_total(3, 1) = diff(i, j, k  , 3) * ( flux_solute(i, j, k  , 3)  +  nFRT * flux_potential(i, j, k  , 3) * solute_face(i, j, k  , 3) )
                    flux_total(3, 2) = diff(i, j, k+1, 3) * ( flux_solute(i, j, k+1, 3)  +  nFRT * flux_potential(i, j, k+1, 3) * solute_face(i, j, k+1, 3) )
                    laplacian = laplacian + flux_total(3, 2) - flux_total(3, 1)
#endif
                    laplacian = laplacian / dx(1)
                    !end laplacian
                    output(i, j, k) = (laplacian - migration) * dt
                    migration = phi_dt(i, j, k) *  76.4
                    !end migration
                    result(i, j, k) = solute(i, j, k) + (laplacian - migration) * dt
                    if (result(i, j, k) .lt. 0d0) then
                        result(i, j, k) = 0.d0
                    end if
                    
                end do ! i
            end do ! j
        end do ! k
        
    end subroutine advance_solute

end module advance_mod