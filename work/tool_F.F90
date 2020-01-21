module tool_mod

    use amrex_fort_module, only : amrex_spacedim, amrex_real
    use amrex_bc_types_module

    implicit none

    public
    
contains
    ! grad mf[size=1, cell-center] -> flux[size=3, face-center] 
    subroutine compute_flux (lo, hi, dom_lo, dom_hi, &
        phi, philo, phihi, &
        flux, fxlo, fxhi, &
        dx, bc) bind(C, name="compute_flux")
        integer lo(3), hi(3), dom_lo(3), dom_hi(3)
        integer philo(3), phihi(3), fxlo(3), fxhi(3)
        real(amrex_real), intent(in)    :: phi  (philo(1):phihi(1), philo(2):phihi(2), philo(3):phihi(3))
        real(amrex_real), intent(inout) :: flux (fxlo(1): fxhi(1), fxlo(2): fxhi(2), fxlo(3): fxhi(3), amrex_spacedim)
        real(amrex_real), intent(in)    :: dx(3)
        integer bc(amrex_spacedim, 2, 1)  ! (dim,lohi,ncomp)

        ! local variables
        integer i,j,k

        ! x-fluxes
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)+1
                    flux(i,j,k,1) = ( phi(i, j, k) - phi(i-1, j, k) ) / dx(1)
                end do
            end do
        end do

        ! lo-x boundary, ghost cell contains value on boundary
        if (dom_lo(1) .eq. lo(1) .and. (bc(1,1,1) .eq. amrex_bc_foextrap .or. bc(1,1,1) .eq. amrex_bc_ext_dir) ) then
            i = lo(1)
            do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                    flux(i,j,k,1) = ( phi(i, j, k) - phi(i-1, j, k) ) / (0.5d0*dx(1))
                end do
            end do
        end if

        ! hi-x boundary, ghost cell contains value on boundary

        if (dom_hi(1) .eq. hi(1) .and. (bc(1,2,1) .eq. amrex_bc_foextrap .or. bc(1,2,1) .eq. amrex_bc_ext_dir) ) then

            i = hi(1) + 1
            do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                    flux(i,j,k,1) = ( phi(i, j, k) - phi(i-1, j, k) ) / (0.5d0*dx(1))
                end do
            end do
        end if

#if (AMREX_SPACEDIM >= 2)
        ! y-fluxes
        do k = lo(3), hi(3)
            do j = lo(2), hi(2)+1
                do i = lo(1), hi(1)
                    flux(i,j,k,2) = ( phi(i,j,k) - phi(i,j-1,k) ) / dx(2)
                end do
            end do
        end do

        ! lo-y boundary, ghost cell contains value on boundary
        if (dom_lo(2) .eq. lo(2) .and. (bc(2,1,1) .eq. amrex_bc_foextrap .or. bc(2,1,1) .eq. amrex_bc_ext_dir) ) then
            j = lo(2)
            do k = lo(3), hi(3)
                do i = lo(1), hi(1)
                    flux(i,j,k,2) = ( phi(i, j, k) - phi(i, j-1, k) ) / (0.5d0*dx(2))
                end do
            end do
        end if

        ! hi-y boundary, ghost cell contains value on boundary
        if (dom_hi(2) .eq. hi(2) .and. (bc(2,2,1) .eq. amrex_bc_foextrap .or. bc(2,2,1) .eq. amrex_bc_ext_dir) ) then
            j = hi(2) + 1
            do k = lo(3), hi(3)
                do i = lo(1), hi(1)
                    flux(i,j,k,2) = ( phi(i, j, k) - phi(i, j-1, k) ) / (0.5d0*dx(2))
                end do
            end do
        end if
#endif

#if (AMREX_SPACEDIM == 3)
        ! z-fluxes
        do k = lo(3), hi(3)+1
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    flux(i,j,k,3) = ( phi(i,j,k) - phi(i,j,k-1) ) / dx(3)
                end do
            end do
        end do

        ! lo-z boundary, ghost cell contains value on boundary
        if (dom_lo(3) .eq. lo(3) .and. (bc(3,1,1) .eq. amrex_bc_foextrap .or. bc(3,1,1) .eq. amrex_bc_ext_dir) ) then
            k = lo(3)
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    flux(i,j,k,3) = ( phi(i,j,k) - phi(i,j,k-1) ) / (0.5d0*dx(3))
                end do
            end do
        end if

        ! hi-z boundary, ghost cell contains value on boundary
        if (dom_hi(3) .eq. hi(3) .and. (bc(3,2,1) .eq. amrex_bc_foextrap .or. bc(3,2,1) .eq. amrex_bc_ext_dir) ) then
            k = hi(3) + 1
                do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    flux(i,j,k,3) = ( phi(i,j,k) - phi(i,j,k-1) ) / (0.5d0*dx(3))
                end do
            end do
        end if

#endif

    end subroutine compute_flux


    ! mf[Domain[dir]] = val
    subroutine fill_physical_boundary_dir(lo, hi, dom_lo, dom_hi, mf, mf_lo, mf_hi, dir, val) &
        bind(C, name="fill_physical_boundary_dir")
        integer lo(3), hi(3), dom_lo(3), dom_hi(3)
        integer mf_lo(3), mf_hi(3)
        integer dir
        real(amrex_real) val
        real(amrex_real), intent(inout) :: mf(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3))

        ! dir: 
        ! sign: "+": hi, "-": lo
        ! num: 1: x, 2:y, 3: z 
        
        ! local var
        integer i, j, k

        select case(dir)
        case(1)
            if(dom_hi(1) .eq. hi(1)) then
                i = hi(1)
                do k = lo(3), hi(3)
                    do j = lo(2), hi(2)
                        mf(i+1, j, k) = val
                    end do
                end do
            end if

        case(-1)
            if(dom_lo(1) .eq. lo(1)) then
                i = lo(1)
                do k = lo(3), hi(3)
                    do j = lo(2), hi(2)
                        mf(i-1, j, k) = val                 
                    end do
                end do
            end if

        case(2)
            if(dom_hi(2) .eq. hi(2)) then
                j = hi(2)
                do k = lo(3), hi(3)
                    do i = lo(i), hi(1)
                        mf(i, j+1, k) = val
                    end do
                end do
            end if

        case(-2)
            if(dom_lo(2) .eq. lo(2)) then
                j = lo(2)
                do k = lo(3), hi(3)
                    do i = lo(i), hi(1)
                        mf(i, j-1, k) = val
                    end do
                end do
            end if

        case(3)
            if(dom_hi(3) .eq. hi(3)) then
                k = hi(3)
                do j = lo(2), hi(2)
                    do i = lo(1), hi(1)
                        mf(i, j, k+1) = val
                    end do
                end do
            end if

        case(-3)
            if(dom_lo(3) .eq. lo(3)) then
                k = lo(3)
                do j = lo(2), hi(2)
                    do i = lo(1), hi(1)
                        mf(i, j, k-1) = val
                    end do
                end do
            end if
        end select
    end subroutine fill_physical_boundary_dir

    subroutine average_smoother(lo, hi, mf, mf_lo, mf_hi)&
        bind(C, name="average_smoother")
        integer lo(3), hi(3)
        integer mf_lo(3), mf_hi(3)
        real(amrex_real), intent(inout) :: mf(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3))

        real(amrex_real)  result(mf_lo(1):mf_hi(1),mf_lo(2):mf_hi(2),mf_lo(3):mf_hi(3))
        integer i, j, k
        
        result = mf
        do k=hi(3),hi(3)
            do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                    result(i, j, k) = 0.5 * mf(i, j, k) + 0.125 * (mf(i+1, j, k) + mf(i-1, j, k) + mf(i, j+1, k) + mf(i, j-1, k))
                end do
            end do
        end do
        mf = result

    end subroutine

    function ran()
        implicit none
        integer, save :: flag = 0
        real :: ran
        if(flag .eq. 0) then
            call random_seed()
            flag = 1
        end if
        call random_number(ran)
    end function ran

end module tool_mod