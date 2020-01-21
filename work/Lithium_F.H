#ifndef LITHIUM_F_H_
#define LITHIUM_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif
    // phi: tanh((itf_position-x)*200)
    void init_phi(const int* lo, const int* hi,
                    BL_FORT_FAB_ARG_3D(mf),
                    const amrex_real* prob_lo, 
                    const amrex_real* prob_hi, 
                    const amrex_real* dx,
                    const amrex_real* itf_position
    );

    // bcoef: interpolation + (1 - interpolation) * cond_liq / cond_sld
    void init_bcoef(const int* lo, const int* hi,
                    BL_FORT_FAB_ARG_3D(mf),
                    BL_FORT_FAB_ARG_3D(phi),
                    const amrex_real* prob_lo, 
                    const amrex_real* prob_hi, 
                    const amrex_real* dx,
                    const amrex_real* cond_liq,
                    const amrex_real* cond_sld
    );
                                       
    void init_potential(const int* lo, const int* hi,
                    BL_FORT_FAB_ARG_3D(mf),
                    BL_FORT_FAB_ARG_3D(phi),
                    const amrex_real* prob_lo, 
                    const amrex_real* prob_hi, 
                    const amrex_real* dx,
                    const amrex_real* itf_position
    );
    
    void update_solute(const int* lo, const int* hi,
                    BL_FORT_FAB_ARG_3D(solute),
                    BL_FORT_FAB_ARG_3D(phi),
                    BL_FORT_FAB_ARG_3D(mu),
                    const amrex_real* prob_lo, 
                    const amrex_real* prob_hi, 
                    const amrex_real* dx,
                    const amrex_real* epsilon_liq,
                    const amrex_real* epsilon_sld

    );

    void fill_physical_boundary_dir(
                    const int* lo, const int* hi,
                    const int* domlo, const int* domhi,
                    BL_FORT_FAB_ARG_3D(mf),
                    const int* dir,
                    const amrex_real* bc_val
    );



    void advance_phase_field(const int* lo, const int* hi,
                    const int* domlo, const int* domhi,
                    BL_FORT_FAB_ARG_3D(phi),
                    BL_FORT_FAB_ARG_3D(phi_dt),
                    BL_FORT_FAB_ARG_3D(solute),
                    BL_FORT_FAB_ARG_3D(potential),
                    BL_FORT_FAB_ARG_3D(result),
                    BL_FORT_FAB_ARG_3D(output),
                    const amrex_real* dx,
                    const amrex_real* dt,
                    const int* bclo,
                    const amrex_real* itf_thickness,
                    const amrex_real* itf_mobi,
                    const amrex_real* nFRT,
                    const amrex_real* alpha,
                    const amrex_real* voltage
    );

    void advance_solute(const int* lo, const int* hi,
                    const int* domlo, const int* domhi,
                    BL_FORT_FAB_ARG_3D(phi),
                    BL_FORT_FAB_ARG_3D(phi_dt),
                    BL_FORT_FAB_ARG_3D(solute),
                    BL_FORT_FAB_ARG_3D(potential),
                    BL_FORT_FAB_ARG_3D(result),
                    BL_FORT_FAB_ARG_3D(output),
                    const amrex_real* dx,
                    const amrex_real* dt,
                    const int* bclo,
                    const amrex_real* c_0,
                    const amrex_real* diff_sld,
                    const amrex_real* diff_liq,
                    const amrex_real* nFRT
    );


    void amrex_user_fab_filcc (amrex_real* q, 
                    const int* qlo, 
                    const int* qhi, 
                    const int* domlo, 
                    const int* domhi,
                    const amrex_real* dx, 
                    const amrex_real* xlo,
                    const int* bclo,
                    const amrex_real* bc_val
                    );
    
    void average_smoother(const int* lo, const int* hi,
                        BL_FORT_FAB_ARG_3D(phi)

                    );

#ifdef __cplusplus
}
#endif

#endif