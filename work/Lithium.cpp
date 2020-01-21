#include <algorithm>
#include <cmath>
#include <iostream>
#include <ctime>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLNodeLaplacian.H>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Utility.H>
#include <AMReX_Cluster.H>
#include <AMReX_BCUtil.H>

#include "Lithium_F.H"
#include "Lithium.H"

using namespace amrex;

Lithium::Lithium()
{       
    ReadParameters();
    ResizeLevelList();   
    InitData();

    WritePlotFile();

    istep[0] += 1;
    // UpdatePotential();

    for (; istep[0] <= max_step; istep[0]++)
    {   
        UpdatePhi();
        if (istep[0] % update_pot_interval == 0 ) {
            UpdatePotential();
        }
        UpdateSolute();

        if ((istep[0] % plot_step == 0 && istep[0] > start_write_plotfile) || (istep[0] == 1))
        {
            WritePlotFile();
        }

        if (istep[0] == switch_step && switch_enable > 0) 
        {
            ChangeVoltage();
        }
    }
}


void Lithium::UpdatePotential()
{   
    const int nlevels = geom.size();
       
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    MLABecLaplacian mlabec(geom, grids, dmap, info);

    mlabec.setMaxOrder(linop_maxorder);
    
    mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                LinOpBCType::Neumann,
                                LinOpBCType::Neumann)},
                    {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                LinOpBCType::Neumann,
                                LinOpBCType::Neumann)});
    for (int lev = 0; lev <= max_level; ++lev)
    { 
        FillDirechletBoundary(potential[lev], geom[lev], bc, bc_val_potential);
        mlabec.setLevelBC(lev, &potential[lev]);
    }
    
    mlabec.setScalars(ascalar, bscalar);
    
    for (int lev = 0; lev <= max_level; ++lev)
    {
        FillDirechletBoundary(phi[lev], geom[lev], bc, bc_val_phi);
        for (MFIter mfi(phi[lev]); mfi.isValid(); ++mfi)
        {
            const Box &bx = mfi.validbox();
        
            init_bcoef(BL_TO_FORTRAN_BOX(bx),
                    BL_TO_FORTRAN_ANYD(bcoef[lev][mfi]),
                    BL_TO_FORTRAN_ANYD(phi[lev][mfi]),
                    geom[lev].ProbLo(), geom[lev].ProbHi(), geom[lev].CellSize(),
                    &cond_liq, &cond_sld
            );
        }
        
        Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(bcoef[lev].boxArray(),
                                                IntVect::TheDimensionVector(idim));
            face_bcoef[idim].define(ba, bcoef[lev].DistributionMap(), 1, 0);
        }
        bcoef[lev].FillBoundary();
        amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef), bcoef[lev], geom[lev]);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            face_bcoef[idim].FillBoundary(geom[0].periodicity());
        }
        
        mlabec.setBCoeffs(lev, amrex::GetArrOfConstPtrs(face_bcoef));
        phi_dt[lev].FillBoundary();
        MultiFab::Copy(rhs[lev], phi_dt[lev], 0, 0, 1, 0);
        rhs[lev].mult(fara_cs);
    }

    MLMG mlmg(mlabec);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setFixedIter(fix_inter);
    mlmg.setBottomTolerance(tol_bottom);
    // mlmg.setBottomSolver(MLMG::BottomSolver::smoother);
#ifdef AMREX_USE_HYPRE
            if (use_hypre) {
                mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
                mlmg.setHypreInterface(hypre_interface);
            }
#endif

    if (first_step_flag){
        mlmg.setFixedIter(30);
        mlmg.solve(GetVecOfPtrs(potential), GetVecOfConstPtrs(rhs), 1e-20, 0);  
        first_step_flag = false;  

    } else {
        mlmg.solve(GetVecOfPtrs(potential), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);    
    }

    // mlmg.getGradSolution(amrex::GetVecOfArrOfPtrs(grad));

    mlmg.compResidual(GetVecOfPtrs(error), GetVecOfPtrs(potential), GetVecOfConstPtrs(rhs));
    error_norm = error[0].norm0();
    AMREX_ALWAYS_ASSERT(error_norm < 1.0);
}

void Lithium::UpdatePhi()
{   

    for (int lev = 0; lev <= max_level; lev++){
        const Box& domain_box = geom[lev].Domain();
        FillDirechletBoundary(phi[lev], geom[lev], bc, bc_val_phi);
        FillDirechletBoundary(potential[lev], geom[lev], bc, bc_val_potential);
        FillDirechletBoundary(solute[lev] , geom[lev], bc, bc_val_solute);

        phi_dt[lev].FillBoundary(geom[lev].periodicity());
        
        for (MFIter mfi(phi[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            advance_phase_field(
                BL_TO_FORTRAN_BOX(bx),
                BL_TO_FORTRAN_BOX(domain_box),
                BL_TO_FORTRAN_ANYD(phi[lev][mfi]),
                BL_TO_FORTRAN_ANYD(phi_dt[lev][mfi]),
                BL_TO_FORTRAN_ANYD(solute[lev][mfi]),
                BL_TO_FORTRAN_ANYD(potential[lev][mfi]),
                BL_TO_FORTRAN_ANYD(phi_new[lev][mfi]),
                BL_TO_FORTRAN_ANYD(output[lev][mfi]),
                geom[lev].CellSize(),
                & dt,
                bc[0].data(),
                & itf_thickness,
                & itf_mobi,
                & nFRT,
                & alpha_asy,
                & voltage
            );
        }
        std::swap(phi[lev], phi_new[lev]);
    }
}


void Lithium::UpdateSolute()
{
    
    for (int lev = 0; lev <= max_level; lev++){
        
        FillDirechletBoundary(phi[lev], geom[lev], bc, bc_val_phi);
        FillDirechletBoundary(potential[lev] , geom[lev], bc, bc_val_potential);
        FillDirechletBoundary(solute[lev] , geom[lev], bc, bc_val_solute);

        const Box& domain_box = geom[lev].Domain();

        for (MFIter mfi(solute[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            advance_solute(
                BL_TO_FORTRAN_BOX(bx),
                BL_TO_FORTRAN_BOX(domain_box),
                BL_TO_FORTRAN_ANYD(phi[lev][mfi]),
                BL_TO_FORTRAN_ANYD(phi_dt[lev][mfi]),
                BL_TO_FORTRAN_ANYD(solute[lev][mfi]),
                BL_TO_FORTRAN_ANYD(potential[lev][mfi]),
                BL_TO_FORTRAN_ANYD(solute_new[lev][mfi]),
                BL_TO_FORTRAN_ANYD(output[lev][mfi]),
                geom[lev].CellSize(), 
                & dt,
                bc[0].data(),
                & c_0,
                & diff_sld,
                & diff_liq,
                & nFRT 
            );           
        }
        std::swap(solute[lev], solute_new[lev]);
        AMREX_ALWAYS_ASSERT(solute[lev].min(0) >= 0.0);

    }
}

void Lithium::FillPhyBndDir(Vector<MultiFab>& mf, int dir, Real bc_val)
{
    for (int lev = 0; lev <= max_level; lev++){
        const Box& domain_box = geom[lev].Domain();

        for (MFIter mfi(mf[lev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();

            fill_physical_boundary_dir(BL_TO_FORTRAN_BOX(bx),
                BL_TO_FORTRAN_BOX(domain_box),
                BL_TO_FORTRAN_ANYD(mf[lev][mfi]),
                & dir,
                & bc_val);
        }
    }
    
}

void Lithium::SetLoValue(int dir, const amrex::Real bc_val_single, Vector<amrex::Real>& this_val){
    AMREX_ALWAYS_ASSERT(this_val.size() == AMREX_SPACEDIM * 2);
    AMREX_ALWAYS_ASSERT(dir > 0 and dir <= AMREX_SPACEDIM);
    this_val[dir-1] = bc_val_single;

}

void Lithium::SetHiValue(int dir, const amrex::Real bc_val_single, Vector<amrex::Real>& this_val){
    AMREX_ALWAYS_ASSERT(this_val.size() == AMREX_SPACEDIM * 2);
    AMREX_ALWAYS_ASSERT(dir > 0 and dir <= AMREX_SPACEDIM);
    this_val[dir + AMREX_SPACEDIM-1] = bc_val_single;
}

void Lithium::FillDirechletBoundary (MultiFab& mf, const Geometry& geom, const Vector<BCRec>& bc, 
                                    const Vector<amrex::Real>& bc_val, const int nComp)
{
    if (Geometry::isAllPeriodic()) return;
    if (mf.nGrow() == 0) return;
    AMREX_ALWAYS_ASSERT(mf.ixType().cellCentered());
    AMREX_ALWAYS_ASSERT(mf.nComp() >= nComp + 1);
    mf.FillBoundary(geom.periodicity());

    const Box& domain_box = geom.Domain();
    Box grown_domain_box = domain_box;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            grown_domain_box.grow(idim,mf.nGrow());
        }
    }
    // Inside grown_domain_box, we have good data.

    const Real* dx = geom.CellSize();
    const Real* prob_lo = Geometry::ProbLo();
    
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = mf[mfi];
        const Box& fab_box = fab.box(); // including ghost cells

        if (! grown_domain_box.contains(fab_box))
        {
            amrex_user_fab_filcc(BL_TO_FORTRAN_N_ANYD(fab,nComp),
                BL_TO_FORTRAN_BOX(domain_box),
                dx, 
                prob_lo,
                bc[0].data(),
                bc_val.dataPtr());
        }
    }
}


void Lithium::ResizeLevelList()
{
    int nlev_max = max_level + 1;

    phi         .resize(nlev_max);
    phi_new     .resize(nlev_max);
    phi_dt      .resize(nlev_max);
    mu          .resize(nlev_max);
    solute_new  .resize(nlev_max);
    solute      .resize(nlev_max);
    potential   .resize(nlev_max);

    rhs         .resize(nlev_max);
    acoef       .resize(nlev_max);
    bcoef       .resize(nlev_max);
    error       .resize(nlev_max);
    grad        .resize(nlev_max);
    flux        .resize(nlev_max);
    output      .resize(nlev_max);

    t_new       .resize(nlev_max, 0.0);
    t_old       .resize(nlev_max, 0.0);

    istep       .resize(nlev_max, 0);
    grad        .resize(1);

    bc          .resize(AMREX_SPACEDIM * 2);
    bc_lo       .resize(AMREX_SPACEDIM);
    bc_hi       .resize(AMREX_SPACEDIM);

}


void Lithium::ReadParameters()
{
    {
        ParmParse pp("li"); // init lithium parameters
        pp.get("itf_position", itf_position);
        pp.get("voltage", voltage);
        pp.get("voltage_positive", voltage_positive);
        pp.get("itf_mobi", itf_mobi);
        pp.get("itf_thickness", itf_thickness);
        pp.get("diff_sld", diff_sld);
        pp.get("diff_liq", diff_liq);
        pp.get("cond_sld", cond_sld);
        pp.get("cond_liq", cond_liq);
        pp.get("temperature", temperature);

        pp.get("tol_abs", tol_abs);
        pp.get("tol_rel", tol_rel);
        pp.get("verbose", verbose);
        pp.get("max_iter", max_iter);
        pp.get("fix_inter", fix_inter);
        pp.get("max_coarsening_level", max_coarsening_level);
        pp.get("tol_bottom", tol_bottom);
        pp.get("bottom_verbose", bottom_verbose);
        pp.get("presmooth", presmooth);
        pp.get("postsmooth", postsmooth);
        pp.get("update_pot_interval", update_pot_interval);

        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);

    }
    {
        ParmParse pp; // parameters controlling the calculating
        pp.get("max_step", max_step);
        pp.get("plot_step", plot_step);
        pp.get("dt", dt);
        pp.get("dt_init", dt_init);
        pp.get("chemical_ratio", chemical_ratio);
        pp.get("start_write_plotfile", start_write_plotfile);
        pp.get("plot_pre_step", plot_pre_step);
        pp.get("pre_step", pre_step);
        pp.get("switch_step", switch_step);
        pp.get("switch_enable", switch_enable);
    }

    // in case changing temperature
    nFRT = ntrans * faraday / (gas * temperature);
    reciprocal_c_0 = 1.0 + exp(epsilon_liq);
    // barrier_height = 12.0 * surf_tension * itf_thickness;
    // grad_energy_coef = 3.0 / 2.0 * surf_tension / itf_thickness;

    for (int n = 0; n < 1; ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            // lo-side BCs
            if (bc_lo[idim] == INT_DIR) {
                bc[n].setLo(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (bc_lo[idim] == FOEXTRAP) {
                bc[n].setLo(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (bc_lo[idim] == EXT_DIR) {
                bc[n].setLo(idim, BCType::ext_dir);  // external Dirichlet
            }
            else {
                amrex::Abort("Invalid bc_lo");
            }

            // hi-side BCs
            if (bc_hi[idim] == INT_DIR) {
                bc[n].setHi(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (bc_hi[idim] == FOEXTRAP) {
                bc[n].setHi(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (bc_hi[idim] == EXT_DIR) {
                bc[n].setHi(idim, BCType::ext_dir);  // external Dirichlet
            }
            else {
                amrex::Abort("Invalid bc_hi");
            }

        }
    }
    bc_val_potential.resize(AMREX_SPACEDIM * 2, 0);
    SetHiValue(1, voltage_positive, bc_val_potential);
    SetLoValue(1, voltage, bc_val_potential);

    bc_val_phi.resize(AMREX_SPACEDIM * 2, 0);
    SetHiValue(1,       0, bc_val_phi);
    SetLoValue(1,       1, bc_val_phi);

    bc_val_solute.resize(AMREX_SPACEDIM * 2, 0);
    SetHiValue(1,       1, bc_val_solute);
    SetLoValue(1,       0, bc_val_solute);    
    
    bc_val_mu.resize(AMREX_SPACEDIM * 2, 0);
    SetHiValue(1,       0, bc_val_mu);
    SetLoValue(1,       0, bc_val_mu);

}

void Lithium::InitData()
{

    // temporary init for single level only
    grids[0] = MakeBaseGrids();
    DistributionMapping dm(grids[0]);
    SetDistributionMap(0, dm);

    for (int lev = 0; lev <= max_level; lev++)
    {
        phi        [lev].define(grids[lev], dmap[lev], 1, 1);
        phi_new    [lev].define(grids[lev], dmap[lev], 1, 1);
        phi_dt     [lev].define(grids[lev], dmap[lev], 1, 1);
        mu         [lev].define(grids[lev], dmap[lev], 1, 1);
        solute_new [lev].define(grids[lev], dmap[lev], 1, 1);
        solute     [lev].define(grids[lev], dmap[lev], 1, 1);
        potential  [lev].define(grids[lev], dmap[lev], 1, 1);
        rhs        [lev].define(grids[lev], dmap[lev], 1, 0);
        acoef      [lev].define(grids[lev], dmap[lev], 1, 0);
        bcoef      [lev].define(grids[lev], dmap[lev], 1, 1);
        error      [lev].define(grids[lev], dmap[lev], 1, 0);
        output     [lev].define(grids[lev], dmap[lev], 1, 0);

        for(int idim = 0; idim< AMREX_SPACEDIM; ++idim)
        {
            grad[lev][idim].define(amrex::convert(grids[0],IntVect::TheDimensionVector(idim)),
                                        dmap[0], 1, 0, MFInfo());
        }

    }

    //
    for (int lev = 0; lev <= max_level; lev++)
    {
        InitDataOnLevel(lev);
    }
    
    // pass the boundary condition to bc
    for (int n = 0; n < bc.size(); ++n)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            
            // lo-side BCs
            if (bc_lo[idim] == INT_DIR) {
                bc[n].setLo(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (bc_lo[idim] == FOEXTRAP) {
                bc[n].setLo(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (bc_lo[idim] == EXT_DIR) {
                bc[n].setLo(idim, BCType::ext_dir);  // external Dirichlet
            }
            else {
                amrex::Abort("Invalid bc_lo");
            }

            // hi-side BCs
            if (bc_hi[idim] == INT_DIR) {
                bc[n].setHi(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (bc_hi[idim] == FOEXTRAP) {
                bc[n].setHi(idim, BCType::foextrap); // first-order extrapolation
            }
            else if (bc_hi[idim] == EXT_DIR) {
                bc[n].setHi(idim, BCType::ext_dir);  // external Dirichlet
            }
            else {
                amrex::Abort("Invalid bc_hi");
            }

        }
    }

}

void Lithium::InitDataOnLevel(int lev)
{


    bcoef       [lev].setVal(0);
    rhs         [lev].setVal(0);
    phi         [lev].setVal(0);
    solute      [lev].setVal(0);
    output      [lev].setVal(0);
    acoef       [lev].setVal(0);
    phi_dt      [lev].setVal(0);
    mu          [lev].setVal(0);
    potential   [lev].setVal(0);
    error       [lev].setVal(0);
    phi_new     [lev].setVal(0);
    solute_new  [lev].setVal(0);

    for (MFIter mfi(phi[lev]); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.validbox();
        const Box &gbx = mfi.fabbox();

        init_phi(BL_TO_FORTRAN_BOX(bx),
            BL_TO_FORTRAN_ANYD(phi[lev][mfi]),
            geom[lev].ProbLo(), geom[lev].ProbHi(), geom[lev].CellSize(),
            &itf_position);

        update_solute(BL_TO_FORTRAN_BOX(bx),
            BL_TO_FORTRAN_ANYD(solute[lev][mfi]),
            BL_TO_FORTRAN_ANYD(phi[lev][mfi]),
            BL_TO_FORTRAN_ANYD(mu[lev][mfi]),
            geom[lev].ProbLo(), 
            geom[lev].ProbHi(), 
            geom[lev].CellSize(),
            & epsilon_liq,
            & epsilon_sld
        );
    }
}


Vector<std::string>
Lithium::PlotFileVarNames() const
{
    Vector<std::string> vname; // corresponding variable names

    vname.push_back("phi");
    vname.push_back("mu");
    vname.push_back("solute");
    vname.push_back("potential");
    vname.push_back("rhs");
    vname.push_back("phi_dt");
    vname.push_back("bcoef");
    vname.push_back("phi_new");
    vname.push_back("error");
    vname.push_back("output");

    return vname;
}


void
Lithium::ChangeVoltage()
{   
    voltage = 0.1;
    for (int lev = 0; lev <= max_level; lev++)
    {
        potential[lev].setVal(0);
        phi_dt[lev].setVal(0);
    }
    SetLoValue(1, voltage, bc_val_potential);
    first_step_flag = true;
    // verbose = 1;
    UpdatePotential();

    plot_step = 1000;
}

std::string
Lithium::PlotFileName(int lev) const
{
    return amrex::Concatenate(plot_file, lev, 7);
}

void Lithium::WritePlotFile()
{
    const std::string &plotfilename = PlotFileName(istep[0]);
    // const auto& mf = PlotFileMF();
    const auto &varnames = PlotFileVarNames();

    Vector<MultiFab> plotData(max_level + 1);

    // copy the components to plotData to write to plotfile
    for (int lev = 0; lev <= max_level; lev++)
    {

        plotData[lev].define(grids[lev], dmap[lev], varnames.size(), 0);
        MultiFab::Copy(plotData[lev], phi[lev],       0, 0, 1, 0);
        MultiFab::Copy(plotData[lev], mu[lev],        0, 1, 1, 0);
        MultiFab::Copy(plotData[lev], solute[lev],    0, 2, 1, 0);
        MultiFab::Copy(plotData[lev], potential[lev], 0, 3, 1, 0);
        MultiFab::Copy(plotData[lev], rhs[lev],       0, 4, 1, 0);
        MultiFab::Copy(plotData[lev], phi_dt[lev],    0, 5, 1, 0);
        MultiFab::Copy(plotData[lev], bcoef[lev],     0, 6, 1, 0);
        MultiFab::Copy(plotData[lev], phi_new[lev],   0, 7, 1, 0);
        MultiFab::Copy(plotData[lev], error[lev],     0, 8, 1, 0);
        MultiFab::Copy(plotData[lev], output[lev],    0, 9, 1, 0);
    }
    // istep[0] = istep[0] + 1;
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    // amrex::WriteMultiLevelPlotfile(plotfilename, max_level+1, mf, varnames,
    //         Geom(), t_new[0], istep, refRatio());
    amrex::WriteMultiLevelPlotfile(plotfilename,
                                   max_level + 1,
                                   GetVecOfConstPtrs(plotData),
                                   varnames,
                                   Geom(),
                                   t_new[0],
                                   istep,
                                   refRatio());
}


// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
// overrides the pure virtual function in AmrCore
void Lithium::MakeNewLevelFromCoarse(int lev, Real time, const BoxArray &ba,
                                    const DistributionMapping &dm)
{
}

// Remake an existing level using provided BoxArray and DistributionMapping and
// fill with existing fine and coarse data.
// overrides the pure virtual function in AmrCore
void Lithium::RemakeLevel(int lev, Real time, const BoxArray &ba,
                        const DistributionMapping &dm)
{
}

// Delete level data
// overrides the pure virtual function in AmrCore
void Lithium::ClearLevel(int lev)
{
}

// Make a new level from scratch using provided BoxArray and DistributionMapping.
// Only used during initialization.
// overrides the pure virtual function in AmrCore
void Lithium::MakeNewLevelFromScratch(int lev, Real time, const BoxArray &ba,
                                      const DistributionMapping &dm)
{
}

// tag all cells for refinement
// overrides the pure virtual function in AmrCore
void Lithium::ErrorEst(int lev, TagBoxArray &tags, Real time, int ngrow)
{
}
