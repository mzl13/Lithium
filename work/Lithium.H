#ifndef LITHIUM_H_
#define LITHIUM_H_

#include <map>
#include <istream>
#include <AMReX_AmrCore.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_Print.H>

using namespace amrex;

class Lithium : public AmrCore{

public:

    Lithium();
    virtual ~Lithium() {}

    // Make a new level using provided BoxArray and DistributionMapping and 
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
            const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and 
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
            const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
            const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

private:
    void ReadParameters(); 
    void ResizeLevelList();
    void WritePlotFile();

    // init cell values for multifabs in each level
    void InitData();
    void InitDataOnLevel(int lev);

    // calculate electrical potential distribution in whole domain
    void UpdatePotential();
    void UpdatePhi();
    void UpdateSolute();
    void UpdateChemicalPotential();

    void ChangeVoltage();

    // Fillboundary for related mf
    void UpdateGhostCells();

    // Fill physical boundary in a single direction
    void FillPhyBndDir(amrex::Vector<MultiFab>& mf, int dir, amrex::Real db_val);

    // do FillBoundary(periority) for Multifab and FillDomainBoundary, including direchlet boundary with bc_val for single nCmp
    void FillDirechletBoundary (MultiFab& phi, const Geometry& geom, const amrex::Vector<BCRec>& bc, const amrex::Vector<amrex::Real>& bc_val, const int nComp=0);

    // set the value for direchlet on DomainBoundary value
    void SetHiValue(int dir, const amrex::Real bc_val_single, amrex::Vector<amrex::Real>& this_val);
    void SetLoValue(int dir, const amrex::Real bc_val_single, amrex::Vector<amrex::Real>& this_val);

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // set plotfile variables names
    amrex::Vector<std::string> PlotFileVarNames () const;

    // write plotfile to disk
    void WritePlotFile () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // data holder
    amrex::Vector<MultiFab>  phi, phi_new, phi_dt;  // phase field, new, old and dphi / dt
    amrex::Vector<MultiFab>  solute, solute_new;            // chemical potential, latin mu
    amrex::Vector<MultiFab>  potential;     // Lithium ion density in electrolyte and overpotetial
    amrex::Vector<MultiFab>  mu;

    // parameters for Linearsolver: alpha * a - beta * (del dot b grad)) phi = rhs
    amrex::Vector<MultiFab> rhs;    
    amrex::Vector<MultiFab> acoef;  // is not necessary: set alpha = 0 auto setting a = 0
    amrex::Vector<MultiFab> bcoef;  // b cell-center value, conductivity in this case
    amrex::Vector<MultiFab> error;  // 
    amrex::Vector<MultiFab> output; // 

    amrex::Real ascalar = 0;        // alpha
    amrex::Real bscalar = -1.0;     // beta
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> grad;  // b grad phi
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> flux;  // grad phi
    // amrex::Array<MultiFab,AMREX_SPACEDIM> face_bcoef;  // face-center value for b

    int switch_step;         // switch voltage step interval
    int switch_enable  = 0;  // 0 False; 1 True


    // mode for linear solver
    int  ls_verbose           = 2;
    int  bottom_verbose       = 2;
    int  max_iter             = 100;
    int  fix_inter            = 99;     // end iter even non-convergence result
    int  max_fmg_iter         = 0;
    int  linop_maxorder       = 2;
    bool agglomeration        = true;
    bool consolidation        = true;
    int  max_coarsening_level = 30;
    bool use_hypre            = false;
    bool use_petsc            = false;
    bool first_step_flag      = true;
    int  chemical_ratio       = 100;
    int  presmooth            = 8;
    int  postsmooth           = 8;
    int  update_pot_interval  = 100;    // potential will update for every $ steps
    amrex:: Real error_norm;

    // normally only changing the two parameters 
    amrex::Real tol_rel ;
    amrex::Real tol_abs ;
    amrex::Real tol_bottom;
#ifdef AMREX_USE_HYPRE
    int hypre_interface_i = 3;  // 1. structed, 2. semi-structed, 3. ij
    amrex::Hypre::Interface hypre_interface = amrex::Hypre::Interface::structed;
#endif

    // control the output and steps
    int max_step, pre_step;              
    int plot_step, plot_pre_step;   // write data if istep % plot_step == 0
    amrex::Vector<int> istep;        
    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Real dt;
    amrex::Real dt_init;
    int start_write_plotfile;
    std::string plot_file {"CAL_DATA/plt"};  // relative path and base name of Boxlib

    // geometry parameters
    amrex::Real itf_position       = 0.5;        // initial interface position
    // physical constants
    const  amrex::Real gas         = 8.314;      // J K-1 mol-1
    const  amrex::Real charge      = 1.602e-16;
    const  amrex::Real faraday     = 96485.332;  // s A mol-1
    // physical parameters 
    amrex:: Real voltage            = -0.3;     // voltage on x low end (x_hi) = 0,
    amrex:: Real voltage_positive   = 0.0;     // voltage on x low end (x_hi) = 0,
    amrex:: Real itf_mobi           = 6.25;     // interfacial mobility
    amrex:: Real rec_const          = 1e-3;     // Reaction constant
    amrex:: Real diff_sld           = 317.9;    // Diffusion coeff. in electrolyte
    amrex:: Real diff_liq           = 0.3;      // Diffusion coeff. in electrode
    amrex:: Real cond_sld           = 1e9;      // conductivity in electrode
    amrex:: Real cond_liq           = 1e2;      // conductivity in electrolyte
    amrex:: Real epsilon_liq        = 2.63;     // initial chemical potential difference of lithium &
    amrex:: Real epsilon_sld        = -13.8;    // and neutral components on the electrolyte/electrode
    amrex:: Real energy_density_liq = 6.95e-2;  // free energy density of the electrolyte phases
    amrex:: Real energy_density_sld = 13.8;     // free energy density of the electrode phases
    amrex:: Real temperature        = 300;      // temperatures
    amrex:: Real alpha_asy          = 0.5;      // asymmetry factor
    amrex:: Real surf_tension       = 0.22;
    amrex:: Real itf_thickness      = 5;
    amrex:: Real reciprocal_c_0;                // initial lithium-ion molar ratio
    amrex:: Real c_0 = 1.0;


    int    ntrans                  = 1;          // number of electrons transferred
    // intermedia parameters for simplicity
    amrex::Real nFRT;                            // n * F / R * T, R = 8.314 J mol-1 K-1, F = 
    amrex::Real fara_cs;                         // n * F * c_s / cond_sld
    amrex::Real site_density_ratio = 5.31;       // the ratio of site density in electro -de/-lyte
    amrex::Real grad_energy_coef   = 0.33;       // 3/2 surf_tension * thickness
    amrex::Real barrier_height     = 2.64;       // Interface thickness = 12 * surf_tension / thickness

    // boundary condition

    amrex::Vector<int> bc_lo;
    amrex::Vector<int> bc_hi;
    amrex::Vector<BCRec> bc{1};
    // specify the value on boundaries: direchlet
    amrex::Vector<amrex::Real>  bc_val_phi, bc_val_potential, bc_val_mu, bc_val_solute, bc_val_bcoef; 

};

#endif