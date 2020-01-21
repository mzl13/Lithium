#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include "Lithium.H"

using namespace amrex;

int main(int argc, char* argv[])
{

    amrex::system::verbose = 0;
    amrex::Initialize(argc,argv);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const double strt_time = amrex::second();

    {

        Lithium li;

        const double cal_time = amrex::second() - strt_time;

        Print()<< "--> Total Time: "<< cal_time << " seconds\n";
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
  
    return 0;
}
