#define HARD_REOPT
#define PROFILE
#include "petar.hpp"

#include <cassert>
#include <fstream>
#include <string>

int main(int argc, char** argv) {
    PeTar::initialFDPS(argc, argv);
    {
        PeTar petar;
        // Match the normal --help cleanup path; FDPS tree destructors expect
        // their internal storage to have been initialized.
        petar.tree_nb.initialize(1, 0.3, 20, 512);
        petar.tree_soft.initialize(1, 0.3, 20, 512);
        petar.input_parameters.print_flag=false;
        petar.input_parameters.fname_snp.value="hard_reopt_logic";
        petar.file_header.nfile=7;
        petar.stat.time=1.25;
        petar.dn_loop=4;
        petar.profile.total.time=10.0;
        petar.profile.hard_isolated.time=9.0;

        // All -1 is disabled even for an otherwise triggering profile.
        petar.checkHardReoptimization();
        assert(!petar.isHardReoptRequested());

        petar.input_parameters.hard_reopt_ratio.value=0.9;
        petar.input_parameters.hard_reopt_windows.value=3;
        petar.input_parameters.hard_reopt_cooldown.value=10;
        petar.hard_reopt_cooldown_remaining=10;

        // A restarted run ignores exactly ten completed profile windows.
        for (int i=0; i<10; ++i) petar.checkHardReoptimization();
        assert(petar.hard_reopt_cooldown_remaining==0);
        assert(petar.hard_reopt_high_window_count==0);
        assert(!petar.isHardReoptRequested());

        // A low window resets the consecutive-high counter.
        petar.checkHardReoptimization();
        petar.checkHardReoptimization();
        assert(petar.hard_reopt_high_window_count==2);
        petar.profile.hard_isolated.time=8.0;
        petar.checkHardReoptimization();
        assert(petar.hard_reopt_high_window_count==0);

        // Three new consecutive high windows create the request.
        petar.profile.hard_isolated.time=9.0;
        petar.checkHardReoptimization();
        petar.checkHardReoptimization();
        assert(!petar.isHardReoptRequested());
        petar.checkHardReoptimization();
        assert(petar.isHardReoptRequested());
        assert(PeTar::getHardReoptExitCode()==42);

        std::ifstream request("hard_reopt_logic.hard_reopt");
        assert(request.good());
        std::string contents((std::istreambuf_iterator<char>(request)),
                             std::istreambuf_iterator<char>());
        assert(contents.find("snapshot=hard_reopt_logic.7")!=std::string::npos);
        assert(contents.find("exit_code=42")!=std::string::npos);
        std::remove("hard_reopt_logic.hard_reopt");
    }
    PeTar::finalizeFDPS();
    return 0;
}
