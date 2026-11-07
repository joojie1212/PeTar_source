#include <iostream>
#include <fstream>
#include <string>
#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1024
#include "io.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"
#include "static_variables.hpp"


#include <cassert>
#include <cstdio>
#include <sstream>
#include <vector>
#include <unistd.h>

struct TestOrbit {
    double semi=1.0, ecc=0.5, period=1.0;
    void calcOrbit(double) {}
};

int main() {
    const std::string path = "/tmp/petar-change-" + std::to_string(getpid());
    ARInteraction interaction;
    interaction.gravitational_constant = 1.0;
    interaction.bse_manager.vscale = 1.0;
    interaction.change_output.openChange(path, false, true, 17);
    PtclHard p[2];
    for (int i=0; i<2; ++i) {
        p[i].id=i+1;
        p[i].mass=i+1;
        p[i].pos=PS::F64vec(i, 2*i, 3*i);
        p[i].vel=PS::F64vec(4*i, 5*i, 6*i);
    }
    TestOrbit bin;
    // Initial prediction is not a change event, even if requested by caller.
    interaction.refreshBHMergerPrediction(bin, p[0], p[1], 1.0, true);
    const double old_deadline=p[0].time_merger;
    assert(!interaction.isBHMergerOrbitClearlyPerturbed(bin,p[0],p[1]));
    bin.semi=0.5;
    assert(interaction.isBHMergerOrbitClearlyPerturbed(bin,p[0],p[1]));
    interaction.refreshBHMergerPrediction(bin,p[0],p[1],2.0,true);
    assert(p[0].merger_record_time==2.0);
    PtclHard copied;
    copied.DataCopy(p[0]);
    assert(copied.merger_record_time==2.0);
    // A replacement companion must not produce a false before/after pair.
    p[1].id=3;
    interaction.refreshBHMergerPrediction(bin,p[0],p[1],3.0,true);
    interaction.change_output.close();
    interaction.change_output.openChange(path,true,true,17);
    bin.ecc=0.9;
    interaction.refreshBHMergerPrediction(bin,p[0],p[1],4.0,true);
    // An unbound orbit records the invalidated deadline as -1.
    bin.semi=-1.0;
    interaction.refreshBHMergerPrediction(bin,p[0],p[1],5.0,true);
    assert(p[0].merger_record_time==-1.0);
    interaction.change_output.close();
    std::ifstream input(path);
    std::string line;
    std::getline(input,line);
    assert(line.find("time_old time_new id1 id2")==0);
    int rows=0;
    while (std::getline(input,line)) {
        std::istringstream stream(line);
        std::vector<double> values;
        double value;
        while(stream>>value) values.push_back(value);
        assert(values.size()==18);
        if (rows==0) {
            assert(values[0]==1.0 && values[1]==2.0);
            assert(values[2]==1 && values[3]==2);
            assert(values[4]==1 && values[5]==2);
            assert(values[6]==1.0 && values[7]==0.5);
            assert(values[8]==0.5 && values[9]==0.5);
            assert(values[10]==old_deadline);
            assert(values[11]<old_deadline);
        }
        if(rows==1) assert(values[0]==3 && values[1]==4 && values[3]==3);
        if(rows==2) assert(values[8]==-1 && values[11]==-1);
        for(int k=0;k<6;++k) assert(values[12+k]==-(k+1));
        ++rows;
    }
    assert(rows==3);
    std::remove(path.c_str());
    std::cout << "orbit change regression passed\n";
}
