#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include "bhmerger_output.hpp"
extern "C" ssize_t __real_write(int,const void*,size_t);
static int fault=0;
extern "C" ssize_t __wrap_write(int fd,const void* p,size_t n) {
    if(fault==1) { fault=2; return __real_write(fd,p,n>5?5:n); }
    if(fault==2) { fault=0; errno=EIO; return -1; }
    return __real_write(fd,p,n);
}
std::string contents(const std::filesystem::path& p) {
    std::ifstream f(p); return {std::istreambuf_iterator<char>(f),{}};
}
int main() {
    auto file=std::filesystem::temp_directory_path()/("astra-merger-log-"+std::to_string(getpid()));
    BHMergerOutput out;
    out.open(file.string(),true,true,20,17);
    assert(out.isOpen());
    const auto header=contents(file); assert(header.find("before1.id")!=std::string::npos);
    out.open(file.string(),true,true,20,17); assert(contents(file)==header);
    BHMergerOutput::ParticleState p{};
    fault=1; out.write(0,p,p,0,1,p,20);
    assert(contents(file)==header); // injected partial write was rolled back
    out.write(0,p,p,0,1,p,20); assert(contents(file).size()>header.size());
    out.open((file/"missing").string(),false,true,20,17); assert(!out.isOpen());
    std::filesystem::remove(file);
    out.openChange(file.string(),true,true,17);
    assert(out.isOpen() && contents(file).find("time_old")!=std::string::npos);
    std::filesystem::remove(file);
    // Each rank/thread owns a persistent stream; reopening appends one header.
    struct ChangeParticle { long long id=1; double mass=1.; double pos[3]={}; double vel[3]={}; } cp;
    for (int rank=0;rank<2;++rank) {
        for (int round=0;round<2;++round) {
            out.openChangeSharded(file.string(),round!=0,rank,4,17);
            assert(out.isOpen());
            auto copied=out; // AR interaction copies share the same stream owners
#pragma omp parallel num_threads(4)
            {
                for(int j=0;j<100;++j)
                    copied.writeChange(j,j+1,cp,cp,1,.1,2,.2,3,4);
            }
            out.flushChanges();
        }
        out.close();
        for(int t=0;t<4;++t) {
            const auto shard=file.string()+"."+std::to_string(rank)+"."+std::to_string(t);
            std::ifstream input(shard); std::string line;
            std::getline(input,line); assert(line.find("time_old")==0);
            int rows=0;
            while(std::getline(input,line)) {
                std::istringstream row(line); double value; int columns=0;
                while(row>>value) ++columns;
                assert(columns==19); ++rows;
            }
            assert(rows==200);
            std::filesystem::remove(shard);
        }
    }
}
