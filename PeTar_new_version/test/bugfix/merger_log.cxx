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
}
