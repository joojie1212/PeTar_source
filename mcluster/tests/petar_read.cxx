// Optional interoperability check using PeTar's actual ASCII readers.
#include "particle_simulator.hpp"
#include "io.hpp"
#include "soft_ptcl.hpp"
#include "darkmatter.hpp"
int main(int argc, char **argv) {
    if (argc != 3) return 1;
    FILE *fs=fopen(argv[1],"r"), *fd=fopen(argv[2],"r");
    if (!fs || !fd) return 2;
    FileHeader sh;
    DarkMatterFileHeader dh;
    sh.readAscii(fs); dh.readAscii(fd);
    if (sh.nfile!=dh.nfile || sh.time!=dh.time || !(dh.eps>0)) return 3;
    for (long long i=1;i<=sh.n_body;++i) {
        FPSoft s; s.readAscii(fs);
        if (s.id!=i || s.star.kw!=14 || s.star.mc!=s.mass || s.star.mt!=s.mass) return 4;
    }
    for (long long i=1;i<=dh.n_body;++i) {
        DarkMatterParticle d; d.readAscii(fd);
        if (d.id!=i) return 5;
    }
    int c;
    while ((c=fgetc(fs))!=EOF) if (!isspace(c)) return 6;
    while ((c=fgetc(fd))!=EOF) if (!isspace(c)) return 7;
    fclose(fs); fclose(fd);
    return 0;
}
