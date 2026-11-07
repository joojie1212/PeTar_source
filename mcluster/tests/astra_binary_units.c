#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int get_binaries(int,double**,double,double,int,int*,int,double,double,double,int,int,double,double,int,int,double);
void instar_(void);
extern struct { int ktype[15][15]; } types_;
static int calls;
static int reject_first;
void __wrap_evolv2_(int *kw,double *mass,double *mt,double *r,double *lum,
                   double *mc,double *rc,double *menv,double *renv,double *ospin,
                   double *epoch,double *tms,double *tphys,double *tphysf,
                   double *dtp,double *z,double *zpars,double *tb,double *ecc,
                   double *vkick) {
    (void)kw;(void)r;(void)lum;(void)mc;(void)rc;(void)menv;(void)renv;
    (void)ospin;(void)epoch;(void)tms;(void)tphys;(void)tphysf;(void)dtp;
    (void)z;(void)zpars;(void)ecc;
    assert(fabs(mass[0]-1.)<1e-12 && fabs(mass[1]-1.)<1e-12);
    assert(fabs(*tb-365.25/sqrt(2.))<1e-10);
    mt[0]=mt[1]=0.5; // Keep the period; test inverse conversion with evolved masses.
    vkick[0]=(reject_first && calls==0) ? 1e9 : 0.;
    vkick[1]=0.;
    ++calls;
}
int main(void) {
    instar_();
    assert(types_.ktype[0][0]==1 && types_.ktype[10][10]==15);
    double data[2][15]={{0}}, *star[2]={data[0],data[1]};
    for(int i=0;i<2;++i) { star[i][0]=1./500.; star[i][7]=1.; }
    int n=2;
    srand48(42);
    const double au=1./206264.806;
    get_binaries(1,star,500.,1.,0,&n,0,au,au,1.,0,1,1e-6,0.02,0,0,1e6);
    assert(calls==1 && n==2);
    double r2=0.,v2=0.;
    for(int k=1;k<=3;++k) {
        r2+=pow(star[0][k]-star[1][k],2);
        v2+=pow(star[0][k+3]-star[1][k+3],2);
    }
    const double a=1./(2./sqrt(r2)-v2/(star[0][0]+star[1][0]));
    assert(fabs(a/(au*cbrt(0.5))-1.)<1e-6);
    for(int i=0;i<2;++i) {
        for(int j=0;j<15;++j) star[i][j]=0.;
        star[i][0]=1./500.; star[i][7]=1.;
    }
    calls=0; reject_first=1;
    get_binaries(1,star,500.,1.,0,&n,0,au,au,1.,0,1,1e-6,0.02,1,0,1e6);
    assert(calls==2);
    for(int k=1;k<=6;++k) assert(fabs(star[0][k]+star[1][k])<1e-10);
    puts("mcluster BSE units, collision table and single-component kick rejection passed");
}
