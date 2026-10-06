/* Numerical regression tests independent of McLuster's legacy main. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../darkmatter.h"
#define CHECK(x) do { if (!(x)) { fprintf(stderr,"FAILED line %d: %s\n",__LINE__,#x); exit(1); } } while (0)

static void test_shell(void) {
    long double k,d,kp,dp,km,dm;
    dm_shell(0.7,2.1,0.3,&k,&d);
    dm_shell(0.700001,2.1,0.3,&kp,&dp);
    dm_shell(0.699999,2.1,0.3,&km,&dm);
    CHECK(fabsl(d-(kp-km)/0.000002L)<1e-10L);
    long double numerical=0;
    const int n=100000;
    for (int j=0;j<n;++j) {
        long double mu=-1+2*(j+0.5L)/n;
        numerical+=1/sqrtl(0.7L*0.7L+2.1L*2.1L-2*0.7L*2.1L*mu+0.3L*0.3L)/n;
    }
    CHECK(fabsl(k-numerical)<1e-10L);
}

static void test_plummer(void) {
    DMOptions o={0};
    o.grid=2048; o.rs=2; o.rt=20; o.eps=0.3;
    DMModel m={0};
    double rh=1/sqrt(pow(2.0,2.0/3.0)-1);
    dm_build(&m,&o,1,rh); /* DM mass zero: analytic self-consistent Plummer limit */
    dm_invert(&m,0); dm_speeds(&m,0);
    double maxerr=0;
    for (int i=0;i<m.n;++i) if (m.r[i]>0.02 && m.r[i]<100) {
        double analytic=24*sqrt(2.0)/(7*pow(DM_PI,3)*pow(DM_G,5))
                        *pow(m.psi[0][i],3.5);
        double err=fabs(m.df[0][i]/analytic-1);
        if (err>maxerr) maxerr=err;
    }
    printf("analytic Plummer DF max relative error: %.6g\n",maxerr);
    CHECK(maxerr<0.01);
    dm_free(&m);
}
static void test_combined(void) {
    DMOptions o={0};
    o.grid=2048; o.rs=10; o.rt=100; o.eps=1; o.mass=10000;
    DMModel m={0}, heavy={0};
    dm_build(&m,&o,1024,1);
    dm_build(&heavy,&o,2048,1);
    for (int c=0;c<2;++c) {
        CHECK(heavy.psi[c][m.n/2]>m.psi[c][m.n/2]);
        dm_invert(&m,c); dm_speeds(&m,c);
        for (int i=0;i<m.n;++i) CHECK(m.df[c][i]>=0 && isfinite(m.df[c][i]));
        /* Compare the sampled second moment with the isotropic Jeans integral.
         * This tests velocities, not merely positivity or density normalization.
         */
        for (int i=0;i<m.n-2;++i) if (m.r[i]>0.2 && m.r[i]<20 && i%31==0) {
            double pressure=0;
            for (int j=i;j<m.n-1;++j) {
                double r=m.r[j], r1=m.r[j+1], s0,s1;
                if (c==0) {
                    double u=r/m.a, v=r1/m.a;
                    s0=-5*u*u/(1+u*u); s1=-5*v*v/(1+v*v);
                } else {
                    double x=r/o.rs,y=r/o.rt, x1=r1/o.rs,y1=r1/o.rt;
                    s0=-1-2*x/(1+x)-4*y*y/(1+y*y);
                    s1=-1-2*x1/(1+x1)-4*y1*y1/(1+y1*y1);
                }
                pressure-=0.5*m.step*(m.rho[c][j]*m.rho[c][j]*s0/m.q[c][j]
                    +m.rho[c][j+1]*m.rho[c][j+1]*s1/m.q[c][j+1]);
            }
            double *cdf=m.speed[c]+(size_t)i*(DM_NV+1), v2=0;
            for (int k=1;k<=DM_NV;++k) {
                double x=(dm_speed_fraction(k)+dm_speed_fraction(k-1))/2;
                v2+=(cdf[k]-cdf[k-1])*2*m.psi[c][i]*x*x;
            }
            CHECK(fabs(v2/(3*pressure/m.rho[c][i])-1)<0.025);
        }
    }
    /* Shared mass affects DM velocities: compare the local second moment. */
    dm_invert(&heavy,1); dm_speeds(&heavy,1);
    int i=dm_bracket(m.r,m.n,1);
    double moment[2]={0};
    DMModel *models[2]={&m,&heavy};
    for (int a=0;a<2;++a) for (int k=1;k<=DM_NV;++k) {
        double x=dm_speed_fraction(k), *cdf=models[a]->speed[1]+(size_t)i*(DM_NV+1);
        moment[a]+=(cdf[k]-cdf[k-1])*x*x*2*models[a]->psi[1][i];
    }
    CHECK(moment[1]>moment[0]*1.15);
    dm_free(&m); dm_free(&heavy);
}
int main(void) {
    test_shell(); test_plummer(); test_combined();
    puts("darkmatter numerical tests passed");
    return 0;
}
