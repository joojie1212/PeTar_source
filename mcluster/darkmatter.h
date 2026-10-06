/* Optional, isolated spherical NFW + Plummer initial-condition generator.
 * All masses: Msun; lengths: pc; velocities: pc/Myr.
 * No calls into the legacy position/velocity or virial-rescaling routines.
 */
#ifndef MCLUSTER_DARKMATTER_H
#define MCLUSTER_DARKMATTER_H
#ifdef DARKMATTER
#include <errno.h>
#include <limits.h>

#define DM_PI 3.141592653589793238462643383279502884
#define DM_G 0.00449830997959438
#define DM_NV 384

typedef struct {
    int enabled, specified, number, grid, bse;
    double mass, rs, rt, eps;
} DMOptions;

typedef struct {
    int n;
    double mass, a, norm, step;
    double *r, *rho[2], *psi[2], *q[2], *df[2], *cdf[2], *speed[2];
    double max_density_error[2];
} DMModel;

static void dm_fail(const char *message) {
    fprintf(stderr, "NFW initial conditions: %s\n", message);
    exit(EXIT_FAILURE);
}
static void *dm_alloc(size_t n, size_t size) {
    void *p = calloc(n, size);
    if (!p) dm_fail("memory allocation failed");
    return p;
}
static double dm_positive(const char *s) {
    char *end;
    errno = 0;
    double v = strtod(s, &end);
    if (errno || end == s || *end || !isfinite(v) || !(v > 0))
        dm_fail("expected a finite positive numeric parameter");
    return v;
}
static int dm_integer(const char *s) {
    char *end;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno || end == s || *end || v < 1 || v > INT_MAX)
        dm_fail("expected a positive integer parameter");
    return (int)v;
}
static void dm_validate(DMOptions *o) {
    if (!o->number || !o->mass || !o->rs || !o->rt || !o->eps)
        dm_fail("--dm-nfw requires --dm-number, --dm-mass, --dm-rs, --dm-rt, --dm-softening");
    if (o->rt < o->rs) dm_fail("--dm-rt must be >= --dm-rs");
    if (!o->grid) o->grid = 2048;
    if (o->grid < 512 || o->grid > 8192)
        dm_fail("--dm-grid must lie between 512 and 8192");
}
static void dm_help(void) {
    printf("\nOptional NFW + Plummer mode (DARKMATTER build):\n"
           "  --dm-nfw                  generate paired PeTar ASCII initial conditions\n"
           "  --dm-number N             equal-mass dark matter particle count\n"
           "  --dm-mass M               total tapered halo mass [Msun], NOT M200\n"
           "  --dm-rs R                 NFW scale radius [pc]\n"
           "  --dm-rt R                 smooth taper scale [pc], >= rs\n"
           "  --dm-softening E          PeTar DM Plummer softening [pc], >0\n"
           "  --dm-grid N               radial/energy grid, 512..8192 (default 2048)\n"
           "  --petar-interrupt none|bse stellar row format (default none)\n"
           "  Requires isolated unsegregated Plummer: -P 0 -t 0 -S 0 -D 3,\n"
           "  no binaries/gas/evolved population; -u 1 -C 3 -Q 0.5.\n"
           "  Writes PREFIX.star, PREFIX.dm, PREFIX.dm.profile, PREFIX.dm.info.\n");
}
static double dm_shape(double r, const DMOptions *o) {
    double x = r/o->rs, y = r/o->rt;
    return 1.0/(x*(1+x)*(1+x)*(1+y*y)*(1+y*y));
}
static double dm_star_rho(double r, const DMModel *m) {
    double x = r/m->a;
    return 3*m->mass/(4*DM_PI*m->a*m->a*m->a)*pow(1+x*x, -2.5);
}
/* Angular average of the Plummer kernel is 2/(A+B).  Unlike the
 * equivalent (A-B)/(2rs), this has no subtraction of close roots.
 * Long double keeps the small central force accurate.
 */
static void dm_shell(long double r, long double s, long double eps,
                     long double *kernel, long double *derivative) {
    long double A = sqrtl((r+s)*(r+s)+eps*eps);
    long double B = sqrtl((r-s)*(r-s)+eps*eps);
    *kernel = 2/(A+B);
    *derivative = -2*((r+s)/A+(r-s)/B)/((A+B)*(A+B));
}
static void dm_build(DMModel *m, const DMOptions *o, double mass, double rh) {
    m->n = o->grid;
    m->mass = mass;
    m->a = rh*sqrt(pow(2.0, 2.0/3.0)-1.0);
    double small = fmin(m->a, fmin(o->rs, o->eps));
    double large = fmax(m->a, o->rt);
    if (!isfinite(mass) || !(mass > 0) || !isfinite(rh) || !(rh > 0))
        dm_fail("stellar mass and half-mass radius must be finite and positive");
    if (large/small > 1e6)
        dm_fail("scale ratio exceeds 1e6; choose less extreme scales");
    double rmin = 1e-3*small, rmax = 1e4*large;
    m->step = log(rmax/rmin)/(m->n-1);
    if (m->step > 0.065)
        dm_fail("radial grid too coarse for these scales; increase --dm-grid");
    m->r = dm_alloc(m->n, sizeof(double));
    for (int c=0; c<2; ++c) {
        m->rho[c] = dm_alloc(m->n, sizeof(double));
        m->psi[c] = dm_alloc(m->n, sizeof(double));
        m->q[c] = dm_alloc(m->n, sizeof(double));
        m->df[c] = dm_alloc(m->n, sizeof(double));
        m->cdf[c] = dm_alloc(m->n, sizeof(double));
        m->speed[c] = dm_alloc((size_t)m->n*(DM_NV+1), sizeof(double));
    }
    double *sr = dm_alloc(m->n-1, sizeof(double));
    double *ms = dm_alloc(m->n-1, sizeof(double));
    double *md = dm_alloc(m->n-1, sizeof(double));
    double sumd = 0, sums = 0;
    for (int i=0; i<m->n; ++i) m->r[i] = rmin*exp(i*m->step);
    for (int i=0; i<m->n-1; ++i) {
        double r = sr[i] = sqrt(m->r[i]*m->r[i+1]);
        md[i] = 4*DM_PI*r*r*r*dm_shape(r,o)*m->step;
        ms[i] = 4*DM_PI*r*r*r*dm_star_rho(r,m)*m->step;
        sumd += md[i]; sums += ms[i];
    }
    m->norm = o->mass/sumd;
    for (int i=0; i<m->n-1; ++i) {
        md[i] *= m->norm;
        ms[i] *= mass/sums;
        m->cdf[0][i+1] = m->cdf[0][i]+ms[i]/mass;
        m->cdf[1][i+1] = m->cdf[1][i]+md[i]/o->mass;
    }
    m->cdf[0][m->n-1] = m->cdf[1][m->n-1] = 1;
    for (int i=0; i<m->n; ++i) {
        double r = m->r[i], u = r/m->a, x = r/o->rs, y = r/o->rt;
        m->rho[0][i] = dm_star_rho(r,m);
        m->rho[1][i] = m->norm*dm_shape(r,o);
        long double pd = 0, ps = 0, dd = 0, ds = 0;
        for (int j=0; j<m->n-1; ++j) {
            long double k, d;
            dm_shell(r, sr[j], o->eps, &k, &d);
            pd += md[j]*k; dd += md[j]*d;
            ps += ms[j]*k; ds += ms[j]*d;
        }
        /* Match PeTar: stars see Newtonian stars + softened DM;
         * DM sees softened stars + softened DM. Both contain BOTH masses.
         */
        m->psi[0][i] = DM_G*(double)(pd+mass/sqrt(r*r+m->a*m->a));
        m->psi[1][i] = DM_G*(double)(pd+ps);
        double deriv[2] = {
            DM_G*(double)(dd-mass*r/pow(r*r+m->a*m->a,1.5)),
            DM_G*(double)(dd+ds)
        };
        double slope[2] = {-5*u*u/(1+u*u),
                          -1-2*x/(1+x)-4*y*y/(1+y*y)};
        for (int c=0; c<2; ++c) {
            if (!(deriv[c] < 0) || !isfinite(m->psi[c][i]))
                dm_fail("non-monotone or non-finite total potential");
            m->q[c][i] = m->rho[c][i]*slope[c]/(r*deriv[c]);
            if (i && !(m->psi[c][i] < m->psi[c][i-1]))
                dm_fail("unresolved central potential; reduce scale contrast");
        }
    }
    free(sr); free(ms); free(md);
}
/* Eddington: f(E) = 1/(sqrt(8) pi^2) int_0^E q'(Psi)/sqrt(E-Psi)dPsi,
 * q = d rho/d Psi. Integrate each linear q segment analytically.
 * The outer q(0)=0 boundary is included, not silently dropped.
 */
static void dm_invert(DMModel *m, int c) {
    int n = m->n;
    for (int i=0; i<n; ++i) {
        double E = m->psi[c][i];
        long double sum = 0, absolute = 0;
        for (int j=n-1; j>=i; --j) {
            double lo = j==n-1 ? 0 : m->psi[c][j+1];
            double hi = m->psi[c][j];
            double qlo = j==n-1 ? 0 : m->q[c][j+1];
            double slope = (m->q[c][j]-qlo)/(hi-lo);
            /* Rationalized sqrt difference removes far-bin cancellation. */
            long double term = 2*slope*(hi-lo)/
                (sqrt(E-lo)+sqrt(fmax(0,E-hi)));
            sum += term; absolute += fabsl(term);
        }
        if (sum < -1e-8L*absolute) {
            fprintf(stderr, "NFW: negative %s DF at r=%.9g pc, E=%.9g: %.9Lg.\n",
                    c ? "DM" : "stellar", m->r[i], E, sum);
            dm_fail("no nonnegative isotropic DF for this density/potential pair; "
                    "try a finer grid to distinguish numerical error, or change physical parameters");
        }
        m->df[c][i] = (double)fmaxl(0,sum)/(sqrt(8.0)*DM_PI*DM_PI);
    }
}
static int dm_bracket(const double *x, int n, double v) {
    int lo = 0, hi = n-1;
    while (hi-lo > 1) {
        int mid = (hi+lo)/2;
        if (x[mid] <= v) lo = mid; else hi = mid;
    }
    return lo;
}
static double dm_df(const DMModel *m, int c, double E) {
    if (!(E > 0)) return 0;
    int n = m->n;
    if (E <= m->psi[c][n-1])
        return m->df[c][n-1]*sqrt(E/m->psi[c][n-1]);
    int lo = 0, hi = n-1;
    while (hi-lo > 1) {
        int k = (lo+hi)/2;
        if (m->psi[c][k] >= E) lo = k; else hi = k;
    }
    double t = (E-m->psi[c][hi])/(m->psi[c][lo]-m->psi[c][hi]);
    t = fmax(0,fmin(1,t));
    double a = m->df[c][hi], b = m->df[c][lo];
    return a>0 && b>0 ? exp((1-t)*log(a)+t*log(b)) : (1-t)*a+t*b;
}
/* Speed grid clustered at v=0 for the cold, softened NFW cusp.
 * Store a conditional CDF at every radius, avoiding a guessed rejection bound.
 */
static double dm_speed_fraction(int k) {
    return k ? exp(-16.0*(1.0-(double)k/DM_NV)) : 0;
}
static void dm_speeds(DMModel *m, int c) {
    for (int i=0; i<m->n; ++i) {
        double *cdf = m->speed[c]+(size_t)i*(DM_NV+1);
        double psi = m->psi[c][i], previous = 0, xv = 0;
        for (int k=1; k<=DM_NV; ++k) {
            double x = dm_speed_fraction(k);
            double value = x*x*dm_df(m,c,psi*(1-x*x));
            cdf[k] = cdf[k-1]+0.5*(value+previous)*(x-xv);
            previous=value; xv=x;
        }
        if (!(cdf[DM_NV]>0) || !isfinite(cdf[DM_NV]))
            dm_fail("invalid conditional speed distribution");
        double reconstructed = 4*DM_PI*pow(2*psi,1.5)*cdf[DM_NV];
        /* Ignore unresolved inner boundary and the negligible outer tail. */
        if (m->cdf[c][i]>1e-5 && m->cdf[c][i]<1-1e-5) {
            double error = fabs(reconstructed/m->rho[c][i]-1);
            if (error>m->max_density_error[c]) m->max_density_error[c]=error;
        }
        for (int k=1; k<=DM_NV; ++k) cdf[k]/=cdf[DM_NV];
    }
    printf("NFW %s DF density reconstruction: max relative error %.5g\n",
           c ? "DM" : "stellar",m->max_density_error[c]);
    if (m->max_density_error[c]>0.03)
        dm_fail("DF density reconstruction error exceeds 3%; increase --dm-grid");
}
static void dm_direction(double length, double *out) {
    double z=2*drand48()-1, az=2*DM_PI*drand48();
    double xy=length*sqrt(fmax(0,1-z*z));
    out[0]=xy*cos(az); out[1]=xy*sin(az); out[2]=length*z;
}
static void dm_sample(const DMModel *m, int c, double *row) {
    double u=drand48();
    int i=dm_bracket(m->cdf[c],m->n,u);
    double t=(u-m->cdf[c][i])/(m->cdf[c][i+1]-m->cdf[c][i]);
    double r=m->r[i]*exp(t*m->step);
    dm_direction(r,row+1);
    /* Mixture of the two adjacent normalized CDFs; same energy convention. */
    double cdf[DM_NV+1];
    for (int k=0;k<=DM_NV;++k)
        cdf[k]=(1-t)*m->speed[c][(size_t)i*(DM_NV+1)+k]
              +t*m->speed[c][(size_t)(i+1)*(DM_NV+1)+k];
    double w=drand48();
    int k=dm_bracket(cdf,DM_NV+1,w);
    double frac=(w-cdf[k])/(cdf[k+1]-cdf[k]);
    double x=dm_speed_fraction(k)+frac*(dm_speed_fraction(k+1)-dm_speed_fraction(k));
    double psi=(1-t)*m->psi[c][i]+t*m->psi[c][i+1];
    dm_direction(x*sqrt(2*psi),row+4);
}
static FILE *dm_open(const char *prefix, const char *suffix) {
    size_t n=strlen(prefix)+strlen(suffix)+1;
    char *path=dm_alloc(n,1);
    snprintf(path,n,"%s%s",prefix,suffix);
    FILE *fp=fopen(path,"w");
    if (!fp) { perror(path); free(path); dm_fail("cannot open output file"); }
    free(path);
    return fp;
}
static void dm_close(FILE *fp) {
    int bad=ferror(fp);
    if (fclose(fp) || bad) dm_fail("failed to write output file");
}
static void dm_free(DMModel *m) {
    free(m->r);
    for (int c=0;c<2;++c) {
        free(m->rho[c]); free(m->psi[c]); free(m->q[c]);
        free(m->df[c]); free(m->cdf[c]); free(m->speed[c]);
    }
}
static int dm_generate(const DMOptions *o, int nstar, double **star,
                       double rh, const char *prefix, int blackhole, unsigned seed) {
    if (nstar<2) dm_fail("at least two stars are required");
    double mass=0;
    for (int i=0;i<nstar;++i) mass+=star[i][0];
    DMModel m={0};
    printf("Building NFW + Plummer model in combined potentials (grid=%d).\n",o->grid);
    dm_build(&m,o,mass,rh);
    for (int c=0;c<2;++c) { dm_invert(&m,c); dm_speeds(&m,c); }
    double (*dark)[7]=dm_alloc(o->number,sizeof(*dark));
    for (int i=0;i<nstar;++i) dm_sample(&m,0,star[i]);
    for (int i=0;i<o->number;++i) {
        dark[i][0]=o->mass/o->number;
        dm_sample(&m,1,dark[i]);
    }
    /* One common translation and boost only. Independent recentring would
     * move the two components relative to their shared potential centre.
     */
    double cm[6]={0};
    for (int i=0;i<nstar;++i)
        for (int j=0;j<6;++j) cm[j]+=star[i][0]*star[i][j+1];
    for (int i=0;i<o->number;++i)
        for (int j=0;j<6;++j) cm[j]+=dark[i][0]*dark[i][j+1];
    for (int j=0;j<6;++j) cm[j]/=mass+o->mass;
    FILE *fs=dm_open(prefix,".star"), *fd=dm_open(prefix,".dm");
    fprintf(fs,"0 %d 0\n",nstar);
    fprintf(fd,"0 %d 0 %.17g\n",o->number,o->eps);
    for (int i=0;i<nstar;++i) {
        fprintf(fs,"%.17g",star[i][0]);
        for (int j=0;j<6;++j) fprintf(fs," %.17g",star[i][j+1]-cm[j]);
        fprintf(fs," 0"); /* binary_state */
        if (o->bse) {
            /* New, unevolved stars: kw=1 (mass>=0.7) or kw=0.
             * PeTar/BSE initializes radii on its first stellar evolution call.
             */
            int kw=blackhole ? 14 : (star[i][0]>=0.7 ? 1 : 0);
            double core=blackhole ? star[i][0] : 0;
            fprintf(fs," 0 0 0 0 %d %.17g %.17g 0 %.17g 0 0 0 0 0",
                    kw,star[i][0],star[i][0],core);
        }
        /* r_search id group[2] changeover[2] acc[3] pot_tot pot_soft n_ngb */
        fprintf(fs," 0 %d 0 0 0 0 0 0 0 0 0 0\n",i+1);
    }
    for (int i=0;i<o->number;++i) {
        fprintf(fd,"%.17g",dark[i][0]);
        for (int j=0;j<6;++j) fprintf(fd," %.17g",dark[i][j+1]-cm[j]);
        fprintf(fd," %d 0 0 0 0 0\n",i+1);
    }
    dm_close(fs); dm_close(fd);
    FILE *profile=dm_open(prefix,".dm.profile");
    fprintf(profile,"# r_pc rho_star rho_dm Psi_star Psi_dm f_star f_dm CDF_star CDF_dm\n");
    for (int i=0;i<m.n;++i)
        fprintf(profile,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                m.r[i],m.rho[0][i],m.rho[1][i],m.psi[0][i],m.psi[1][i],
                m.df[0][i],m.df[1][i],m.cdf[0][i],m.cdf[1][i]);
    dm_close(profile);
    FILE *info=dm_open(prefix,".dm.info");
    fprintf(info,"model=NFW/(1+(r/rt)^2)^2 + Plummer\n"
            "units=Msun pc Myr pc/Myr\nG=%.17g\nseed=%u\n"
            "Nstar=%d\nMstar=%.17g\nRh=%.17g\nNdm=%d\nMdm_total=%.17g\n"
            "rs=%.17g\nrt=%.17g\nsoftening=%.17g\ngrid=%d\npetar_interrupt=%s\n"
            "DF_density_error_star=%.17g\nDF_density_error_dm=%.17g\n",
            DM_G,seed,nstar,mass,rh,o->number,o->mass,o->rs,o->rt,o->eps,o->grid,
            o->bse ? "bse" : "none",m.max_density_error[0],m.max_density_error[1]);
    fprintf(info,"subtracted_common_position=%.17g %.17g %.17g\n"
            "subtracted_common_velocity=%.17g %.17g %.17g\n",cm[0],cm[1],cm[2],cm[3],cm[4],cm[5]);
    dm_close(info);
    printf("PeTar inputs: %s.star and %s.dm (use -u 1 -i 1; interrupt=%s).\n",
           prefix,prefix,o->bse ? "bse" : "none");
    free(dark); dm_free(&m);
    return 0;
}
#endif /* DARKMATTER */
#endif /* MCLUSTER_DARKMATTER_H */
