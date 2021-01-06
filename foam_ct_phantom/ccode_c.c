// -----------------------------------------------------------------------
// Copyright 2019 Centrum Wiskunde & Informatica, Amsterdam

// Author: Daniel M. Pelt
// Contact: D.M.Pelt@cwi.nl
// Website: http://dmpelt.github.io/foam_ct_phantom/
// License: MIT

// This file is part of foam_ct_phantom, a Python package for generating
// foam-like phantoms for CT.
// -----------------------------------------------------------------------

#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#ifdef _MSC_VER

#define DECLDIR __declspec(dllexport)

#else

#define DECLDIR

#endif

#ifndef __MTWISTER_H
#define __MTWISTER_H

#define STATE_VECTOR_LENGTH 624
#define STATE_VECTOR_M      397 /* changes to STATE_VECTOR_LENGTH also require changes to this */

typedef struct tagMTRand {
  unsigned long mt[STATE_VECTOR_LENGTH];
  int index;
} MTRand;

MTRand seedRand(unsigned long seed);
unsigned long genRandLong(MTRand* rand);
float genRand(MTRand* rand);
void m_seedRand(MTRand* rand, unsigned long seed);

#endif /* #ifndef __MTWISTER_H */


// OpenMP set number of threads
DECLDIR void set_threads(const unsigned int nthrd){
    omp_set_num_threads(nthrd);
}

MTRand randgen;
DECLDIR void setseed(const unsigned int seed){
    m_seedRand(&randgen, seed);
}

DECLDIR void drawnewpositions(float * const pos3, float * const ds, const unsigned int ntrials, const float zrange){
    for(unsigned int i=0; i<ntrials; i++){
        float x,y,z;
        while(1){
            x = 2*genRand(&randgen)-1;
            y = 2*genRand(&randgen)-1;
            if(x*x + y*y<=1){
                z = zrange*(2*genRand(&randgen)-1);
                break;
            }
        }
        pos3[3*i] = x;
        pos3[3*i+1] = y;
        pos3[3*i+2] = z;
        ds[i] = sqrtf(x*x + y*y)-1;
    }
}

DECLDIR unsigned int newsphere(float * const pos3, float * const ds, const float * const spheres, const unsigned int ntrials, const unsigned int nspheres, const float zrange, unsigned int * const updated){
    const float nx = spheres[(nspheres-1)*5];
    const float ny = spheres[(nspheres-1)*5+1];
    const float nz = spheres[(nspheres-1)*5+2];
    const float nsz = spheres[(nspheres-1)*5+3];
    
    const unsigned int nthr = omp_get_max_threads();
    unsigned int * seeds = (unsigned int *) malloc(sizeof(unsigned int)*nthr);
    int j;
    for(j=0;j<nthr;j++){
        seeds[j] = (unsigned int) UINT_MAX*genRand(&randgen);
    }

    unsigned int nupdated = 0;
    unsigned int * tobeplaced = (unsigned int *) malloc(sizeof(unsigned int) * ntrials);
    unsigned int ntobeplaced = 0;

    int i;
    #pragma omp parallel for shared(nupdated, ntobeplaced) private(i)
    for(i=0; i<ntrials; i++){
        float x = pos3[3*i];
        float y = pos3[3*i+1];
        float z = pos3[3*i+2];
        float dsv = sqrtf((nx-x)*(nx-x) + (ny-y)*(ny-y) + (nz-z)*(nz-z)) - nsz;
        if(dsv<0){
           #pragma omp critical
            {
                tobeplaced[ntobeplaced] = i;
                ntobeplaced++;
            }
        }else{
            if(dsv<-ds[i]){
                ds[i] = -dsv;
                #pragma omp critical
                {
                    updated[nupdated] = i;
                    nupdated++;
                }
            }
        }
    }

    volatile unsigned int cur=0;
    #pragma omp parallel shared(nupdated, cur) firstprivate(ntobeplaced)
    {
        MTRand lrandgen;
        m_seedRand(&lrandgen, seeds[omp_get_thread_num()]);
        float dsv, x, y, z;
        char cont = 1;
        while(cont){
            do{
                if (cur>=ntobeplaced) break;
                x = 2*genRand(&lrandgen)-1;
                y = 2*genRand(&lrandgen)-1;
                dsv = 1 - sqrtf(x*x + y*y);
                if(dsv<0){
                    continue;
                }
                z = zrange*(2*genRand(&lrandgen)-1);
                int j;
                for(j=0; j<nspheres; j++){
                    const float dsn = sqrtf((spheres[5*j]-x)*(spheres[5*j]-x) + (spheres[5*j+1]-y)*(spheres[5*j+1]-y) + (spheres[5*j+2]-z)*(spheres[5*j+2]-z)) - spheres[5*j+3];
                    if(dsn<dsv){
                        dsv = dsn;
                        if (dsv<0) break;
                    }
                }
            }while(dsv<0);

            #pragma omp critical
            {
                if (cur>=ntobeplaced){
                    cont=0;
                }else{
                    int k = tobeplaced[cur];
                    pos3[3*k] = x;
                    pos3[3*k+1] = y;
                    pos3[3*k+2] = z;
                    ds[k] = -dsv;
                    updated[nupdated] = k;
                    nupdated++;
                    cur++;
                }
            }
        }
    }
    free(seeds);
    return nupdated;
}

// DECLDIR unsigned int newsphere(float * const pos3, float * const ds, const float * const spheres, const unsigned int ntrials, const unsigned int nspheres, const float zrange, unsigned int * const updated){
//     const float nx = spheres[(nspheres-1)*5];
//     const float ny = spheres[(nspheres-1)*5+1];
//     const float nz = spheres[(nspheres-1)*5+2];
//     const float nsz = spheres[(nspheres-1)*5+3];
//     unsigned int nupdated = 0;
//     unsigned int * tobeplaced = (unsigned int *) malloc(sizeof(unsigned int) * ntrials);
//     unsigned int ntobeplaced = 0;
//     int i;
//     #pragma omp parallel for shared(nupdated, ntobeplaced) private(i)
//     for(i=0; i<ntrials; i++){
//         float x = pos3[3*i];
//         float y = pos3[3*i+1];
//         float z = pos3[3*i+2];
//         float dsv = sqrtf((nx-x)*(nx-x) + (ny-y)*(ny-y) + (nz-z)*(nz-z)) - nsz;
//         if(dsv<0){
//            #pragma omp critical
//             {
//                 tobeplaced[ntobeplaced] = i;
//                 ntobeplaced++;
//             }
//         }else{
//             if(dsv<-ds[i]){
//                 ds[i] = -dsv;
//                 #pragma omp critical
//                 {
//                     updated[nupdated] = i;
//                     nupdated++;
//                 }
//             }
//         }
//     }

//     unsigned int j=0;
//     #pragma omp parallel shared(nupdated, j)
//     while(1){
//         float x, y, z;
        
//     }


//     free(tobeplaced);
//     return nupdated;
// }



DECLDIR void genvol(const float * const spheres, const unsigned int nspheres, float * const vol, const unsigned int * const n, const float voxsize, const float * const c, const unsigned int supersampling, const unsigned int iz){
    const unsigned int ninslc = n[0]*n[1];
    int i;
    #pragma omp parallel for private(i)
    for(i=0; i<ninslc; i++){
        vol[iz*ninslc+i] = 0;
    }
    for(unsigned int sz=0; sz<supersampling; sz++){                      
        #pragma omp parallel
        {
            unsigned int * const temp = (unsigned int *) malloc(nspheres*sizeof(unsigned int));
            long double * const dzs = (long double *) malloc(nspheres*sizeof(long double));
            long double * const s2s = (long double *) malloc(nspheres*sizeof(long double));
            for(unsigned int j=0; j<nspheres; j++){
                s2s[j] = spheres[5*j+3]*spheres[5*j+3];
            }
            long double z = c[2] + (iz+0.5)*voxsize - n[2]*voxsize/2 + sz*voxsize/supersampling - voxsize/2 + voxsize/(2*supersampling);
            unsigned int ntocheck=0;
            for(unsigned int j=0; j<nspheres; j++){
                dzs[j] = (z-spheres[5*j+2])*(z-spheres[5*j+2]);
                if(dzs[j]<s2s[j]){
                    temp[ntocheck] = j;
                    ntocheck++;
                }
            }
            #pragma omp for schedule(dynamic) private(i)
            for(i=0; i<ninslc; i++){
                const unsigned int iy = i/n[0];
                const unsigned int ix = i%n[0];
                for(unsigned int sx=0; sx<supersampling; sx++){
                    const long double x = c[0] + (ix+0.5)*voxsize - n[0]*voxsize/2 + sx*voxsize/supersampling - voxsize/2 + voxsize/(2*supersampling);
                    for(unsigned int sy=0; sy<supersampling; sy++){
                        const long double y = c[1] + (iy+0.5)*voxsize - n[1]*voxsize/2 + sy*voxsize/supersampling - voxsize/2 + voxsize/(2*supersampling);
                        if(sqrtl(x*x+y*y)>1){
                            continue;
                        }                    
                        unsigned char found=0;
                        long double s2,dx,dy,dz;
                        for(unsigned int q=0; q<ntocheck; q++){
                            const unsigned int j = temp[q];
                            dx = (x-spheres[5*j])*(x-spheres[5*j]);
                            if(dx+dzs[j] >= s2s[j]) continue;
                            dy = (y-spheres[5*j+1])*(y-spheres[5*j+1]);
                            if(dx+dy+dzs[j]<s2s[j]){
                                vol[iz*ninslc+i]+=spheres[5*j+4];
                                found=1;
                                break;
                            }
                        }
                        if(found==0){
                            vol[iz*ninslc+i]+=1;
                        }
                    }
                }
            }
            free(dzs);
            free(temp);
            free(s2s);
        }
    }
    #pragma omp parallel for private(i)
    for(i=0; i<ninslc; i++){
        vol[iz*ninslc+i]/=supersampling*supersampling*supersampling;
    }
}

DECLDIR void genparproj(const float * const spheres, const unsigned int nspheres, float * const proj, const unsigned int * const n, const float pixsize, const float * const c, const float angle, const float * const rotc){
    const unsigned int ntotal = n[0]*n[1];
    const long double ca = cosl(angle);
    const long double sa = sinl(angle);
    int i;
    #pragma omp parallel for private(i)
    for(i=0; i<ntotal; i++){
        const unsigned int ix = i % n[0];
        const long double x = c[0] + (ix+0.5)*pixsize - n[0]*pixsize/2;
        const long double px = rotc[0] * ca + rotc[1] * sa;
        const long double dx = (x-px)*(x-px);
        if(dx>=1){
            proj[i]=0;
        }else{
            proj[i] = 2*sqrtl(1 - dx);
        }
    }

    const unsigned int nthr = omp_get_max_threads();
    long double * const temp = (long double *) calloc(nthr*ntotal,sizeof(long double));

    #pragma omp parallel
    {
        const unsigned int tidx = ntotal*omp_get_thread_num();
        #pragma omp for schedule(dynamic) private(i)
        for(i=0; i<nspheres; i++){
            const long double s2 = spheres[5*i+3]*spheres[5*i+3];
            const long double py = spheres[5*i+2];
            const long double px = (rotc[0] + spheres[5*i]) * ca + (rotc[1] + spheres[5*i+1]) * sa;
            const long double pyi = (py - c[1])/pixsize + 0.5*(n[1]-1);
            const long double pxi = (px - c[0])/pixsize + 0.5*(n[0]-1);
            const unsigned int sz = spheres[5*i+3]/pixsize + 1;
            int l = pxi-sz;
            int r = pxi+sz;
            int u = pyi-sz;
            int d = pyi+sz;
            if(l<0) l=0;
            if(r>=n[0]) r=n[0]-1;
            if(u<0) u=0;
            if(d>=n[1]) d=n[1]-1;
            for(unsigned int j=u; j<=d;j++){
                const long double y = c[1] + (j+0.5)*pixsize - n[1]*pixsize/2;
                const long double dy = (y-spheres[5*i+2])*(y-spheres[5*i+2]);
                if(dy >= s2) continue;
                for(unsigned int k=l; k<=r; k++){
                    const long double x = c[0] + (k+0.5)*pixsize - n[0]*pixsize/2;
                    const long double dx = (x-px)*(x-px);
                    if(dx+dy<s2){
                        temp[tidx+j*n[0]+k] -= 2*(1-spheres[5*i+4])*sqrtl(s2 - dx - dy);
                    }
                }
            }
        }
    }

    #pragma omp parallel for private(i)
    for(i=0; i<ntotal; i++){
        long double tmpf = 0;
        for(unsigned int j=0;j<nthr;j++){
            tmpf += temp[j*ntotal + i];
        }
        proj[i] += tmpf;
    }
    free(temp);
}

DECLDIR void gen3dproj(const float * const spheres, const unsigned int nspheres, float * const proj, const unsigned int * const n, const float pixsize, const float maxz, const float * const rot, const unsigned int cutout, const float cutoff){
    const unsigned int ntotal = n[0]*n[1];
    int i;

    const unsigned int nthr = omp_get_max_threads();
    long double * const temp = (long double *) calloc(nthr*ntotal,sizeof(long double));
    long double * const mint = (long double *) malloc(nthr*ntotal*sizeof(long double));
    
    #pragma omp parallel for private(i)
    for(i=0;i<nthr*ntotal;i++){
        mint[i] = INFINITY;
    }

    #pragma omp parallel
    {
        const unsigned int tidx = ntotal*omp_get_thread_num();
        #pragma omp for schedule(dynamic) private(i)
        for(i=0; i<nspheres; i++){
            if(spheres[5*i+4]<cutoff) continue;
            const long double s2 = spheres[5*i+3]*spheres[5*i+3];
            const long double px = spheres[5*i];
            const long double pz = spheres[5*i+1];
            const long double py = spheres[5*i+2];
            const long double pyi = py/pixsize + 0.5*(n[1]-1);
            const long double pxi = px/pixsize + 0.5*(n[0]-1);
            const unsigned int sz = spheres[5*i+3]/pixsize + 1;
            int l = pxi-sz;
            int r = pxi+sz;
            int u = pyi-sz;
            int d = pyi+sz;
            if(l<0) l=0;
            if(r>=n[0]) r=n[0]-1;
            if(u<0) u=0;
            if(d>=n[1]) d=n[1]-1;
            for(unsigned int j=u; j<=d;j++){
                const long double y = (j+0.5)*pixsize - n[1]*pixsize/2;
                const long double dy = (y-spheres[5*i+2])*(y-spheres[5*i+2]);
                if(dy >= s2) continue;
                for(unsigned int k=l; k<=r; k++){
                    const long double x = (k+0.5)*pixsize - n[0]*pixsize/2;
                    const long double dx = (x-px)*(x-px);
                    if(dx+dy<s2){
                        long double t = pz - sqrtl(s2 - dx - dy);
                        long double rx = rot[0]*x + rot[1]*t + rot[2]*y;
                        long double ry = rot[3]*x + rot[4]*t + rot[5]*y;
                        long double rz = rot[6]*x + rot[7]*t + rot[8]*y;
                        long double fac = 1;
                        char reject = fabsl(rz)>maxz;
                        if(!reject && cutout>0){
                            reject = rx<0 && ry<0 && rz<0;
                        }
                        if (reject){
                            fac = 0.75;
                            t = pz + sqrtl(s2 - dx - dy);
                            rx = rot[0]*x + rot[1]*t + rot[2]*y;
                            ry = rot[3]*x + rot[4]*t + rot[5]*y;
                            rz = rot[6]*x + rot[7]*t + rot[8]*y;
                            reject = fabsl(rz)>maxz;
                            if(!reject && cutout>0){
                                reject = rx<0 && ry<0 && rz<0;
                            }
                            if (reject){
                                continue;
                            }
                        }
                        if(t < mint[tidx+j*n[0]+k]){
                            mint[tidx+j*n[0]+k] = t;
                            temp[tidx+j*n[0]+k] = fac*sqrtl(s2 - dx - dy)/spheres[5*i+3];
                        }
                    }
                }
            }
        }
    }

    #pragma omp parallel for private(i)
    for(i=0; i<ntotal; i++){
        long double tmpf = temp[i];
        long double cmint = mint[i];
        for(unsigned int j=1;j<nthr;j++){
            if(mint[j*ntotal + i]<cmint){
                cmint = mint[j*ntotal + i];
                tmpf = temp[j*ntotal + i];
            }
            
        }
        proj[i] = tmpf;
    }
    free(temp);
    free(mint);
}

DECLDIR void genconeproj(const float * const spheres, const unsigned int nspheres, float * const proj, const unsigned int * const n, const float pixsize, const float zoff, const float sod, const float sdd){
    const unsigned int ntotal = n[0]*n[1];
    int i;
    #pragma omp parallel for schedule(dynamic) private(i)
    for(i=0; i<ntotal; i++){
        long double tmp=0;
        const unsigned int ix = i % n[0];
        const unsigned int iy = i / n[0];
        const long double x = (ix+0.5)*pixsize - n[0]*pixsize/2;
        const long double y = (iy+0.5)*pixsize - n[1]*pixsize/2;
        
        // https://math.stackexchange.com/questions/2613781/line-cylinder-intersection
        // long double bxl, bxr, byl, byr, bzl, bzr;
        const long double x0 = -sod;
        const long double k = x/sdd;
        const long double l = y/sdd;
        const long double df = 1 - (x0*x0 - 1)*k*k;
        if(df>0){
            const long double t1 = -(sqrtl(df)+x0)/(k*k+1);
            const long double t2 = (sqrtl(df)-x0)/(k*k+1);
            const long double dx = t2 - t1;
            const long double dy = k*dx;
            const long double dz = l*dx;
            tmp = sqrtl(dx*dx+dy*dy+dz*dz);
        }else{
            proj[i]=0;
            continue;
        }

        const long double sz = sqrtl(1+k*k+l*l);
        const long double tx = 1/sz;
        const long double ty = k/sz;
        const long double tz = l/sz;

        // https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
        for(unsigned int j=0; j<nspheres; j++){
            const long double sz = spheres[5*j+2]-zoff;
            const long double sx = spheres[5*j];
            const long double sy = spheres[5*j+1];
            const long double s2 = spheres[5*j+3]*spheres[5*j+3];
            const long double psz = (x0-sx)*tx - sy*ty - sz*tz;
            const long double pdx = x0-sx - psz*tx;
            const long double pdy = -sy - psz*ty;
            const long double pdz = -sz - psz*tz;
            const long double dist = pdx*pdx+pdy*pdy+pdz*pdz;
            if(dist<s2){
                tmp -= 2*(1-spheres[5*j+4])*sqrtl(s2 - dist);
            }
        }
        proj[i] = tmp;
    }
}
    

DECLDIR void average2d(const float * const invol, float * const outvol, const unsigned int nx, const unsigned int ny, const unsigned int ss){
    const unsigned int npix = nx*ny;
    int i;
    #pragma omp parallel for private(i)
    for(i=0; i<npix; i++){
        const unsigned int ix = i % nx;
        const unsigned int iy = i / nx;
        long double tmp = 0;
        for(unsigned int y2=ss*iy;y2<ss*(iy+1);y2++){
            for(unsigned int x2=ss*ix;x2<ss*(ix+1);x2++){
                tmp += invol[y2*ss*nx + x2];
            }
        }
        outvol[i] = tmp/(ss*ss);
    }
}

DECLDIR unsigned int gettouching(const float * const spheres, const unsigned int nspheres, const unsigned int i, const float cutoff, unsigned int * const indices){
    const unsigned int nthr = omp_get_max_threads();
    unsigned int * const js = (unsigned int *) malloc(nthr*nspheres*sizeof(unsigned int));
    unsigned int * const njs = (unsigned int *) malloc(nthr*sizeof(unsigned int));
    for(unsigned int j=0;j<nthr;j++) njs[j]=0;

    const float x = spheres[5*i];
    const float y = spheres[5*i+1];
    const float z = spheres[5*i+2];
    const float s = spheres[5*i+3];

    #pragma omp parallel
    {
        unsigned int tid = omp_get_thread_num();
        unsigned int nj = 0;
        int j;
        #pragma omp for private(j)
        for(j=0; j<nspheres; j++){
            if(i==j) continue;
            const float dst = sqrtf((spheres[5*j]-x)*(spheres[5*j]-x)+(spheres[5*j+1]-y)*(spheres[5*j+1]-y)+(spheres[5*j+2]-z)*(spheres[5*j+2]-z)) - spheres[5*j+3] - s;
            if(dst<cutoff){
                js[tid*nspheres + nj] = j;
                nj++;
            }
        }
        njs[tid] = nj;
    }
    unsigned int nindices = 0;
    for(unsigned int j=0; j<nthr; j++){
        for(unsigned int k=0; k<njs[j]; k++){
            indices[nindices] = js[j*nspheres + k];
            nindices++;
        }
    }
    free(js);
    free(njs);
    return nindices;
}

// Adapted from https://www.johndcook.com/blog/cpp_random_number_generation/
unsigned int poissonsmall(const float lambda){
    float p = 1.0;
    const float L = expf(-lambda);
	unsigned int k = 0;
	do
	{
		k++;
		p *= genRand(&randgen);
	}
	while (p > L);
    return k - 1;
}

unsigned int poissonlarge(const float lambda){
	// "Rejection method PA" from "The Computer Generation of Poisson Random Variables" by A. C. Atkinson
	// Journal of the Royal Statistical Society Series C (Applied Statistics) Vol. 28, No. 1. (1979)
	// The article is on pages 29-35. The algorithm given here is on page 32.

	float c = 0.767 - 3.36/lambda;
	float beta = M_PI/sqrtf(3.0*lambda);
	float alpha = beta*lambda;
	float k = logf(c) - lambda - logf(beta);

	for(;;)
	{
		float u = genRand(&randgen);
		float x = (alpha - logf((1.0 - u)/u))/beta;
		int n = (int) floor(x + 0.5);
		if (n < 0)
			continue;
		float v = genRand(&randgen);
		float y = alpha - beta*x;
        float temp = 1.0 + expf(y);
		float lhs = y + logf(v/(temp*temp));
		float rhs = k + n*logf(lambda) - lgammaf(n+1);
		if (lhs <= rhs)
			return n;
	}
}

float poisson(const float lambda){
    return (lambda < 30.0) ? poissonsmall(lambda) : poissonlarge(lambda);
}




DECLDIR void applypoisson(float * const proj, const unsigned int npix, const float flux, const float factor){
    int i;
    #pragma omp parallel for private(i)
    for(i=0; i<npix; i++){
        float tmp = poisson(flux*expf(-proj[i]*factor));
        if(tmp<=0) tmp = 1;
        proj[i] = -logf(tmp/flux)/factor;
    }
}

/* An implementation of the MT19937 Algorithm for the Mersenne Twister
 * by Evan Sultanik.  Based upon the pseudocode in: M. Matsumoto and
 * T. Nishimura, "Mersenne Twister: A 623-dimensionally
 * equidistributed uniform pseudorandom number generator," ACM
 * Transactions on Modeling and Computer Simulation Vol. 8, No. 1,
 * January pp.3-30 1998.
 *
 * http://www.sultanik.com/Mersenne_twister
 */

#define UPPER_MASK		0x80000000
#define LOWER_MASK		0x7fffffff
#define TEMPERING_MASK_B	0x9d2c5680
#define TEMPERING_MASK_C	0xefc60000

inline void m_seedRand(MTRand* rand, unsigned long seed) {
  /* set initial seeds to mt[STATE_VECTOR_LENGTH] using the generator
   * from Line 25 of Table 1 in: Donald Knuth, "The Art of Computer
   * Programming," Vol. 2 (2nd Ed.) pp.102.
   */
  rand->mt[0] = seed & 0xffffffff;
  for(rand->index=1; rand->index<STATE_VECTOR_LENGTH; rand->index++) {
    rand->mt[rand->index] = (6069 * rand->mt[rand->index-1]) & 0xffffffff;
  }
}

/**
* Creates a new random number generator from a given seed.
*/
MTRand seedRand(unsigned long seed) {
  MTRand rand;
  m_seedRand(&rand, seed);
  return rand;
}

/**
 * Generates a pseudo-randomly generated long.
 */
unsigned long genRandLong(MTRand* rand) {

  unsigned long y;
  static unsigned long mag[2] = {0x0, 0x9908b0df}; /* mag[x] = x * 0x9908b0df for x = 0,1 */
  if(rand->index >= STATE_VECTOR_LENGTH || rand->index < 0) {
    /* generate STATE_VECTOR_LENGTH words at a time */
    int kk;
    if(rand->index >= STATE_VECTOR_LENGTH+1 || rand->index < 0) {
      m_seedRand(rand, 4357);
    }
    for(kk=0; kk<STATE_VECTOR_LENGTH-STATE_VECTOR_M; kk++) {
      y = (rand->mt[kk] & UPPER_MASK) | (rand->mt[kk+1] & LOWER_MASK);
      rand->mt[kk] = rand->mt[kk+STATE_VECTOR_M] ^ (y >> 1) ^ mag[y & 0x1];
    }
    for(; kk<STATE_VECTOR_LENGTH-1; kk++) {
      y = (rand->mt[kk] & UPPER_MASK) | (rand->mt[kk+1] & LOWER_MASK);
      rand->mt[kk] = rand->mt[kk+(STATE_VECTOR_M-STATE_VECTOR_LENGTH)] ^ (y >> 1) ^ mag[y & 0x1];
    }
    y = (rand->mt[STATE_VECTOR_LENGTH-1] & UPPER_MASK) | (rand->mt[0] & LOWER_MASK);
    rand->mt[STATE_VECTOR_LENGTH-1] = rand->mt[STATE_VECTOR_M-1] ^ (y >> 1) ^ mag[y & 0x1];
    rand->index = 0;
  }
  y = rand->mt[rand->index++];
  y ^= (y >> 11);
  y ^= (y << 7) & TEMPERING_MASK_B;
  y ^= (y << 15) & TEMPERING_MASK_C;
  y ^= (y >> 18);
  return y;
}

/**
 * Generates a pseudo-randomly generated double in the range [0..1].
 */
float genRand(MTRand* rand) {
  return((float)genRandLong(rand) / (unsigned long)0xffffffff);
}



/**
 * Skip list implementation
 * Follows http://www.mathcs.emory.edu/~cheung/Courses/323/Syllabus/Map/skip-list-impl.html
 */

struct sknode {
    float val;
    unsigned int idx;
    struct sknode * left;
    struct sknode * right;
    struct sknode * up;
    struct sknode * down;
};


void init_sknode(struct sknode * c){
    c->right = NULL;
    c->left = NULL;
    c->up = NULL;
    c->down = NULL;
}


struct sknode * skhead = NULL;
struct sknode * sktail = NULL;
unsigned int skh = 1;
unsigned int nskpool = 0;
struct sknode * skpool = NULL;

void clear_skiplist(){
    struct sknode * cur;
    struct sknode * next;
    if (skpool) {
        for (unsigned int i=0; i<nskpool; i++){
            cur = skpool[i].up;
            while (cur){
                next = cur->up;
                free(cur);
                cur = next;
            }
        }
        free(skpool);
        skpool = NULL;
    }
    if (skhead) {
        cur = skhead;
        while (cur){
            next = cur->down;
            free(cur);
            cur = next;
        }
        skhead = NULL;
    }
    if (sktail) {
        cur = sktail;
        while (cur){
            next = cur->down;
            free(cur);
            cur = next;
        }
        sktail = NULL;
    }
}

MTRand skrand;

struct sksort {
    unsigned int idx;
    float val;
};

int skcompare(const void *a, const void *b){
    struct sksort *aa = (struct sksort *)a;
    struct sksort *bb = (struct sksort *)b;
    return (aa->val > bb->val) - (aa->val < bb->val);
}

void insert_node_to_skiplist(struct sknode * node, struct sknode * prev){
    while (1){
        struct sknode * right = prev->right;
        prev->right = node;
        node->left = prev;
        node->right = right;
        right->left = node;
        if (genRand(&skrand)>0.5){
            struct sknode * newnode = (struct sknode *) malloc(sizeof(struct sknode));
            init_sknode(newnode);
            newnode->val = node->val;
            newnode->idx = node->idx;
            newnode->down = node;
            node->up = newnode;
            node = newnode;
            while(prev->up == NULL){
                if (prev->left){
                    prev = prev->left;
                }else{
                    struct sknode * newhead = (struct sknode *) malloc(sizeof(struct sknode));
                    init_sknode(newhead);
                    newhead->val = -INFINITY;
                    newhead->down = skhead;
                    skhead->up = newhead;
                    skhead = newhead;
                    
                    struct sknode * newtail = (struct sknode *) malloc(sizeof(struct sknode));
                    init_sknode(newtail);
                    newtail->val = INFINITY;
                    newtail->down = sktail;
                    sktail->up = newtail;
                    sktail = newtail;

                    skhead->right = sktail;
                    sktail->left = skhead;
                    break;
                }
            }
            prev = prev->up;
        }else{
            break;
        }
    }
}

void add_idx_to_skiplist(const unsigned int * const idx, const unsigned int nidx){
    struct sksort * sortlist = (struct sksort *) malloc(nidx * sizeof(struct sksort));
    for (unsigned int i=0; i<nidx; i++){
        sortlist[i].idx = idx[i];
        sortlist[i].val = skpool[idx[i]].val;
    }
    qsort (sortlist, nidx, sizeof(struct sksort), skcompare);
    struct sknode * cur = skhead;
    unsigned char up_possible = 0;
    for (unsigned int i=0; i<nidx; i++){
        float cv = sortlist[i].val;
        while (1){
            while (cur->right->val < cv){
                if (up_possible && cur->up){
                    cur = cur->up;
                }else{
                    cur = cur->right;
                }
            }
            up_possible = 0;
            if (cur->down){
                cur = cur->down;
            }else{
                break;
            }
        }
        insert_node_to_skiplist(skpool + sortlist[i].idx, cur);
        up_possible=1;
    }
    free(sortlist);
}

void remove_idx_from_skiplist(const unsigned int * const idx, const unsigned int nidx){
    struct sknode * cur;
    struct sknode * next;
    for (unsigned int i=0; i<nidx; i++){
        cur = skpool + idx[i];
        cur->left->right = cur->right;
        cur->right->left = cur->left;
        cur = cur->up;
        if (cur) cur->down->up = NULL;
        while(cur){
            cur->left->right = cur->right;
            cur->right->left = cur->left;
            next = cur->up;
            free(cur);
            cur = next;
        }
    }
}

DECLDIR void update_skiplist(const float * const sizes, const unsigned int * const idx, const unsigned int nidx){
    remove_idx_from_skiplist(idx, nidx);
    for (unsigned int i=0; i<nidx; i++){
        skpool[idx[i]].val = sizes[idx[i]];
    }
    add_idx_to_skiplist(idx, nidx);
}

struct sknode * skiplist_iter;
DECLDIR void reset_iter_skiplist(){
    skiplist_iter = skhead;
    while(skiplist_iter->down){
        skiplist_iter = skiplist_iter->down;
    }
    skiplist_iter = skiplist_iter->right;
}

DECLDIR int iter_skiplist(){
    if (skiplist_iter->right){
        skiplist_iter = skiplist_iter->right;
        return skiplist_iter->left->idx;
    }else{
        return -1;
    }
}

DECLDIR void check_skiplist(){
    float cur = -INFINITY;
    reset_iter_skiplist();
    int idx = iter_skiplist();
    while(idx!=-1){
        if (skpool[idx].val<cur){
            printf("SKIPLIST ERROR! %f < %f\n", skpool[idx].val, cur);
        }
        idx = iter_skiplist();
    }
}

DECLDIR void init_skiplist(float * sizes, unsigned int ns){
    clear_skiplist();
    m_seedRand(&skrand, 12345);
    skhead = (struct sknode * ) malloc(sizeof(struct sknode));
    init_sknode(skhead);
    skhead->val = -INFINITY;
    sktail = (struct sknode * ) malloc(sizeof(struct sknode));
    init_sknode(sktail);
    sktail->val = INFINITY;
    skhead->right = sktail;
    sktail->left = skhead;
    skpool = (struct sknode * ) malloc(ns*sizeof(struct sknode));
    nskpool = ns;
    for (unsigned int i=0; i<ns; i++){
        skpool[i].idx = i;
        skpool[i].val = sizes[i];
        init_sknode(skpool + i);
    }

    unsigned int *idx = (unsigned int *) malloc(ns*sizeof(unsigned int));
    for (unsigned int i=0; i<ns ; i++){
        idx[i] = i;
    }

    add_idx_to_skiplist(idx, ns);
    
    free(idx);
}
