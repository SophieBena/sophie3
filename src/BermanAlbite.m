//
//  BermanAlbite.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/19/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "BermanAlbite.h"

@implementation BermanAlbite

#define NS 2

-(id)init {
	if ((self = [super initWithH:-3921618.0
							  S:224.412
							 k0:393.64
							 k1:-24.155E2
							 k2:-78.928E5
							 k3:107.064E7
							 v0:10.083
							 v1:-1.945E-6
							 v2:4.861E-12
							 v3:26.307E-6
							 v4:32.407E-10])) {
		[self setPhaseFormula:@"NaAlSi3O8"];
		[self setPhaseName:@"albite"];
		tOld = -9999.0;
		pOld = -9999.0;
		for (NSUInteger i=0; i<NS; i++) sOld[i] = 2.0;
		gDis    = 0.0;
		hDis    = 0.0;
		sDis    = 0.0;
		cpDis   = 0.0;
		dcpdt   = 0.0;
		vDis    = 0.0;
		dvdt    = 0.0;
		dvdp    = 0.0;
		d2vdt2  = 0.0;
		d2vdtdp = 0.0;
		d2vdp2  = 0.0;
	}
	return self;
}

#define SQUARE(x) ((x)*(x))
#define CUBE(x)   ((x)*(x)*(x))

#define MAX_ITER 200

static const double a0   =     5.479;
static const double b	 =  6854.0;
static const double aod0 =    41.620;
static const double bod  = -9301.0;
static const double cod  = 43600.0;
static const double d0   =    -2.171;
static const double d1   =    -3.043;
static const double d2   =    -0.001569;
static const double d3   =     0.000002109;
static const double tc   =  1251.0;
static const double tod  =   824.1;

/* q = sOrd[0], qod = sOrd[1] */

#define S   -(0.5*a0*sOrd[0]*sOrd[0] + 0.5*aod0*sOrd[1]*sOrd[1] \
+ (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[0]*sOrd[1])
#define H   - 0.5*a0*tc*sOrd[0]*sOrd[0] + 0.25 *b*sOrd[0]*sOrd[0]*sOrd[0]*sOrd[0] \
- 0.5*aod0*tod*sOrd[1]*sOrd[1] + 0.25*bod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1] \
+ (cod/6.0)*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1] \
+ (d0 - d2*t*t - 2.0*d3*t*t*t)*sOrd[0]*sOrd[1]
#define V   (H)/335282.925
#define G   (H) - t*(S) + (p-1.0)*(V)

/*----------------------------------------------------------------------------*/

#define DGDS0 (-a0*tc*sOrd[0] + b*sOrd[0]*sOrd[0]*sOrd[0] \
+ (d0-d2*t*t-2.0*d3*t*t*t)*sOrd[1])*(1.0+(p-1.0)/335282.925) \
+ t*(a0*sOrd[0] + (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[1])
#define DGDS1 (-aod0*tod*sOrd[1] + bod*sOrd[1]*sOrd[1]*sOrd[1] \
+ cod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1] \
+ (d0-d2*t*t-2.0*d3*t*t*t)*sOrd[0])*(1.0+(p-1.0)/335282.925) \
+ t*(aod0*sOrd[1] + (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[0])
#define DGDT  -(S) + (-2.0*d2*t-6.0*d3*t*t)*sOrd[0]*sOrd[1]*(p-1.0)/335282.925
#define DGDP  (V)

/*----------------------------------------------------------------------------*/

#define D2GDS0S0 (-a0*tc + 3.0*b*sOrd[0]*sOrd[0])*(1.0+(p-1.0)/335282.925) + t*a0
#define D2GDS0S1 (d0-d2*t*t-2.0*d3*t*t*t)*(1.0+(p-1.0)/335282.925) \
+ t*(d1+2.0*d2*t+3.0*d3*t*t)
#define D2GDS0DT -2.0*(d2*t+3.0*d3*t*t)*sOrd[1]*(p-1.0)/335282.925 \
+ a0*sOrd[0] + (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[1]
#define D2GDS0DP (-a0*tc*sOrd[0] + b*sOrd[0]*sOrd[0]*sOrd[0] \
+ (d0-d2*t*t-2.0*d3*t*t*t)*sOrd[1])/335282.925

#define D2GDS1S1 (-aod0*tod + 3.0*bod*sOrd[1]*sOrd[1] + 5.0*cod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]) \
*(1.0+(p-1.0)/335282.925) + t*aod0
#define D2GDS1DT -2.0*(d2*t+3.0*d3*t*t)*sOrd[0]*(p-1.0)/335282.925 \
+ aod0*sOrd[1] + (d1+2.0*d2*t+3.0*d3*t*t)*sOrd[0]
#define D2GDS1DP (-aod0*tod*sOrd[1] + bod*sOrd[1]*sOrd[1]*sOrd[1] \
+ cod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1] \
+ (d0-d2*t*t-2.0*d3*t*t*t)*sOrd[0])/335282.925

#define D2GDT2   -2.0*(d2+6.0*d3*t)*sOrd[0]*sOrd[1]*(p-1.0)/335282.925 \
+ (2.0*d2+6.0*d3*t)*sOrd[0]*sOrd[1]
#define D2GDTDP  -2.0*(d2*t+3.0*d3*t*t)*sOrd[0]*sOrd[1]/335282.925
#define D2GDP2   0.0

/*----------------------------------------------------------------------------*/

#define D3GDS0S0S0 6.0*b*sOrd[0]*(1.0+(p-1.0)/335282.925)
#define D3GDS0S0S1 0.0
#define D3GDS0S0DT a0
#define D3GDS0S0DP (-a0*tc + 3.0*b*sOrd[0]*sOrd[0])/335282.925
#define D3GDS0S1S1 0.0
#define D3GDS0S1DT (-2.0*d2*t-6.0*d3*t*t)*(1.0+(p-1.0)/335282.925) \
+ d1 + 2.0*d2*t + 3.0*d3*t*t + t*(2.0*d2+6.0*d3*t)
#define D3GDS0S1DP (d0-d2*t*t-2.0*d3*t*t*t)/335282.925
#define D3GDS0DT2  -2.0*(d2+6.0*d3*t)*sOrd[1]*(p-1.0)/335282.925 \
+ (2.0*d2+6.0*d3*t)*sOrd[1]
#define D3GDS0DTDP -2.0*(d2*t+3.0*d3*t*t)*sOrd[1]/335282.925
#define D3GDS0DP2  0.0

#define D3GDS1S1S1 (6.0*bod*sOrd[1] + 20.0*cod*sOrd[1]*sOrd[1]*sOrd[1]) \
*(1.0+(p-1.0)/335282.925)
#define D3GDS1S1DT aod0
#define D3GDS1S1DP (-aod0*tod + 3.0*bod*sOrd[1]*sOrd[1] \
+ 5.0*cod*sOrd[1]*sOrd[1]*sOrd[1]*sOrd[1])/335282.925
#define D3GDS1DT2  -2.0*(d2+6.0*d3*t)*sOrd[0]*(p-1.0)/335282.925 \
+ (2.0*d2+6.0*d3*t)*sOrd[0]
#define D3GDS1DTDP -2.0*(d2*t+3.0*d3*t*t)*sOrd[0]/335282.925
#define D3GDS1DP2  0.0

#define D3GDT3     -12.0*d3*sOrd[0]*sOrd[1]*(p-1.0)/335282.925 + 6.0*d3*sOrd[0]*sOrd[1]
#define D3GDT2DP   -2.0*(d2+6.0*d3*t)*sOrd[0]*sOrd[1]/335282.925
#define D3GDTDP2   0.0
#define D3GDP3     0.0

#define fillD2GDS2 \
d2gds2[0][0] = D2GDS0S0;     d2gds2[0][1] = D2GDS0S1; \
d2gds2[1][0] = d2gds2[0][1]; d2gds2[1][1] = D2GDS1S1;

#define fillD2GDSDT \
d2gdsdt[0] = D2GDS0DT;  d2gdsdt[1] = D2GDS1DT;

#define fillD2GDSDP \
d2gdsdp[0] = D2GDS0DP;  d2gdsdp[1] = D2GDS1DP;

#define fillD3GDS3 \
d3gds3[0][0][0] = D3GDS0S0S0;      d3gds3[0][0][1] = D3GDS0S0S1; \
d3gds3[0][1][0] = d3gds3[0][0][1]; d3gds3[0][1][1] = D3GDS0S1S1; \
d3gds3[1][0][0] = d3gds3[0][0][1]; d3gds3[1][0][1] = d3gds3[0][1][1]; \
d3gds3[1][1][0] = d3gds3[0][1][1]; d3gds3[1][1][1] = D3GDS1S1S1;

#define fillD3GDS2DT \
d3gds2dt[0][0] = D3GDS0S0DT;     d3gds2dt[0][1] = D3GDS0S1DT; \
d3gds2dt[1][0] = d3gds2dt[0][1]; d3gds2dt[1][1] = D3GDS1S1DT;

#define fillD3GDS2DP \
d3gds2dp[0][0] = D3GDS0S0DP;     d3gds2dp[0][1] = D3GDS0S1DP; \
d3gds2dp[1][0] = d3gds2dp[0][1]; d3gds2dp[1][1] = D3GDS1S1DP;

#define fillD3GDSDT2 \
d3gdsdt2[0] = D3GDS0DT2; d3gdsdt2[1] = D3GDS1DT2;

#define fillD3GDSDTDP \
d3gdsdtdp[0] = D3GDS0DTDP; d3gdsdtdp[1] = D3GDS1DTDP;

#define fillD3GDSDP2 \
d3gdsdp2[0] = D3GDS0DP2; d3gdsdp2[1] = D3GDS1DP2;

/**
 Matrix inversion routine converted from Numerical recipies to have zero array indexing
 */

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

-(void)gaussj:(double [NS][NS])a {
	int indxc[NS], indxr[NS], ipiv[NS];
	int i, icol = -1, irow = -1, j, k, l,ll;
	double big, dum, pivinv, temp;

	for (j=0; j<NS; j++) ipiv[j]=0;
	for (i=0; i<NS; i++) {
		big=0.0;
		for (j=0; j<NS; j++)
			if (ipiv[j] != 1)
				for (k=0; k<NS; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					} else if (ipiv[k] > 1) return;
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0; l<NS; l++) SWAP(a[irow][l],a[icol][l])
				}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0) return;
		pivinv = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for (l=0; l<NS; l++) a[icol][l] *= pivinv;
		for (ll=0; ll<NS; ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l=0; l<NS; l++) a[ll][l] -= a[icol][l]*dum;
			}
	}
	for (l=(NS-1); l>=0; l--) {
		if (indxr[l] != indxc[l])
			for (k=0; k<NS; k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
}

/*
 *=============================================================================
 * Local function to compute ordering state and associated derivatives
 */

-(void)order:(double)t
		   p:(double)p
		   s:(double [NS])sOrd        // s[NS]                BINARY MASK: 0000000001
		  dt:(double [NS])dt          // ds[NS]/dt            BINARY MASK: 0000000100
		  dp:(double [NS])dp          // ds[NS]/dp            BINARY MASK: 0000001000
		 dt2:(double [NS])dt2         // d2s[NS]/dt2          BINARY MASK: 0010000000
		 dtp:(double [NS])dtp         // d2s[NS]/dtp          BINARY MASK: 0100000000
		 dp2:(double [NS])dp2         // d2s[NS]/dp2          BINARY MASK: 1000000000
{
	int i, j, k, l, iter=0;
	double d2gdsdt[NS], d2gdsdp[NS], d3gds3[NS][NS][NS], d3gds2dt[NS][NS],
    d3gdsdt2[NS], temp[NS], d3gds2dp[NS][NS], d3gdsdtdp[NS], d3gdsdp2[NS];

	/* look-up or compute the current ordering state */
	if ( (t != tOld) || (p != pOld) ) {
		double dgds[NS], sNew[NS];

		for (i=0; i<NS; i++) sOld[i] = 2.0;

		sNew[0] = 0.6;
		sNew[1] = 0.9;

		while (   ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) || (fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON))
			   && (iter < MAX_ITER)) {
			double sOrd[NS], deltaS[NS];

			for (i=0; i<NS; i++) sOrd[i] = sNew[i];

			dgds[0] = DGDS0;
			dgds[1] = DGDS1;

			invd2gds2[0][0] = D2GDS0S0;
			invd2gds2[0][1] = D2GDS0S1;
			invd2gds2[1][0] = invd2gds2[0][1];
			invd2gds2[1][1] = D2GDS1S1;

			for (i=0; i<NS; i++) sOld[i] = sOrd[i];

			[self gaussj:invd2gds2];

			for (i=0; i<NS; i++) {
				for(j=0; j<NS; j++) sOrd[i] += - invd2gds2[i][j]*dgds[j];
				deltaS[i] = sOrd[i] - sOld[i];
			}

			for (i=0; i<NS; i++) sNew[i] = sOrd[i];
			iter++;
		}
		tOld = t;
		pOld = p;
	}

	/* s */
	for (i=0; i<NS; i++) sOrd[i] = sOld[i];

	/* dsdt */
	fillD2GDSDT

	for (i=0; i<NS; i++) {
		dt[i] = 0.0;
		for (j=0; j<NS; j++) dt[i] += - invd2gds2[i][j]*d2gdsdt[j];
	}

	/* dsdp */
	fillD2GDSDP

	for (i=0; i<NS; i++) {
		dp[i] = 0.0;
		for (j=0; j<NS; j++) dp[i] += - invd2gds2[i][j]*d2gdsdp[j];
	}

	/* d2sdt2 */
	fillD2GDSDT
	fillD3GDS3
	fillD3GDS2DT
	fillD3GDSDT2

	for (i=0; i<NS; i++) {
		for (j=0; j<NS; j++) {
			temp[j] = d3gdsdt2[j];
			for (k=0; k<NS; k++) {
				temp[j] +=  2.0*d3gds2dt[j][k]*dt[k];
				for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dt[k]*dt[l];
			}
		}
		dt2[i] = 0.0;
		for (j=0; j<NS; j++) dt2[i] += - invd2gds2[i][j]*temp[j];
	}

	/* d2sdtdp */
	fillD2GDSDT
	fillD2GDSDP
	fillD3GDS3
	fillD3GDS2DT
	fillD3GDS2DP
	fillD3GDSDTDP

	for (i=0; i<NS; i++) {
		for (j=0; j<NS; j++) {
			temp[j] = d3gdsdtdp[j];
			for (k=0; k<NS; k++) {
				temp[j] += d3gds2dt[j][k]*dp[k] + d3gds2dp[j][k]*dt[k];
				for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dt[k]*dp[l];
			}
		}
		dtp[i] = 0.0;
		for (j=0; j<NS; j++) dtp[i] += - invd2gds2[i][j]*temp[j];
	}

	/* d2sdp2 */
	fillD2GDSDP
	fillD3GDS3
	fillD3GDS2DP
	fillD3GDSDP2

	for (i=0; i<NS; i++) {
		for (j=0; j<NS; j++) {
			temp[j] = d3gdsdp2[j];
			for (k=0; k<NS; k++) {
				temp[j] +=  2.0*d3gds2dp[j][k]*dp[k];
				for (l=0; l<NS; l++) temp[j] += d3gds3[j][k][l]*dp[k]*dp[l];
			}
		}
		dp2[i] = 0.0;
		for (j=0; j<NS; j++) dp2[i] += - invd2gds2[i][j]*temp[j];
	}
}

-(void)albiteSalje:(double)p
				 t:(double)t
{
	double sOrd[NS], dsdt[NS], dsdp[NS], d2sdt2[NS], d2sdtdp[NS], d2sdp2[NS];

	[self order:t p:p s:sOrd dt:dsdt dp:dsdp dt2:d2sdt2 dtp:d2sdtdp dp2:d2sdp2];
	int skip = (fabs(sOrd[0]) < sqrt(DBL_EPSILON)) && (fabs(sOrd[1]) < sqrt(DBL_EPSILON));

	if (skip) {
		gDis    = 0.0;
		hDis    = 0.0;
		sDis    = 0.0;
		cpDis   = 0.0;
		dcpdt   = 0.0;
		vDis    = 0.0;
		dvdt    = 0.0;
		dvdp    = 0.0;
		d2vdt2  = 0.0;
		d2vdtdp = 0.0;
		d2vdp2  = 0.0;
		return;
	} else {
		int i, j, k;
		double d2gds2[NS][NS], d2gdsdt[NS], d2gdsdp[NS], d2gdt2, d2gdtdp, d2gdp2;
		double d3gds3[NS][NS][NS], d3gds2dt[NS][NS], d3gds2dp[NS][NS],
        d3gdsdt2[NS], d3gdsdtdp[NS], d3gdsdp2[NS], d3gdt3, d3gdt2dp, d3gdtdp2, d3gdp3;
		double temp;

		/* ---------- */
		gDis = (G);
		/* ---------- */
		hDis = gDis - t*(DGDT);
		/* ---------- */
		sDis = -(DGDT);
		/* ---------- */
		fillD2GDS2
		fillD2GDSDT
		d2gdt2  = D2GDT2;

		cpDis = d2gdt2;
		for (i=0; i<NS; i++) {
			cpDis += 2.0*d2gdsdt[i]*dsdt[i];
			for (j=0; j<NS; j++) cpDis += d2gds2[i][j]*dsdt[i]*dsdt[j];
		}
		temp = cpDis;
		cpDis *= -t;
		/* ---------- */
		fillD3GDS3
		fillD3GDS2DT
		fillD3GDSDT2
		d3gdt3 = D3GDT3;

		dcpdt = d3gdt3;
		for (i=0; i<NS; i++) {
			dcpdt += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
			for (j=0; j<NS; j++) {
				dcpdt += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j] + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
				for (k=0; k<NS; k++) dcpdt += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
			}
		}
		dcpdt = -t*dcpdt - temp;
		/* ---------- */
		vDis = (V);
		/* ---------- */
		fillD2GDSDP
		d2gdtdp = D2GDTDP;

		dvdt = d2gdtdp;
		for (i=0; i<NS; i++) {
			dvdt += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
			for (j=0; j<NS; j++) dvdt += d2gds2[i][j]*dsdt[i]*dsdp[j];
		}
		/* ---------- */
		d2gdp2 = D2GDP2;

		dvdp = d2gdp2;
		for (i=0; i<NS; i++) {
			dvdp += 2.0*d2gdsdp[i]*dsdp[i];
			for (j=0; j<NS; j++) dvdp += d2gds2[i][j]*dsdp[i]*dsdp[j];
		}
		/* ---------- */
		fillD3GDS2DP
		fillD3GDSDTDP
		d3gdt2dp = D3GDT2DP;

		d2vdt2 = d3gdt2dp;
		for (i=0; i<NS; i++) {
			d2vdt2 += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i] + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
			for (j=0; j<NS; j++) {
				d2vdt2 += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
				for (k=0; k<NS; k++) d2vdt2
					+= d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
			}
		}
		/* ---------- */
		fillD3GDSDP2
		d3gdtdp2 = D3GDTDP2;

		d2vdtdp = d3gdtdp2;
		for (i=0; i<NS; i++) {
			d2vdtdp += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i] + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
			for (j=0; j<NS; j++) {
				d2vdtdp += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j] + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
				+ 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j] + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
				for (k=0; k<NS; k++) d2vdtdp
					+= d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
			}
		}
		/* ---------- */
		d3gdp3 = D3GDP3;

		d2vdp2 = d3gdp3;
		for (i=0; i<NS; i++) {
			d2vdp2 += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
			for (j=0; j<NS; j++) {
				d2vdp2 += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j] + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
				for (k=0; k<NS; k++) d2vdp2 += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
			}
		}
	}

}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getEnthalpyFromT:t andP:p] + hDis;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getEntropyFromT:t andP:p] + sDis;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getHeatCapacityFromT:t andP:p] + cpDis;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getDcpDtFromT:t andP:p] + dcpdt;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getVolumeFromT:t andP:p] + vDis;
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getDvDtFromT:t andP:p] + dvdt;
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getDvDpFromT:t andP:p] + dvdp;
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getD2vDt2FromT:t andP:p] + d2vdt2;
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getD2vDtDpFromT:t andP:p] + d2vdtdp;
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	if ( (t != tOld) || (p != pOld) ) [self albiteSalje:p t:t];
	return [super getD2vDp2FromT:t andP:p] + d2vdp2;
}

@end
