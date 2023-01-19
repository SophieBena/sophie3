//
//  LiquidMeltsGenericEM.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/18/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "BermanProperties.h"

@interface LiquidMeltsGenericEM : BermanProperties {
@private
	double gLiqAtTfusion, hLiqAtTfusion, sLiqAtTfusion;
}

@property (readwrite) double vLiq;
@property (readwrite) double dvdtLiq;
@property (readwrite) double dvdpLiq;
@property (readwrite) double d2vdtdpLiq;
@property (readwrite) double d2vdp2Liq;
@property (readwrite) double tFusion;
@property (readwrite) double sFusion;
@property (readwrite) double cpLiq;

-(id)initWithH:(double)hIn
			 S:(double)sIn
			k0:(double)k0In
			k1:(double)k1In
			k2:(double)k2In
			k3:(double)k3In
			l1:(double)l1In
			l2:(double)l2In
			Tt:(double)TtIn
		deltaH:(double)deltaHIn
		  vLiq:(double)vLiqIn
	   dvdtLiq:(double)dvdtLiqIn
	   dvdpLiq:(double)dvdpLiqIn
	d2vdtdpLiq:(double)d2vdtdpLiqIn
	 d2vdp2Liq:(double)d2vdp2LiqIn
	   tFusion:(double)tFusionIn
	   sFusion:(double)sFusionIn
		 cpLiq:(double)cpLiqIn;

@end
