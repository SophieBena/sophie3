//
//  BermanProperties.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/15/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "BermanProperties.h"
#ifndef __APPLE__
#include <dispatch/dispatch.h>
#endif

@interface BermanPropertiesSettings : NSObject
@property (atomic, assign, getter = isGibbsFreeEnergyReferenceStateUsed) BOOL gibbsFreeEnergyReferenceStateUsed;
@end

@implementation BermanPropertiesSettings
- (instancetype)init {
    self = [super init];
    if (self) {
        _gibbsFreeEnergyReferenceStateUsed = NO;
    }
    return self;
}
@end

@implementation BermanProperties

+(BermanPropertiesSettings *)bermanPropertySettings {
    static BermanPropertiesSettings *bermanPropertySettings = nil;
    static dispatch_once_t onceToken;
    dispatch_once(&onceToken, ^{
        bermanPropertySettings = [[BermanPropertiesSettings alloc] init];
    });
    return bermanPropertySettings;
}

+(void)enableGibbsFreeEnergyReferenceStateUsed {
    @synchronized (self.bermanPropertySettings) {
        if (!self.bermanPropertySettings.isGibbsFreeEnergyReferenceStateUsed) {
            self.bermanPropertySettings.gibbsFreeEnergyReferenceStateUsed = YES;
        }
    }
}

+(void)disableGibbsFreeEnergyReferenceStateUsed {
    @synchronized (self.bermanPropertySettings) {
        if (self.bermanPropertySettings.isGibbsFreeEnergyReferenceStateUsed) {
            self.bermanPropertySettings.gibbsFreeEnergyReferenceStateUsed = NO;
        }
    }
}

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
			v0:(double)v0In
			v1:(double)v1In
			v2:(double)v2In
			v3:(double)v3In
			v4:(double)v4In {
	if ((self = [super init])) {
		h = hIn;
		s = sIn;
		k0 = k0In;
		k1 = k1In;
		k2 = k2In;
		k3 = k3In;
		l1 = l1In;
		l2 = l2In;
		Tt = TtIn;
		deltaH = deltaHIn;
		v0 = v0In;
		v1 = v1In;
		v2 = v2In;
		v3 = v3In;
		v4 = v4In;

		tr = 298.15;
		pr = 1.0;
		trl = 298.15;
	}
	return self;
}

-(id)initWithH:(double)hIn
			 S:(double)sIn
			k0:(double)k0In
			k1:(double)k1In
			k2:(double)k2In
			k3:(double)k3In
			v0:(double)v0In
			v1:(double)v1In
			v2:(double)v2In
			v3:(double)v3In
			v4:(double)v4In {
	return [self initWithH:hIn S:sIn k0:k0In k1:k1In k2:k2In k3:k3In l1:0.0 l2:0.0 Tt:0.0 deltaH:0.0 v0:v0In v1:v1In v2:v2In v3:v3In v4:v4In];
}

-(id)initWithH:(double)hIn
			 S:(double)sIn
			k0:(double)k0In
			k1:(double)k1In
			k2:(double)k2In
			k3:(double)k3In {
	return [self initWithH:hIn S:sIn k0:k0In k1:k1In k2:k2In k3:k3In l1:0.0 l2:0.0 Tt:0.0 deltaH:0.0 v0:0.0 v1:0.0 v2:0.0 v3:0.0 v4:0.0];
}

-(void)setTr:(double)trIn {
	tr = trIn;
}

-(void)setPr:(double)prIn {
	pr = prIn;
}

-(void)setTrl:(double)trlIn {
	trl = trlIn;
}

// Methods inherited from the stoichiometric phase protocol.

-(double)getGibbsFreeEnergyFromT:(double)t andP:(double)p {
	return [self getEnthalpyFromT:t andP:p] - t*[self getEntropyFromT:t andP:p];
}

-(double)getEnthalpyFromT:(double)t andP:(double)p {
	double result = h + k0*(t-tr) + 2.0*k1*(sqrt(t)-sqrt(tr)) - k2*(1.0/t-1.0/tr) - 0.5*k3*(1.0/(t*t)-1.0/(tr*tr));
	if(Tt > 0.0) {
		if (t > Tt) result += deltaH + 0.5*l1*l1*(Tt*Tt-trl*trl) + (2.0/3.0)*l1*l2*(Tt*Tt*Tt-trl*trl*trl) + 0.25*l2*l2*(Tt*Tt*Tt*Tt-trl*trl*trl*trl);
		else        result += 0.5*l1*l1*(t*t-trl*trl) + (2.0/3.0)*l1*l2*(t*t*t-trl*trl*trl) + 0.25*l2*l2*(t*t*t*t-trl*trl*trl*trl);
	}
	result += v0*((v1/2.0-v2)*(p*p-pr*pr)+v2*(p*p*p-pr*pr*pr)/3.0 + (1.0-v1+v2+v3*(t-tr)+v4*(t-tr)*(t-tr))*(p-pr)) - t*v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
    if (BermanProperties.bermanPropertySettings.isGibbsFreeEnergyReferenceStateUsed) result += 298.15*self.entropyFromRobieEtAl1979;
	return result;
}

-(double)getEntropyFromT:(double)t andP:(double)p {
	double result = s + k0*log(t/tr) - 2.0*k1*(1.0/sqrt(t)-1.0/sqrt(tr)) - 0.5*k2*(1.0/(t*t)-1.0/(tr*tr)) - (1.0/3.0)*k3*(1.0/(t*t*t)-1.0/(tr*tr*tr));
	if(Tt > 0.0) {
		if(t > Tt) result += deltaH/Tt + l1*l1*(Tt-trl) + l1*l2*(Tt*Tt-trl*trl) + (1.0/3.0)*l2*l2*(Tt*Tt*Tt-trl*trl*trl);
		else       result += l1*l1*(t-trl) + l1*l2*(t*t-trl*trl) + (1.0/3.0)*l2*l2*(t*t*t-trl*trl*trl);
	}
	result += -v0*(v3 + 2.0*(t-tr)*v4)*(p-pr);
	return result;
}

-(double)getHeatCapacityFromT:(double)t andP:(double)p {
	double result = k0 + k1/sqrt(t) + k2/(t*t) + k3/(t*t*t);
	if (t < Tt) result += t*(l1+l2*t)*(l1+l2*t);
	result += -t*v0*2.0*v4*(p-pr);
	return result;
}

-(double)getDcpDtFromT:(double)t andP:(double)p {
	double result = - 0.5*k1/pow(t, (double) 1.5) - 2.0*k2/(t*t*t) - 3.0*k3/(t*t*t*t);
	if (t < Tt) result += (l1+l2*t)*(l1+l2*t) + t*2.0*(l1+l2*t)*l2;
	result += -v0*2.0*v4*(p-pr);
	return result;
}

-(double)getVolumeFromT:(double)t andP:(double)p {
	return v0*(1.0 + (p-pr)*v1 + (p-pr)*(p-pr)*v2 + (t-tr)*v3 + (t-tr)*(t-tr)*v4);
}

-(double)getDvDtFromT:(double)t andP:(double)p {
	return v0*(v3 + 2.0*(t-tr)*v4);
}

-(double)getDvDpFromT:(double)t andP:(double)p {
	return v0*(v1 + 2.0*(p-pr)*v2);
}

-(double)getD2vDt2FromT:(double)t andP:(double)p {
	return v0*2.0*v4;
}

-(double)getD2vDtDpFromT:(double)t andP:(double)p {
	return 0.0;
}

-(double)getD2vDp2FromT:(double)t andP:(double)p {
	return v0*2.0*v2;
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kh      = @"h";
static NSString *ks      = @"s";
static NSString *kk0     = @"k0";
static NSString *kk1     = @"k1";
static NSString *kk2     = @"k2";
static NSString *kk3     = @"k3";
static NSString *kl1     = @"l1";
static NSString *kl2     = @"l2";
static NSString *kTt     = @"Tt";
static NSString *kdeltaH = @"deltaH";
static NSString *kv0     = @"v0";
static NSString *kv1     = @"v1";
static NSString *kv2     = @"v2";
static NSString *kv3     = @"v3";
static NSString *kv4     = @"v4";
static NSString *ktr     = @"tr";
static NSString *kpr     = @"pr";
static NSString *ktrl    = @"trl";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super initWithCoder:aDecoder])) {
        h      = [aDecoder decodeDoubleForKey:kh];
        s      = [aDecoder decodeDoubleForKey:ks];
        k0     = [aDecoder decodeDoubleForKey:kk0];
        k1     = [aDecoder decodeDoubleForKey:kk1];
        k2     = [aDecoder decodeDoubleForKey:kk2];
        k3     = [aDecoder decodeDoubleForKey:kk3];
        l1     = [aDecoder decodeDoubleForKey:kl1];
        l2     = [aDecoder decodeDoubleForKey:kl2];
        Tt     = [aDecoder decodeDoubleForKey:kTt];
        deltaH = [aDecoder decodeDoubleForKey:kdeltaH];
        v0     = [aDecoder decodeDoubleForKey:kv0];
        v1     = [aDecoder decodeDoubleForKey:kv1];
        v2     = [aDecoder decodeDoubleForKey:kv2];
        v3     = [aDecoder decodeDoubleForKey:kv3];
        v4     = [aDecoder decodeDoubleForKey:kv4];
        tr     = [aDecoder decodeDoubleForKey:ktr];
        pr     = [aDecoder decodeDoubleForKey:kpr];
        trl    = [aDecoder decodeDoubleForKey:ktrl];

    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    [super encodeWithCoder:aCoder];
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeDouble:h      forKey:kh];
        [aCoder encodeDouble:s      forKey:ks];
        [aCoder encodeDouble:k0     forKey:kk0];
        [aCoder encodeDouble:k1     forKey:kk1];
        [aCoder encodeDouble:k2     forKey:kk2];
        [aCoder encodeDouble:k3     forKey:kk3];
        [aCoder encodeDouble:l1     forKey:kl1];
        [aCoder encodeDouble:l2     forKey:kl2];
        [aCoder encodeDouble:Tt     forKey:kTt];
        [aCoder encodeDouble:deltaH forKey:kdeltaH];
        [aCoder encodeDouble:v0     forKey:kv0];
        [aCoder encodeDouble:v1     forKey:kv1];
        [aCoder encodeDouble:v2     forKey:kv2];
        [aCoder encodeDouble:v3     forKey:kv3];
        [aCoder encodeDouble:v4     forKey:kv4];
        [aCoder encodeDouble:tr     forKey:ktr];
        [aCoder encodeDouble:pr     forKey:kpr];
        [aCoder encodeDouble:trl    forKey:ktrl];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class BermanProperties only supports NSKeyedArchiver coders."];
}

@synthesize h, s, k0, k1, k2, k3, l1, l2, Tt, deltaH, v0, v1, v2, v3, v4;

@end
