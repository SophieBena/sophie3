//
//  StixrudeEndmembers.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/9/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface StixrudeEndmembers : NSObject { }

+(NSDictionary *)calculateThermodynamicPropertiesOfEndmemberAtT:(double)t andAtP:(double)p withParameters:(NSArray *)parameters;

@end
