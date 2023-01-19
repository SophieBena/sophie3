//
//  IntVector.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/14/12.
//  Copyright (c) 2012 OFM Research Inc. All rights reserved.
//

#import "IntegerVector.h"

@implementation IntegerVector

-(id)init {
    return [self initWithSize:1 andInitialValue:0];
}

-(id)initWithSize:(NSUInteger)sizeIn {
    if (sizeIn < 1) return nil;
    return [self initWithSize:sizeIn andInitialValue:0];
}

-(id)initWithSize:(NSUInteger)sizeIn andInitialValue:(NSInteger)value {
    if (sizeIn < 1) return nil;
    if ((self = [super init])) {
        object = [NSMutableData dataWithCapacity:sizeof(NSInteger)*sizeIn];
        NSInteger *pointer = [object mutableBytes];
        for (NSUInteger i=0; i<sizeIn; i++) pointer[i] = value;
    }
    return self;
}

-(NSInteger *)pointerToInteger {
    return [object mutableBytes];
}

@end
