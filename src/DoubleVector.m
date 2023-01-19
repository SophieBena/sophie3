//
//  DoubleVector.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/16/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "DoubleVector.h"


@implementation DoubleVector

-(id)init {
    _size = 1;
    return [self initWithSize:1 andInitialValue:0.0];
}

-(id)initWithSize:(NSUInteger)sizeIn {
    if (sizeIn < 1) return nil;
    return [self initWithSize:sizeIn andInitialValue:0.0];
}

-(id)initWithSize:(NSUInteger)sizeIn andInitialValue:(double)value {
    if (sizeIn < 1) return nil;
    if ((self = [super init])) {
        _size = sizeIn;
        object = [NSMutableData dataWithCapacity:sizeof(double)*sizeIn];
        double *pointer = [object mutableBytes];
        for (NSUInteger i=0; i<sizeIn; i++) pointer[i] = value;
    }
    return self;
}

-(double *)pointerToDouble {
    return [object mutableBytes];
}

-(double)valueAtIndex:(NSInteger)index {
    if ((index < 0) || (index >= self.size)) return 0.0;
    double *pointer = [object mutableBytes];
    return pointer[index];
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kObject = @"object";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super init])) {
#ifdef __APPLE__
        object = (NSMutableData *) [aDecoder decodeObjectOfClass:[NSMutableData class] forKey:kObject];
#else
        object = (NSMutableData *) [aDecoder decodeObjectForKey:kObject];
#endif
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeObject:object forKey:kObject];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class DoubleVector only supports NSKeyedArchiver coders."];
}

+ (BOOL)supportsSecureCoding {
    return YES;
}

@end
