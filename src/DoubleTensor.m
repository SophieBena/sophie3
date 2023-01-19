//
//  DoubleTensor.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/17/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "DoubleTensor.h"


@implementation DoubleTensor

-(id)init {
    _firstSize = 1;
    _secondSize = 1;
    _thirdSize = 1;
    return [self initWithFirstSize:1 andWithSecondSize:1 andWithThirdSize:1 andInitialValue:0.0];
}

-(id)initWithFirstSize:(NSUInteger)firstSizeIn andWithSecondSize:(NSUInteger)secondSizeIn andWithThirdSize:(NSUInteger)thirdSizeIn {
    if ((firstSizeIn < 1) || (secondSizeIn < 1) || (thirdSizeIn < 1)) return nil;
    return [self initWithFirstSize:firstSizeIn andWithSecondSize:secondSizeIn andWithThirdSize:thirdSizeIn andInitialValue:0.0];
}

-(id)initWithFirstSize:(NSUInteger)firstSizeIn andWithSecondSize:(NSUInteger)secondSizeIn andWithThirdSize:(NSUInteger)thirdSizeIn andInitialValue:(double)value {
    if ((firstSizeIn < 1) || (secondSizeIn < 1) || (thirdSizeIn < 1)) return nil;
    if ((self = [super init])) {
        dataObject = [NSMutableData dataWithCapacity:sizeof(double)*firstSizeIn*secondSizeIn*thirdSizeIn];
        double *pointer = [dataObject mutableBytes];
        for (NSUInteger i=0; i<firstSizeIn*secondSizeIn*thirdSizeIn; i++) pointer[i] = value;

        tensorObject = (double ***) malloc(firstSizeIn*sizeof(double **));
        for (NSUInteger i=0; i<firstSizeIn; i++) {
            tensorObject[i] = (double **) malloc(secondSizeIn*sizeof(double *));
            for (NSUInteger j=0; j<secondSizeIn; j++) {
                tensorObject[i][j] = &pointer[i*secondSizeIn*thirdSizeIn + j*thirdSizeIn];
            }
        }
        _firstSize = firstSizeIn;
        _secondSize = secondSizeIn;
        _thirdSize = thirdSizeIn;
    }
    return self;
}

-(double ***)pointerToPointerToPointerToDouble {
    return tensorObject;
}

-(double)valueAtFirstIndex:(NSInteger)firstIndex andSecondIndex:(NSInteger)secondIndex andThirdIndex:(NSInteger)thirdIndex {
    if ((firstIndex < 0) || (firstIndex >= self.firstSize)) return 0.0;
    if ((secondIndex < 0) || (secondIndex >= self.secondSize)) return 0.0;
    if ((thirdIndex < 0) || (thirdIndex >= self.thirdSize)) return 0.0;
    return tensorObject[firstIndex][secondIndex][thirdIndex];
}

-(void)dealloc {
    for (NSUInteger i=0; i<_firstSize; i++) free(tensorObject[i]);
    free(tensorObject);
}

@end
