//
//  DoubleTensor.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/17/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface DoubleTensor : NSObject {
    NSMutableData *dataObject;
    double ***tensorObject;
}

@property (readonly) NSUInteger firstSize;
@property (readonly) NSUInteger secondSize;
@property (readonly) NSUInteger thirdSize;

-(id)initWithFirstSize:(NSUInteger)firstSizeIn andWithSecondSize:(NSUInteger)secondSizeIn andWithThirdSize:(NSUInteger)thirdSizeIn;
-(id)initWithFirstSize:(NSUInteger)firstSizeIn andWithSecondSize:(NSUInteger)secondSizeIn andWithThirdSize:(NSUInteger)thirdSizeIn andInitialValue:(double)value;
-(double ***)pointerToPointerToPointerToDouble;

-(double)valueAtFirstIndex:(NSInteger)firstIndex andSecondIndex:(NSInteger)secondIndex andThirdIndex:(NSInteger)thirdIndex;

@end
