//
//  DoubleVector.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/16/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface DoubleVector : NSObject <NSSecureCoding> {
    NSMutableData *object;
}

@property (readonly) NSUInteger size;

-(id)initWithSize:(NSUInteger)sizeIn;
-(id)initWithSize:(NSUInteger)sizeIn andInitialValue:(double)value;
-(double *)pointerToDouble;

-(double)valueAtIndex:(NSInteger)index;

@end
