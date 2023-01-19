//
//  IntVector.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 1/14/12.
//  Copyright (c) 2012 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface IntegerVector : NSObject {
    NSMutableData *object;
}

-(id)initWithSize:(NSUInteger)sizeIn;
-(id)initWithSize:(NSUInteger)sizeIn andInitialValue:(NSInteger)value;
-(NSInteger *)pointerToInteger;

@end
