//
//  DoubleMatrix.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/16/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>

@interface DoubleMatrix : NSObject {
    NSMutableData *dataObject;
    double **matrixObject;
}

@property (readonly) NSUInteger rowSize;
@property (readonly) NSUInteger colSize;

-(id)initWithRowSize:(NSUInteger)rowSizeIn andWithColumnSize:(NSUInteger)columnSizeIn;
-(id)initWithRowSize:(NSUInteger)rowSizeIn andWithColumnSize:(NSUInteger)columnSizeIn andInitialValue:(double)value;
-(double **)pointerToPointerToDouble;

-(double)valueAtRowIndex:(NSInteger)rowIndex andColIndex:(NSInteger)colIndex;

@end
