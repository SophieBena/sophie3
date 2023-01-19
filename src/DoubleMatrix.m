//
//  DoubleMatrix.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/16/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "DoubleMatrix.h"


@implementation DoubleMatrix

-(id) init {
    _rowSize = 1;
    _colSize = 1;
    return [self initWithRowSize:1 andWithColumnSize:1 andInitialValue:0.0];
}

-(id) initWithRowSize:(NSUInteger)rowSizeIn andWithColumnSize:(NSUInteger)columnSizeIn {
    if ((rowSizeIn < 1) || (columnSizeIn < 1)) return nil;
    return [self initWithRowSize:rowSizeIn andWithColumnSize:columnSizeIn andInitialValue:0.0];
}

-(id)initWithRowSize:(NSUInteger)rowSizeIn andWithColumnSize:(NSUInteger)columnSizeIn andInitialValue:(double)value {
    if ((rowSizeIn < 1) || (columnSizeIn < 1)) return nil;
    if ((self = [super init])) {
        dataObject = [NSMutableData dataWithCapacity:sizeof(double)*columnSizeIn*rowSizeIn];
        double *pointer = [dataObject mutableBytes];
        for (NSUInteger i=0; i<columnSizeIn*rowSizeIn; i++) pointer[i] = value;

        matrixObject = (double **) malloc(rowSizeIn*sizeof(double *));
        for (NSUInteger i=0; i<rowSizeIn; i++) {
            matrixObject[i] = &pointer[i*columnSizeIn];
        }
        _rowSize = rowSizeIn;
        _colSize = columnSizeIn;
    }
    return self;
}

-(double **)pointerToPointerToDouble {
    return matrixObject;
}

-(double)valueAtRowIndex:(NSInteger)rowIndex andColIndex:(NSInteger)colIndex {
    if ((rowIndex < 0) || (rowIndex >= self.rowSize)) return 0.0;
    if ((colIndex < 0) || (colIndex >= self.colSize)) return 0.0;
    return matrixObject[rowIndex][colIndex];
}

-(void)dealloc {
    free(matrixObject);
}

@end
