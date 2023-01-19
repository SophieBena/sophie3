//
//  Stixrude.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/19/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import "Stixrude.h"

@implementation Stixrude

+(NSString *)convertNumbersToSubscripts:(NSString *)inputString {
    NSString *temp1 = [inputString stringByReplacingOccurrencesOfString:@"2" withString:@"₂"];
    NSString *temp2 = [temp1 stringByReplacingOccurrencesOfString:@"3" withString:@"₃"];
    NSString *temp3 = [temp2 stringByReplacingOccurrencesOfString:@"4" withString:@"₄"];
    NSString *temp4 = [temp3 stringByReplacingOccurrencesOfString:@"5" withString:@"₅"];
    NSString *temp5 = [temp4 stringByReplacingOccurrencesOfString:@"1" withString:@"₁"];
    NSString *temp6 = [temp5 stringByReplacingOccurrencesOfString:@"-" withString:@"_"];
    return temp6;
}

/**
 Provides an array of NSString objects that correspond to the names of oxide composition variables in MELTS.

 @return
 An array of Cocoa strings.
 */
+(NSArray *) oxideNames {
	NSArray *systemOxides = [NSArray arrayWithObjects:@"SiO2", @"Al2O3", @"FeO", @"MgO", @"CaO", @"Na2O", nil];
	NSMutableArray *array = [NSMutableArray arrayWithCapacity:[systemOxides count]];
	for (NSString *entry in systemOxides) {
		[array addObject:[Stixrude convertNumbersToSubscripts:entry]];
	}
	return [NSArray arrayWithArray:array];
}

/**
 Provides an array of NSString objects that correspond to the names of phases known to MELTS.

 @return
 An array of Cocoa strings.
 */
+(NSArray *) phaseNames {
	return [NSArray arrayWithObjects:
			@"Olivine",
			@"Wadsleyite",
			@"Ringwoodite",
			@"Spinel",
			@"Garnet",
			@"Orthopyroxene",
			@"Clinopyroxene",
			@"HP-Clinopyroxene",
			@"Akimotoite",
			@"Perovskite",
			@"Post-Perovskite",
			@"Ferropericlase",
			@"CaPerovskite",
			@"Quartz",
			@"Coesite",
			@"Stishovite",
			@"Seifertite",
			@"CaFerritePhase",
			@"Kyanite",
			@"Nepheline",
			@"Feldspar",
			nil];
}

@end
