//
//  main.m
//  TestEquilibrateClasses
//
//  Created by Mark Ghiorso on 3/2/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>
#include "EquilibrateUsingMELTSv102.h"

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        EquilibrateUsingMELTSv102 *melts = [[EquilibrateUsingMELTSv102 alloc] init];
        NSArray *phases = [melts phaseNames];
        NSLog(@"%lu", phases.count);
        double wt[15];
        wt[0] = 77.5;   // SiO2
        wt[1] =  0.08;  // TiO2
        wt[2] = 12.5;   // Al2O3
        wt[3] =  0.207; // Fe2O3
        wt[4] =  0.0;   // Cr2O3
        wt[5] =  0.473; // FeO
        wt[6] =  0.0;   // MnO
        wt[7] =  0.03;  // MgO
        wt[8] =  0.0;   // NiO
        wt[9] =  0.0;   // CoO
        wt[10] = 0.43;  // CaO
        wt[11] = 3.98;  // Na2O
        wt[12] = 4.88;  // K2O
        wt[13] = 0.0;   // P2O5
        wt[14] = 5.5;   // H2O
        [melts setComposition:wt];
        [melts setTemperature:760.0+273.15];
        [melts setPressure:1750.0];
        NSDictionary *results = [melts execute];
        NSLog(@"%@", [results description]);
        NSString *xml = [melts equilibrateResultsAsXML];
        NSLog(@"%@", xml);
    }
    return 0;
}
