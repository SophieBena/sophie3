//
//  main.m
//  TestEquilibrateClasses
//
//  Created by Mark Ghiorso on 3/2/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "EquilibrateUsingMELTSv102.h"

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        EquilibrateUsingMELTSv102 *melts = [[EquilibrateUsingMELTSv102 alloc] init];
        NSArray *phases = [melts phaseNames];
        NSLog(@"%lu", phases.count);
        double wt[15];
        wt[0] = 77.8;   // SiO2
        wt[1] =  0.09;  // TiO2
        wt[2] = 12.0;   // Al2O3
        wt[3] =  0.196; // Fe2O3
        wt[4] =  0.0;   // Cr2O3
        wt[5] =  0.474; // FeO
        wt[6] =  0.0;   // MnO
        wt[7] =  0.04;  // MgO
        wt[8] =  0.0;   // NiO
        wt[9] =  0.0;   // CoO
        wt[10] = 0.45;  // CaO
        wt[11] = 3.7;   // Na2O
        wt[12] = 5.36;  // K2O
        wt[13] = 0.0;   // P2O5
        wt[14] = 3.74;  // H2O
        [melts setTemperature:770.0+273.15];
        [melts setPressure:1750.0];
        NSDictionary *results = [melts execute];
        NSLog(@"%@", [results description]);
    }
    return 0;
}
