#import <Foundation/Foundation.h>
#import "BermanStoichiometricPhases.h"

int main(int argc, const char *argv[]) {
    @autoreleasepool {
        NSLog(@"Entry ...");
        CoesiteBerman *cb = [[CoesiteBerman alloc] init];
        NSLog(@"Gibbs energy: %lf", [cb getGibbsFreeEnergyFromT:1000.0 andP:1000.0]);
        NSLog(@"Allocated.");
    }
    return 0;
}
