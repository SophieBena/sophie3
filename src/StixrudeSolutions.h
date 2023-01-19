//
//  StixrudeSolutions.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/17/11.
//  Copyright 2011 OFM Research Inc. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "StixrudeSolutionPhase.h"

#define NA 6
#define ND 6
#define NS 16

@interface StixrudeReciprocalSolutionPhase : StixrudeSolutionPhase {
    NSUInteger nd;
    NSArray *dependentEndmembers;
    void (^siteMoleFraction)(double m[NA], double sites[NS]);
    void (^DsiteMoleFraction)(double m[NA], double Dsites[NS][NA]);
    void (^D2siteMoleFraction)(double m[NA], double D2sites[NS][NA][NA]);
    double sites[NS];
    double Dsites[NS][NA];
    double D2sites[NS][NA][NA];
}
@end

@interface FeldsparStixrude : StixrudeSolutionPhase { }
@end
@interface OlivineStixrude : StixrudeSolutionPhase { }
@end
@interface WadsleyiteStixrude : StixrudeSolutionPhase { }
@end
@interface RingwooditeStixrude : StixrudeSolutionPhase { }
@end
@interface PerovskiteStixrude : StixrudeSolutionPhase { }
@end
@interface PostPerovskiteStixrude : StixrudeSolutionPhase { }
@end
@interface OrthopyroxeneStixrude : StixrudeReciprocalSolutionPhase { }
@end
@interface ClinopyroxeneStixrude : StixrudeReciprocalSolutionPhase { }
@end
@interface HPClinopyroxeneStixrude : StixrudeSolutionPhase { }
@end
@interface AkimotoiteStixrude : StixrudeSolutionPhase { }
@end
@interface GarnetStixrude : StixrudeReciprocalSolutionPhase { }
@end
@interface FerropericlaseStixrude : StixrudeSolutionPhase { }
@end
@interface CaFerritePhaseStixrude : StixrudeSolutionPhase { }
@end
@interface SpinelStixrude : StixrudeSolutionPhase { }
@end

#undef NA
#undef ND
#undef NS
