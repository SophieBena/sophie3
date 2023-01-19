//
//  BermanStoichiometricPhases.h
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/15/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "BermanProperties.h"

/**
 Contains a list of interface declarations for all Berman-type stoichiometric phases.
 Convention for naming is PhasenameBerman.
 */

@interface AegirineBerman : BermanProperties {
}
@end

@interface AenigmatiteBerman : BermanProperties {
}
@end

@interface AkermaniteBerman : BermanProperties {
}
@end

@interface AndalusiteBerman : BermanProperties {
}
@end

@interface ApatiteBerman : BermanProperties {
}
@end

@interface ChromiteBerman : BermanProperties {
}
@end

@interface CoesiteBerman : BermanProperties {
}
@end

@interface CorundumBerman : BermanProperties {
}
@end

@interface CristobaliteBerman : BermanProperties {
@private
	BermanProperties *alphaCristobalite;
	BermanProperties *betaCristobalite;
	double dTtdp;
}
-(BOOL)isAlphaPhaseAtT:(double)t andP:(double)p;
@end

@interface FayaliteBerman : BermanProperties {
}
@end

@interface ForsteriteBerman : BermanProperties {
}
@end

@interface GehleniteBerman : BermanProperties {
@private
	double d0, d1, d2, d3, d4, d5;
}
@end

@interface HematiteBerman : BermanProperties {
}
@end

@interface IlmeniteBerman : BermanProperties {
}
@end

@interface KalsiliteBerman : BermanProperties {
}
@end

@interface KyaniteBerman : BermanProperties {
}
@end

@interface LimeBerman : BermanProperties {
}
@end

@interface MagnetiteBerman : BermanProperties {
}
@end

@interface MuscoviteBerman : BermanProperties {
}
@end

@interface NephelineBerman : BermanProperties {
}
@end

@interface PericlaseBerman : BermanProperties {
}
@end

@interface PerovskiteBerman : BermanProperties {
}
@end

@interface PhlogopiteBerman : BermanProperties {
}
@end

@interface QuartzBerman : BermanProperties {
@private
	BermanProperties *alphaQuartz;
	BermanProperties *betaQuartz;
	double dTtdp;
}
+(void)enableQuartzCorrectionUsed;
+(void)disableQuartzCorrectionUsed;
-(BOOL)isAlphaPhaseAtT:(double)t andP:(double)p;
@end

@interface RutileBerman : BermanProperties {
}
@end

@interface SanidineBerman : BermanProperties {
@private
	double d0, d1, d2, d3, d4, d5;
}
@end

@interface SillimaniteBerman : BermanProperties {
}
@end

@interface SpheneBerman : BermanProperties {
}
@end

@interface TridymiteBerman : BermanProperties {
@private
	BermanProperties *alphaTridymite;
	BermanProperties *betaTridymite;
}
-(BOOL)isAlphaPhaseAtT:(double)t andP:(double)p;
@end

@interface WhitlockiteBerman : BermanProperties {
}
@end

// Additional Berman Phases

#define NS 2

@interface AlbiteBerman: BermanProperties {
@private
    double tOld, pOld, sOld[NS], invd2gds2[NS][NS];
    double gDis, hDis, sDis, cpDis, dcpdt, vDis, dvdt, dvdp, d2vdt2, d2vdtdp, d2vdp2;
}

#undef NS

@end

@interface High_AlbiteBerman: BermanProperties { }
@end

@interface Low_AlbiteBerman: BermanProperties { }
@end

@interface AlmandineBerman: BermanProperties { }
@end

@interface AnorthiteBerman: BermanProperties { }
@end

@interface AnthophylliteBerman: BermanProperties { }
@end

@interface AntigoriteBerman: BermanProperties { }
@end

@interface BruciteBerman: BermanProperties { }
@end

@interface Ca_Al_PyroxeneBerman: BermanProperties { }
@end

@interface CalciteBerman: BermanProperties { }
@end

@interface ChrysotileBerman: BermanProperties { }
@end

@interface ClinochloreBerman: BermanProperties { }
@end

@interface CordieriteBerman: BermanProperties { }
@end

@interface DiasporeBerman: BermanProperties { }
@end

@interface DiopsideBerman: BermanProperties { }
@end

@interface DolomiteBerman: BermanProperties {
@private
    double d0, d1, d2, d3, d4, d5;
}
@end

@interface ClinoenstatiteBerman: BermanProperties { }
@end

@interface OrthoenstatiteBerman: BermanProperties { }
@end

@interface ProtoenstatiteBerman: BermanProperties { }
@end

@interface FerrosiliteBerman: BermanProperties { }
@end

@interface GrossularBerman: BermanProperties { }
@end

@interface JadeiteBerman: BermanProperties { }
@end

@interface KaoliniteBerman: BermanProperties { }
@end

@interface LawsoniteBerman: BermanProperties { }
@end

@interface MagnesiteBerman: BermanProperties { }
@end

@interface MargariteBerman: BermanProperties { }
@end

@interface MeioniteBerman: BermanProperties { }
@end

@interface MerwiniteBerman: BermanProperties { }
@end

@interface MonticelliteBerman: BermanProperties { }
@end

@interface ParagoniteBerman: BermanProperties { }
@end

@interface Potassium_FeldsparBerman: BermanProperties {
@private
    double d0, d1, d2, d3, d4, d5;
}
@end

@interface MicroclineBerman: BermanProperties { }
@end

@interface PrehniteBerman: BermanProperties { }
@end

@interface PyropeBerman: BermanProperties { }
@end

@interface PyrophylliteBerman: BermanProperties { }
@end

@interface Mg_Al_SpinelBerman: BermanProperties { }
@end

@interface TalcBerman: BermanProperties { }
@end

@interface TremoliteBerman: BermanProperties { }
@end

@interface WollastoniteBerman: BermanProperties { }
@end

@interface PseudowollastoniteBerman: BermanProperties { }
@end

@interface ZoisiteBerman: BermanProperties { }
@end

@interface ClinozoisiteBerman: BermanProperties { }
@end

@interface Oxygen_GasBerman: BermanProperties { }
@end

@interface Sulfur_GasBerman: BermanProperties { }
@end

@interface Hydrogen_GasBerman: BermanProperties { }
@end
