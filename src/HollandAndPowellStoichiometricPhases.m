//
//  HollandAndPowellStoichiometricPhases.m
//  PhaseMELTSobjC
//
//  Created by Mark Ghiorso on 6/8/17.
//  Copyright © 2017 Mark Ghiorso. All rights reserved.
//

#import "HollandAndPowellStoichiometricPhases.h"

@implementation ForsteriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2172200.0
                               S:95.10000000000001
                               a:233.3
                               b:0.001494
                               c:-603800.0
                               d:-1869.6999999999998
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.366
                              a0:6.13e-05
                               K:1250000.0])) {
        [self setPhaseFormula:@"Mg2SiO4"];
        [self setPhaseName:@"Forsterite"];
    }
    return self;
}

@end

@implementation FayaliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1478150.0
                               S:151.0
                               a:201.1
                               b:0.017329999999999998
                               c:-1960600.0
                               d:-900.9000000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.631
                              a0:5.05e-05
                               K:1330000.0])) {
        [self setPhaseFormula:@"Fe2SiO4"];
        [self setPhaseName:@"Fayalite"];
    }
    return self;
}

@end

@implementation TephroiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1732280.0
                               S:155.9
                               a:219.6
                               b:0.0
                               c:-1292700.0
                               d:-1308.3
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.899
                              a0:5.05e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Mn2SiO4"];
        [self setPhaseName:@"Tephroite"];
    }
    return self;
}

@end

@implementation Larnite_BredigiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2307200.0
                               S:127.6
                               a:247.5
                               b:-0.003206
                               c:0.0
                               d:-2051.8999999999996
                             Tc0:1710.0
                            Smax:10.030000000000001
                            Vmax:0.05
                              v0:5.16
                              a0:5.05e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Ca2SiO4"];
        [self setPhaseName:@"Larnite_Bredigite"];
    }
    return self;
}

@end

@implementation MonticelliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2253050.0
                               S:108.10000000000001
                               a:250.7
                               b:-0.010433
                               c:-797200.0
                               d:-1996.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:5.148
                              a0:5.63e-05
                               K:1120000.0])) {
        [self setPhaseFormula:@"CaMgSiO4"];
        [self setPhaseName:@"Monticellite"];
    }
    return self;
}

@end

@implementation ClinohumiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-9607320.0
                               S:445.0
                               a:1070.0
                               b:-0.016533
                               c:-7899600.0
                               d:-7373.9
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:19.785
                              a0:6.1e-05
                               K:1290000.0])) {
        [self setPhaseFormula:@"Mg9Si4O18H2"];
        [self setPhaseName:@"Clinohumite"];
    }
    return self;
}

@end

@implementation PyropeHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6284720.0
                               S:266.29999999999995
                               a:633.5
                               b:0.0
                               c:-5196100.0
                               d:-4315.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:11.318
                              a0:4.36e-05
                               K:1737000.0])) {
        [self setPhaseFormula:@"Mg3Al2Si3O12"];
        [self setPhaseName:@"Pyrope"];
    }
    return self;
}

@end

@implementation AlmandineHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5263520.0
                               S:340.0
                               a:677.3000000000001
                               b:0.0
                               c:-3772700.0
                               d:-5044.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:11.511
                              a0:4.03e-05
                               K:1690000.0])) {
        [self setPhaseFormula:@"Fe3Al2Si3O12"];
        [self setPhaseName:@"Almandine"];
    }
    return self;
}

@end

@implementation SpessartineHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5646400.0
                               S:367.0
                               a:584.6
                               b:-0.001593
                               c:-7516700.0
                               d:-2750.1000000000004
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:11.792
                              a0:4.62e-05
                               K:1790000.0])) {
        [self setPhaseFormula:@"Mn3Al2Si3O12"];
        [self setPhaseName:@"Spessartine"];
    }
    return self;
}

@end

@implementation GrossularHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6644150.0
                               S:255.0
                               a:626.0
                               b:0.0
                               c:-5779200.0
                               d:-4002.9000000000005
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:12.535
                              a0:3.93e-05
                               K:1680000.0])) {
        [self setPhaseFormula:@"Ca3Al2Si3O12"];
        [self setPhaseName:@"Grossular"];
    }
    return self;
}

@end

@implementation AndraditeHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5768130.0
                               S:318.0
                               a:638.5999999999999
                               b:0.0
                               c:-4955100.0
                               d:-3989.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:13.204
                              a0:3.93e-05
                               K:1590000.0])) {
        [self setPhaseFormula:@"Ca3Fe2Si3O12"];
        [self setPhaseName:@"Andradite"];
    }
    return self;
}

@end

@implementation Osumilite_1HollandAndPowell

-(id)init {
    if ((self = [super initWithH:-14968190.0
                               S:701.0
                               a:1625.8
                               b:-0.035547999999999996
                               c:-8063500.0
                               d:-13490.9
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:37.893
                              a0:7.6e-06
                               K:810000.0])) {
        [self setPhaseFormula:@"KMg2Al5Si10O30"];
        [self setPhaseName:@"Osumilite_1"];
    }
    return self;
}

@end

@implementation Osumilite_2HollandAndPowell

-(id)init {
    if ((self = [super initWithH:-14810340.0
                               S:724.0
                               a:1610.6000000000001
                               b:-0.034457
                               c:-8262100.0
                               d:-13128.8
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:38.44
                              a0:7.6e-06
                               K:810000.0])) {
        [self setPhaseFormula:@"KMg3Al3Si11O30"];
        [self setPhaseName:@"Osumilite_2"];
    }
    return self;
}

@end

@implementation Fe_OsumiliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-14248460.0
                               S:762.0
                               a:1656.0
                               b:-0.034163
                               c:-6497700.0
                               d:-14114.3
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:38.32
                              a0:8e-06
                               K:800000.0])) {
        [self setPhaseFormula:@"KFe2Al5Si10O30"];
        [self setPhaseName:@"Fe_Osumilite"];
    }
    return self;
}

@end

@implementation VesuvianiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-42352670.0
                               S:1890.0
                               a:4488.0
                               b:-0.057952
                               c:-22269300.0
                               d:-33478.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:85.2
                              a0:5e-05
                               K:1670000.0])) {
        [self setPhaseFormula:@"Ca19Mg2Al11Si18O78H9"];
        [self setPhaseName:@"Vesuvianite"];
    }
    return self;
}

@end

@implementation AndalusiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2588800.0
                               S:92.7
                               a:277.3
                               b:-0.006588
                               c:-1914100.0
                               d:-2265.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:5.153
                              a0:4.11e-05
                               K:1334000.0])) {
        [self setPhaseFormula:@"Al2SiO5"];
        [self setPhaseName:@"Andalusite"];
    }
    return self;
}

@end

@implementation KyaniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2593110.0
                               S:83.5
                               a:279.4
                               b:-0.007124
                               c:-2055600.0
                               d:-2289.4
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.414
                              a0:4.04e-05
                               K:1590000.0])) {
        [self setPhaseFormula:@"Al2SiO5"];
        [self setPhaseName:@"Kyanite"];
    }
    return self;
}

@end

@implementation SillimaniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2585680.0
                               S:95.5
                               a:280.2
                               b:-0.0069
                               c:-1375700.0
                               d:-2399.4
                             Tc0:2200.0
                            Smax:4.0
                            Vmax:0.035
                              v0:4.986
                              a0:2.21e-05
                               K:1320000.0])) {
        [self setPhaseFormula:@"Al2SiO5"];
        [self setPhaseName:@"Sillimanite"];
    }
    return self;
}

@end

@implementation Hydroxy_TopazHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2904980.0
                               S:100.5
                               a:387.7
                               b:-0.00712
                               c:-857200.0
                               d:-3744.2000000000003
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:5.339
                              a0:4.04e-05
                               K:1315000.0])) {
        [self setPhaseFormula:@"Al2SiO6H2"];
        [self setPhaseName:@"Hydroxy_Topaz"];
    }
    return self;
}

@end

@implementation Mg_StauroliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-25101490.0
                               S:910.0
                               a:2820.5
                               b:-0.059366
                               c:-13774000.0
                               d:-24126.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:44.26
                              a0:1.2e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Mg4Al18Si7.5O48H4"];
        [self setPhaseName:@"Mg_Staurolite"];
    }
    return self;
}

@end

@implementation Fe_StauroliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-23753930.0
                               S:1010.0
                               a:2880.0
                               b:-0.056595
                               c:-10642000.0
                               d:-25373.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:44.88
                              a0:1.2e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Fe4Al18Si7.5O48H4"];
        [self setPhaseName:@"Fe_Staurolite"];
    }
    return self;
}

@end

@implementation Mn_StauroliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-24203880.0
                               S:1024.0
                               a:2873.3
                               b:-0.089064
                               c:-12688000.0
                               d:-24749.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:45.46
                              a0:1.2e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Mn4Al18Si7.5O48H4"];
        [self setPhaseName:@"Mn_Staurolite"];
    }
    return self;
}

@end

@implementation Mg_ChloritoidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3551420.0
                               S:140.0
                               a:464.4
                               b:-0.012654
                               c:-1147200.0
                               d:-4341.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.875
                              a0:5.42e-05
                               K:1465000.0])) {
        [self setPhaseFormula:@"MgAl2SiO7H2"];
        [self setPhaseName:@"Mg_Chloritoid"];
    }
    return self;
}

@end

@implementation Fe_ChloritoidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3215380.0
                               S:155.0
                               a:484.59999999999997
                               b:-0.013808
                               c:-198900.0
                               d:-4762.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.98
                              a0:5.42e-05
                               K:1465000.0])) {
        [self setPhaseFormula:@"FeAl2SiO7H2"];
        [self setPhaseName:@"Fe_Chloritoid"];
    }
    return self;
}

@end

@implementation Mn_ChloritoidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3329280.0
                               S:166.0
                               a:464.4
                               b:-0.012654
                               c:-1147200.0
                               d:-4341.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:7.175
                              a0:5.42e-05
                               K:1465000.0])) {
        [self setPhaseFormula:@"MnAl2SiO7H2"];
        [self setPhaseName:@"Mn_Chloritoid"];
    }
    return self;
}

@end

@implementation MerwiniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4546420.0
                               S:253.1
                               a:417.5
                               b:0.008117000000000001
                               c:-2923000.0
                               d:-2320.3
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:9.847
                              a0:6.15e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Ca3MgSi2O8"];
        [self setPhaseName:@"Merwinite"];
    }
    return self;
}

@end

@implementation SpurriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5849520.0
                               S:330.0
                               a:614.1
                               b:-0.0035080000000000003
                               c:-2493100.0
                               d:-4168.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.697
                              a0:6.5e-05
                               K:950000.0])) {
        [self setPhaseFormula:@"Ca5Si2CO11"];
        [self setPhaseName:@"Spurrite"];
    }
    return self;
}

@end

@implementation ZoisiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6898610.0
                               S:297.0
                               a:595.7
                               b:0.062297000000000005
                               c:-5921300.0
                               d:-3394.7
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:13.575
                              a0:6.7e-05
                               K:1120000.0])) {
        [self setPhaseFormula:@"Ca2Al3Si3O13H"];
        [self setPhaseName:@"Zoisite"];
    }
    return self;
}

@end

@implementation ClinozoisiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6898150.0
                               S:301.0
                               a:567.0
                               b:0.018063
                               c:-7034000.0
                               d:-2603.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:13.63
                              a0:4.6e-05
                               K:1120000.0])) {
        [self setPhaseFormula:@"Ca2Al3Si3O13H"];
        [self setPhaseName:@"Clinozoisite"];
    }
    return self;
}

@end

@implementation Fe_EpidoteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6002180.0
                               S:357.0
                               a:520.1
                               b:0.031499
                               c:-15426000.0
                               d:218.79999999999998
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.19
                              a0:5.05e-05
                               K:1294000.0])) {
        [self setPhaseFormula:@"Ca2AlFe2Si3O13H"];
        [self setPhaseName:@"Fe_Epidote"];
    }
    return self;
}

@end

@implementation EpidoteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6463210.0
                               S:328.0
                               a:544.6
                               b:0.024780999999999997
                               c:-11230000.0
                               d:-1192.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:13.91
                              a0:5.05e-05
                               K:1233000.0])) {
        [self setPhaseFormula:@"Ca2Al2FeSi3O13H"];
        [self setPhaseName:@"Epidote"];
    }
    return self;
}

@end

@implementation LawsoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4869140.0
                               S:230.0
                               a:687.8
                               b:0.001566
                               c:375900.0
                               d:-7179.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.132
                              a0:5.82e-05
                               K:1014000.0])) {
        [self setPhaseFormula:@"CaAl2Si2O10H4"];
        [self setPhaseName:@"Lawsonite"];
    }
    return self;
}

@end

@implementation PumpellyiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-14389500.0
                               S:629.0
                               a:1720.8000000000002
                               b:-0.024928
                               c:-5998700.0
                               d:-14620.300000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:29.55
                              a0:5e-05
                               K:1615000.0])) {
        [self setPhaseFormula:@"Ca4Al5MgSi6O28H7"];
        [self setPhaseName:@"Pumpellyite"];
    }
    return self;
}

@end

@implementation GehleniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3986880.0
                               S:202.0
                               a:405.7
                               b:-0.007099
                               c:-1188300.0
                               d:-3174.4
                             Tc0:700.0
                            Smax:11.0
                            Vmax:0.097
                              v0:9.024
                              a0:4.17e-05
                               K:1080000.0])) {
        [self setPhaseFormula:@"Ca2Al2SiO7"];
        [self setPhaseName:@"Gehlenite"];
    }
    return self;
}

@end

@implementation AkermaniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3866360.0
                               S:212.5
                               a:385.40000000000003
                               b:0.003209
                               c:-247500.0
                               d:-2889.9
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:9.254
                              a0:5.08e-05
                               K:1420000.0])) {
        [self setPhaseFormula:@"Ca2MgSi2O7"];
        [self setPhaseName:@"Akermanite"];
    }
    return self;
}

@end

@implementation RankiniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3944430.0
                               S:210.0
                               a:372.3
                               b:-0.002893
                               c:-2462400.0
                               d:-2181.2999999999997
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:9.651
                              a0:6.5e-05
                               K:950000.0])) {
        [self setPhaseFormula:@"Ca3Si2O7"];
        [self setPhaseName:@"Rankinite"];
    }
    return self;
}

@end

@implementation TilleyiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6368160.0
                               S:390.0
                               a:741.7
                               b:-0.005345
                               c:-1434600.0
                               d:-5878.5
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:17.039
                              a0:6.5e-05
                               K:950000.0])) {
        [self setPhaseFormula:@"Ca5Si2C2O13"];
        [self setPhaseName:@"Tilleyite"];
    }
    return self;
}

@end

@implementation CordieriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-9163370.0
                               S:407.5
                               a:821.3000000000001
                               b:0.043339
                               c:-8211200.000000001
                               d:-5000.0
                             Tc0:1800.0
                            Smax:20.0
                            Vmax:0.2
                              v0:23.322
                              a0:7.6e-06
                               K:810000.0])) {
        [self setPhaseFormula:@"Mg2Al4Si5O18"];
        [self setPhaseName:@"Cordierite"];
    }
    return self;
}

@end

@implementation Hydrous_CordieriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-9446980.0
                               S:487.3
                               a:869.7
                               b:0.051995
                               c:-7723700.0
                               d:-5251.2
                             Tc0:1800.0
                            Smax:20.0
                            Vmax:0.2
                              v0:23.322
                              a0:7.6e-06
                               K:810000.0])) {
        [self setPhaseFormula:@"Mg2Al4Si5O19H2"];
        [self setPhaseName:@"Hydrous_Cordierite"];
    }
    return self;
}

@end

@implementation Fe_CordieriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-8436070.0
                               S:475.0
                               a:851.5
                               b:0.044724
                               c:-6645000.0
                               d:-5623.400000000001
                             Tc0:1800.0
                            Smax:20.0
                            Vmax:0.2
                              v0:23.71
                              a0:7.6e-06
                               K:810000.0])) {
        [self setPhaseFormula:@"Fe2Al4Si5O18"];
        [self setPhaseName:@"Fe_Cordierite"];
    }
    return self;
}

@end

@implementation Mn_CordieriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-8681180.0
                               S:475.0
                               a:847.7
                               b:0.02849
                               c:-7668200.0
                               d:-5311.4
                             Tc0:1800.0
                            Smax:20.0
                            Vmax:0.2
                              v0:24.027
                              a0:7.6e-06
                               K:810000.0])) {
        [self setPhaseFormula:@"Mn2Al4Si5O18"];
        [self setPhaseName:@"Mn_Cordierite"];
    }
    return self;
}

@end

@implementation Phase_AHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-7130410.0
                               S:350.0
                               a:964.0
                               b:-0.011521
                               c:-4517800.0
                               d:-7724.700000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:15.442
                              a0:8.26e-05
                               K:1450000.0])) {
        [self setPhaseFormula:@"Mg7Si2O14H6"];
        [self setPhaseName:@"Phase_A"];
    }
    return self;
}

@end

@implementation SpheneHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2595550.0
                               S:131.20000000000002
                               a:233.7
                               b:0.004043
                               c:-2306500.0
                               d:-1207.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:5.565
                              a0:4.2e-05
                               K:1100000.0])) {
        [self setPhaseFormula:@"CaTiSiO5"];
        [self setPhaseName:@"Sphene"];
    }
    return self;
}

@end

@implementation ZirconHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2032770.0
                               S:84.02999999999999
                               a:237.0
                               b:-0.01788
                               c:-149600.0
                               d:-2267.7999999999997
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.926
                              a0:2.22e-05
                               K:1160000.0])) {
        [self setPhaseFormula:@"ZrSiO4"];
        [self setPhaseName:@"Zircon"];
    }
    return self;
}

@end

@implementation EnstatiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3090470.0
                               S:132.5
                               a:356.2
                               b:-0.00299
                               c:-596900.0
                               d:-3185.2999999999997
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.262
                              a0:5.05e-05
                               K:1070000.0])) {
        [self setPhaseFormula:@"Mg2Si2O6"];
        [self setPhaseName:@"Enstatite"];
    }
    return self;
}

@end

@implementation FerrosiliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2388630.0
                               S:190.6
                               a:398.7
                               b:-0.006579
                               c:1290100.0
                               d:-4058.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.592
                              a0:6.32e-05
                               K:1010000.0])) {
        [self setPhaseFormula:@"Fe2Si2O6"];
        [self setPhaseName:@"Ferrosilite"];
    }
    return self;
}

@end

@implementation Mg_Tschermak_PyroxeneHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3189320.0
                               S:131.0
                               a:371.40000000000003
                               b:-0.004082
                               c:-398400.0
                               d:-3547.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:5.9
                              a0:5.08e-05
                               K:1144000.0])) {
        [self setPhaseFormula:@"MgAl2SiO6"];
        [self setPhaseName:@"Mg_Tschermak_Pyroxene"];
    }
    return self;
}

@end

@implementation DiopsideHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3202760.0
                               S:142.7
                               a:314.5
                               b:4.1e-05
                               c:-2745900.0
                               d:-2020.0999999999997
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.619
                              a0:5.7e-05
                               K:1223000.0])) {
        [self setPhaseFormula:@"CaMgSi2O6"];
        [self setPhaseName:@"Diopside"];
    }
    return self;
}

@end

@implementation HedenbergiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2844160.0
                               S:174.2
                               a:340.2
                               b:0.000812
                               c:-1047800.0
                               d:-2646.7000000000003
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.795
                              a0:5.7e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"CaFeSi2O6"];
        [self setPhaseName:@"Hedenbergite"];
    }
    return self;
}

@end

@implementation JadeiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3027850.0
                               S:133.5
                               a:301.09999999999997
                               b:0.010143
                               c:-2239300.0
                               d:-2055.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.04
                              a0:4.66e-05
                               K:1284000.0])) {
        [self setPhaseFormula:@"NaAlSi2O6"];
        [self setPhaseName:@"Jadeite"];
    }
    return self;
}

@end

@implementation AcmiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2586650.0
                               S:170.6
                               a:307.09999999999997
                               b:0.016758000000000002
                               c:-1685500.0
                               d:-2125.7999999999997
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.459
                              a0:4.66e-05
                               K:1060000.0])) {
        [self setPhaseFormula:@"NaFeSi2O6"];
        [self setPhaseName:@"Acmite"];
    }
    return self;
}

@end

@implementation Ca_Tschermak_PyroxeneHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3307030.0
                               S:138.0
                               a:347.6
                               b:-0.006973999999999999
                               c:-1781600.0
                               d:-2757.5
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.356
                              a0:4.43e-05
                               K:1140000.0])) {
        [self setPhaseFormula:@"CaAl2SiO6"];
        [self setPhaseName:@"Ca_Tschermak_Pyroxene"];
    }
    return self;
}

@end

@implementation RhodoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1321750.0
                               S:100.5
                               a:138.4
                               b:0.004088
                               c:-1936000.0
                               d:-538.9000000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.494
                              a0:5.08e-05
                               K:1250000.0])) {
        [self setPhaseFormula:@"MnSiO3"];
        [self setPhaseName:@"Rhodonite"];
    }
    return self;
}

@end

@implementation PyroxmangiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1322530.0
                               S:99.3
                               a:138.4
                               b:0.004088
                               c:-1936000.0
                               d:-538.9000000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.472
                              a0:5.08e-05
                               K:1250000.0])) {
        [self setPhaseFormula:@"MnSiO3"];
        [self setPhaseName:@"Pyroxmangite"];
    }
    return self;
}

@end

@implementation WollastoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1634060.0
                               S:82.5
                               a:159.3
                               b:0.0
                               c:-967300.0
                               d:-1075.3999999999999
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.993
                              a0:4.6e-05
                               K:795000.0])) {
        [self setPhaseFormula:@"CaSiO3"];
        [self setPhaseName:@"Wollastonite"];
    }
    return self;
}

@end

@implementation PseudowollastoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1627690.0
                               S:88.2
                               a:157.79999999999998
                               b:0.0
                               c:-967300.0
                               d:-1075.3999999999999
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.008
                              a0:5.39e-05
                               K:1050000.0])) {
        [self setPhaseFormula:@"CaSiO3"];
        [self setPhaseName:@"Pseudowollastonite"];
    }
    return self;
}

@end

@implementation TremoliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-12310380.0
                               S:550.0
                               a:1260.2
                               b:0.0038299999999999996
                               c:-11455000.0
                               d:-8237.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:27.27
                              a0:5.34e-05
                               K:762000.0])) {
        [self setPhaseFormula:@"Ca2Mg5Si8O24H2"];
        [self setPhaseName:@"Tremolite"];
    }
    return self;
}

@end

@implementation FerroactinoliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-10511240.0
                               S:705.0
                               a:1290.0
                               b:0.029991
                               c:-8447500.0
                               d:-8947.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:28.28
                              a0:5.34e-05
                               K:760000.0])) {
        [self setPhaseFormula:@"Ca2Fe5Si8O24H2"];
        [self setPhaseName:@"Ferroactinolite"];
    }
    return self;
}

@end

@implementation TschermakiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-12541200.0
                               S:545.0
                               a:1244.8
                               b:0.024348
                               c:-11965000.0
                               d:-8112.099999999999
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:26.8
                              a0:5.34e-05
                               K:760000.0])) {
        [self setPhaseFormula:@"Ca2Mg3Al4Si6O24H2"];
        [self setPhaseName:@"Tschermakite"];
    }
    return self;
}

@end

@implementation PargasiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-12720650.0
                               S:601.0
                               a:1280.2
                               b:0.022997
                               c:-12359500.0
                               d:-8065.799999999999
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:27.19
                              a0:5.34e-05
                               K:912000.0])) {
        [self setPhaseFormula:@"NaCa2Mg4Al3Si6O24H2"];
        [self setPhaseName:@"Pargasite"];
    }
    return self;
}

@end

@implementation GlaucophaneHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-11969470.0
                               S:535.0
                               a:1717.5
                               b:-0.12107000000000001
                               c:7075000.0
                               d:-19272.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:26.05
                              a0:5.3e-05
                               K:883000.0])) {
        [self setPhaseFormula:@"Na2Mg3Al2Si8O24H2"];
        [self setPhaseName:@"Glaucophane"];
    }
    return self;
}

@end

@implementation Fe_GlaucophaneHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-10889740.0
                               S:624.0
                               a:1762.8999999999999
                               b:-0.118992
                               c:9423700.0
                               d:-20207.100000000002
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:26.59
                              a0:5.3e-05
                               K:890000.0])) {
        [self setPhaseFormula:@"Na2Al2Fe3Si8O24H2"];
        [self setPhaseName:@"Fe_Glaucophane"];
    }
    return self;
}

@end

@implementation RiebeckiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-10035240.0
                               S:691.0
                               a:1746.8999999999999
                               b:-0.113572
                               c:9370300.0
                               d:-19468.699999999997
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:27.49
                              a0:5.3e-05
                               K:890000.0])) {
        [self setPhaseFormula:@"Na2Fe5Si8O24H2"];
        [self setPhaseName:@"Riebeckite"];
    }
    return self;
}

@end

@implementation AnthophylliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-12069200.0
                               S:536.0
                               a:1277.3000000000002
                               b:0.025825
                               c:-9704600.0
                               d:-9074.7
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:26.54
                              a0:5e-05
                               K:700000.0])) {
        [self setPhaseFormula:@"Mg7Si8O24H2"];
        [self setPhaseName:@"Anthophyllite"];
    }
    return self;
}

@end

@implementation Fe_AnthophylliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-9627230.0
                               S:725.0
                               a:1383.1
                               b:0.030669000000000002
                               c:-4224700.0
                               d:-11257.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:27.87
                              a0:5e-05
                               K:700000.0])) {
        [self setPhaseFormula:@"Fe7Si8O24H2"];
        [self setPhaseName:@"Fe_Anthophyllite"];
    }
    return self;
}

@end

@implementation CummingtoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-12057480.0
                               S:540.0
                               a:1277.3000000000002
                               b:0.025825
                               c:-9704600.0
                               d:-9074.7
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:26.33
                              a0:5e-05
                               K:700000.0])) {
        [self setPhaseFormula:@"Mg7Si8O24H2"];
        [self setPhaseName:@"Cummingtonite"];
    }
    return self;
}

@end

@implementation GruneriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-9613200.0
                               S:730.0
                               a:1383.1
                               b:0.030669000000000002
                               c:-4224700.0
                               d:-11257.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:27.84
                              a0:5e-05
                               K:648000.0])) {
        [self setPhaseFormula:@"Fe7Si8O24H2"];
        [self setPhaseName:@"Grunerite"];
    }
    return self;
}

@end

@implementation GedriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-12319430.0
                               S:515.0
                               a:1307.7
                               b:0.023642000000000003
                               c:-9307400.0
                               d:-9799.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:25.8
                              a0:4.8e-05
                               K:770000.0])) {
        [self setPhaseFormula:@"Mg5Al4Si6O24H2"];
        [self setPhaseName:@"Gedrite"];
    }
    return self;
}

@end

@implementation Sapphirine_442HollandAndPowell

-(id)init {
    if ((self = [super initWithH:-11019620.0
                               S:443.0
                               a:1160.3000000000002
                               b:-0.024324
                               c:-7706600.0
                               d:-8974.199999999999
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:19.905
                              a0:4.9e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Mg4Al8Si2O20"];
        [self setPhaseName:@"Sapphirine_442"];
    }
    return self;
}

@end

@implementation Sapphirine_793HollandAndPowell

-(id)init {
    if ((self = [super initWithH:-11065630.0
                               S:448.0
                               a:1167.8999999999999
                               b:-0.02487
                               c:-7607300.0
                               d:-9155.300000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:19.775
                              a0:4.9e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Mg3.5Al9Si1.5O20"];
        [self setPhaseName:@"Sapphirine_793"];
    }
    return self;
}

@end

@implementation Fe_SapphirineHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-9835510.0
                               S:551.0
                               a:1257.8
                               b:-0.022171000000000003
                               c:-1664000.0
                               d:-11348.4
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:20.391
                              a0:4.9e-05
                               K:1200000.0])) {
        [self setPhaseFormula:@"Fe3.5Al9Si1.5O20"];
        [self setPhaseName:@"Fe_Sapphirine"];
    }
    return self;
}

@end

@implementation Mg_CarpholiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4781240.0
                               S:210.0
                               a:667.8
                               b:-0.012558999999999999
                               c:-1167100.0
                               d:-6440.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.59
                              a0:5e-05
                               K:525000.0])) {
        [self setPhaseFormula:@"MgAl2Si2O10H4"];
        [self setPhaseName:@"Mg_Carpholite"];
    }
    return self;
}

@end

@implementation Fe_CarpholiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4413200.0
                               S:255.0
                               a:674.8
                               b:-0.010091999999999999
                               c:-715800.0
                               d:-6554.5
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.69
                              a0:5e-05
                               K:525000.0])) {
        [self setPhaseFormula:@"FeAl2Si2O10H4"];
        [self setPhaseName:@"Fe_Carpholite"];
    }
    return self;
}

@end

@implementation DeeriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-18347650.0
                               S:1650.0
                               a:3164.4
                               b:-0.027883
                               c:-5039100.0
                               d:-26721.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:55.74
                              a0:5e-05
                               K:630000.0])) {
        [self setPhaseFormula:@"Fe18Si12O50H10"];
        [self setPhaseName:@"Deerite"];
    }
    return self;
}

@end

@implementation MuscoviteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5984180.0
                               S:292.0
                               a:756.4
                               b:-0.01984
                               c:-2170000.0
                               d:-6979.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.083
                              a0:5.96e-05
                               K:490000.0])) {
        [self setPhaseFormula:@"KAl3Si3O12H2"];
        [self setPhaseName:@"Muscovite"];
    }
    return self;
}

@end

@implementation CeladoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5842420.0
                               S:290.0
                               a:741.1999999999999
                               b:-0.018748
                               c:-2368800.0
                               d:-6616.900000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:13.957
                              a0:5.96e-05
                               K:700000.0])) {
        [self setPhaseFormula:@"KMgAlSi4O12H2"];
        [self setPhaseName:@"Celadonite"];
    }
    return self;
}

@end

@implementation Fe_CeladoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5477430.0
                               S:329.0
                               a:756.3
                               b:-0.019147
                               c:-1586100.0
                               d:-6928.7
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.07
                              a0:5.96e-05
                               K:700000.0])) {
        [self setPhaseFormula:@"KFeAlSi4O12H2"];
        [self setPhaseName:@"Fe_Celadonite"];
    }
    return self;
}

@end

@implementation ParagoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5946340.0
                               S:276.0
                               a:803.0
                               b:-0.03158
                               c:217000.0
                               d:-8151.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:13.211
                              a0:7.74e-05
                               K:550000.0])) {
        [self setPhaseFormula:@"NaAl3Si3O12H2"];
        [self setPhaseName:@"Paragonite"];
    }
    return self;
}

@end

@implementation MargariteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6241230.0
                               S:267.0
                               a:744.4
                               b:-0.0168
                               c:-2074400.0
                               d:-6783.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:12.964
                              a0:4.87e-05
                               K:1300000.0])) {
        [self setPhaseFormula:@"CaAl4Si2O12H2"];
        [self setPhaseName:@"Margarite"];
    }
    return self;
}

@end

@implementation PhlogopiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6219160.0
                               S:328.0
                               a:770.3
                               b:-0.036939
                               c:-2328900.0
                               d:-6531.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.964
                              a0:5.79e-05
                               K:513000.0])) {
        [self setPhaseFormula:@"KMg3AlSi3O12H2"];
        [self setPhaseName:@"Phlogopite"];
    }
    return self;
}

@end

@implementation AnniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5151670.0
                               S:418.0
                               a:815.6999999999999
                               b:-0.034860999999999996
                               c:19800.0
                               d:-7466.700000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:15.432
                              a0:5.79e-05
                               K:513000.0])) {
        [self setPhaseFormula:@"KFe3AlSi3O12H2"];
        [self setPhaseName:@"Annite"];
    }
    return self;
}

@end

@implementation Mn_BiotiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5462840.0
                               S:433.0
                               a:809.9
                               b:-0.059212999999999995
                               c:-1514400.0
                               d:-6998.700000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:15.264
                              a0:5.79e-05
                               K:513000.0])) {
        [self setPhaseFormula:@"KMn3AlSi3O12H2"];
        [self setPhaseName:@"Mn_Biotite"];
    }
    return self;
}

@end

@implementation EastoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6338170.0
                               S:318.0
                               a:785.5
                               b:-0.038031
                               c:-2130300.0
                               d:-6893.7
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.738
                              a0:5.79e-05
                               K:513000.0])) {
        [self setPhaseFormula:@"KMg2Al3Si2O12H2"];
        [self setPhaseName:@"Eastonite"];
    }
    return self;
}

@end

@implementation Na_PhlogopiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6172820.0
                               S:318.0
                               a:773.5
                               b:-0.040228999999999994
                               c:-2597900.0
                               d:-6512.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.45
                              a0:5.79e-05
                               K:513000.0])) {
        [self setPhaseFormula:@"NaMg3AlSi3O12H2"];
        [self setPhaseName:@"Na_Phlogopite"];
    }
    return self;
}

@end

@implementation ClinochloreHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-8912410.0
                               S:430.5
                               a:1161.8
                               b:0.010133
                               c:-7657300.0
                               d:-9690.9
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:21.09
                              a0:3.98e-05
                               K:870000.0])) {
        [self setPhaseFormula:@"Mg5Al2Si3O18H8"];
        [self setPhaseName:@"Clinochlore"];
    }
    return self;
}

@end

@implementation AmesiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-9034440.0
                               S:410.0
                               a:1177.0
                               b:0.009041
                               c:-7458700.0
                               d:-10053.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:20.52
                              a0:3.98e-05
                               K:870000.0])) {
        [self setPhaseFormula:@"Mg4Al4Si2O18H8"];
        [self setPhaseName:@"Amesite"];
    }
    return self;
}

@end

@implementation Al_Free_ChloriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-8727380.0
                               S:428.0
                               a:1146.6000000000001
                               b:0.011225
                               c:-7855900.0
                               d:-9328.8
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:21.66
                              a0:3.98e-05
                               K:870000.0])) {
        [self setPhaseFormula:@"Mg6Si4O18H8"];
        [self setPhaseName:@"Al_Free_Chlorite"];
    }
    return self;
}

@end

@implementation DaphniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-7134850.0
                               S:565.0
                               a:1237.4
                               b:0.013594
                               c:-3743000.0
                               d:-11250.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:21.34
                              a0:3.98e-05
                               K:870000.0])) {
        [self setPhaseFormula:@"Fe5Al2Si3O18H8"];
        [self setPhaseName:@"Daphnite"];
    }
    return self;
}

@end

@implementation Mn_ChloriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-7666530.0
                               S:595.0
                               a:1227.8
                               b:-0.02699
                               c:-6299800.0
                               d:-10469.4
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:22.59
                              a0:3.98e-05
                               K:870000.0])) {
        [self setPhaseFormula:@"Mn5Al2Si3O18H8"];
        [self setPhaseName:@"Mn_Chlorite"];
    }
    return self;
}

@end

@implementation SudoiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-8626160.0
                               S:404.0
                               a:1436.1
                               b:-0.048748999999999994
                               c:-2748500.0
                               d:-13764.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:20.3
                              a0:3.98e-05
                               K:870000.0])) {
        [self setPhaseFormula:@"Mg2Al4Si3O18H8"];
        [self setPhaseName:@"Sudoite"];
    }
    return self;
}

@end

@implementation Fe_SudoiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-7903900.0
                               S:462.0
                               a:1466.3
                               b:-0.047365
                               c:-1182800.0
                               d:-14388.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:20.4
                              a0:3.98e-05
                               K:870000.0])) {
        [self setPhaseFormula:@"Fe2Al4Si3O18H8"];
        [self setPhaseName:@"Fe_Sudoite"];
    }
    return self;
}

@end

@implementation PyrophylliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5640640.0
                               S:239.4
                               a:784.5
                               b:-0.042948
                               c:1251000.0
                               d:-8495.900000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:12.81
                              a0:7.5e-05
                               K:525000.0])) {
        [self setPhaseFormula:@"Al2Si4O12H2"];
        [self setPhaseName:@"Pyrophyllite"];
    }
    return self;
}

@end

@implementation TalcHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5897100.0
                               S:260.0
                               a:622.1999999999999
                               b:0.0
                               c:-6385500.0
                               d:-3916.3
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:13.625
                              a0:3.7e-05
                               K:480000.0])) {
        [self setPhaseFormula:@"Mg3Si4O12H2"];
        [self setPhaseName:@"Talc"];
    }
    return self;
}

@end

@implementation Fe_TalcHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4802920.0
                               S:352.0
                               a:579.7
                               b:0.039494
                               c:-6459300.0
                               d:-3088.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.225
                              a0:3.7e-05
                               K:480000.0])) {
        [self setPhaseFormula:@"Fe3Si4O12H2"];
        [self setPhaseName:@"Fe_Talc"];
    }
    return self;
}

@end

@implementation Tschermak_TalcHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-5987990.0
                               S:259.0
                               a:549.5
                               b:0.036324
                               c:-8606600.0
                               d:-2515.2999999999997
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:13.51
                              a0:3.7e-05
                               K:480000.0])) {
        [self setPhaseFormula:@"Mg2Al2Si3O12H2"];
        [self setPhaseName:@"Tschermak_Talc"];
    }
    return self;
}

@end

@implementation KaoliniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4122180.0000000005
                               S:203.7
                               a:436.7
                               b:-0.034295000000000006
                               c:-4055900.0
                               d:-2699.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:9.934
                              a0:5.1e-05
                               K:645000.0])) {
        [self setPhaseFormula:@"Al2Si2O9H4"];
        [self setPhaseName:@"Kaolinite"];
    }
    return self;
}

@end

@implementation PrehniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6203180.0
                               S:292.8
                               a:724.9
                               b:-0.013865
                               c:-2059000.0
                               d:-6323.9
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:14.026
                              a0:5.1e-05
                               K:835000.0])) {
        [self setPhaseFormula:@"Ca2Al2Si3O12H2"];
        [self setPhaseName:@"Prehnite"];
    }
    return self;
}

@end

@implementation ChrysotileHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4359030.0
                               S:221.29999999999998
                               a:624.7
                               b:-0.02077
                               c:-1721800.0
                               d:-5619.4
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.746
                              a0:4.7e-05
                               K:525000.0])) {
        [self setPhaseFormula:@"Mg3Si2O9H4"];
        [self setPhaseName:@"Chrysotile"];
    }
    return self;
}

@end

@implementation AntigoriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-71424310.0
                               S:3591.0
                               a:9621.0
                               b:-0.091183
                               c:-35941600.0
                               d:-83034.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:175.48
                              a0:4.7e-05
                               K:525000.0])) {
        [self setPhaseFormula:@"Mg48Si34O147H62"];
        [self setPhaseName:@"Antigorite"];
    }
    return self;
}

@end

@implementation AlbiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3934560.0
                               S:210.10000000000002
                               a:452.0
                               b:-0.013364
                               c:-1275900.0
                               d:-3953.6
                             Tc0:950.0
                            Smax:16.0
                            Vmax:0.124
                              v0:10.006
                              a0:4.56e-05
                               K:593000.0])) {
        [self setPhaseFormula:@"NaAlSi3O8"];
        [self setPhaseName:@"Albite"];
    }
    return self;
}

@end

@implementation High_AlbiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3924800.0
                               S:223.39999999999998
                               a:452.0
                               b:-0.013364
                               c:-1275900.0
                               d:-3953.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.109
                              a0:4.56e-05
                               K:593000.0])) {
        [self setPhaseFormula:@"NaAlSi3O8"];
        [self setPhaseName:@"High_Albite"];
    }
    return self;
}

@end

@implementation MicroclineHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3975110.0
                               S:216.0
                               a:448.79999999999995
                               b:-0.010074999999999999
                               c:-1007300.0
                               d:-3973.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.892
                              a0:3.35e-05
                               K:574000.0])) {
        [self setPhaseFormula:@"KAlSi3O8"];
        [self setPhaseName:@"Microcline"];
    }
    return self;
}

@end

@implementation SanidineHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3964960.0
                               S:230.0
                               a:448.79999999999995
                               b:-0.010074999999999999
                               c:-1007300.0
                               d:-3973.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.9
                              a0:3.35e-05
                               K:574000.0])) {
        [self setPhaseFormula:@"KAlSi3O8"];
        [self setPhaseName:@"Sanidine"];
    }
    return self;
}

@end

@implementation AnorthiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4233480.0
                               S:200.0
                               a:371.59999999999997
                               b:0.012615000000000001
                               c:-4110200.0
                               d:-2038.4000000000003
                             Tc0:2300.0
                            Smax:11.0
                            Vmax:0.05
                              v0:10.079
                              a0:2.38e-05
                               K:919000.0])) {
        [self setPhaseFormula:@"CaAl2Si2O8"];
        [self setPhaseName:@"Anorthite"];
    }
    return self;
}

@end

@implementation QuartzHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-910830.0
                               S:41.5
                               a:110.7
                               b:-0.005189
                               c:0.0
                               d:-1128.3000000000002
                             Tc0:847.0
                            Smax:4.95
                            Vmax:0.1188
                              v0:2.2688
                              a0:6.5e-06
                               K:750000.0])) {
        [self setPhaseFormula:@"SiO2"];
        [self setPhaseName:@"Quartz"];
    }
    return self;
}

@end

@implementation TridymiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-906690.0
                               S:46.1
                               a:97.9
                               b:-0.00335
                               c:-636200.0
                               d:-774.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.7
                              a0:5e-06
                               K:750000.0])) {
        [self setPhaseFormula:@"SiO2"];
        [self setPhaseName:@"Tridymite"];
    }
    return self;
}

@end

@implementation CristobaliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-905990.0
                               S:46.5
                               a:97.9
                               b:-0.00335
                               c:-636200.0
                               d:-774.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.61
                              a0:8.1e-06
                               K:600000.0])) {
        [self setPhaseFormula:@"SiO2"];
        [self setPhaseName:@"Cristobalite"];
    }
    return self;
}

@end

@implementation CoesiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-905470.0
                               S:40.800000000000004
                               a:96.5
                               b:-0.000577
                               c:-444800.0
                               d:-798.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.064
                              a0:1.8e-05
                               K:1000000.0])) {
        [self setPhaseFormula:@"SiO2"];
        [self setPhaseName:@"Coesite"];
    }
    return self;
}

@end

@implementation StishoviteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-875630.0
                               S:24.5
                               a:68.1
                               b:0.00601
                               c:-1978200.0
                               d:-82.10000000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:1.4014
                              a0:2.5e-05
                               K:3160000.0])) {
        [self setPhaseFormula:@"SiO2"];
        [self setPhaseName:@"Stishovite"];
    }
    return self;
}

@end

@implementation NephelineHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2095080.0
                               S:124.39999999999999
                               a:272.7
                               b:-0.012398000000000001
                               c:0.0
                               d:-2763.1
                             Tc0:467.0
                            Smax:10.0
                            Vmax:0.08
                              v0:5.419
                              a0:8.1e-05
                               K:600000.0])) {
        [self setPhaseFormula:@"NaAlSiO4"];
        [self setPhaseName:@"Nepheline"];
    }
    return self;
}

@end

@implementation KalsiliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2121920.0
                               S:134.0
                               a:242.0
                               b:-0.004482
                               c:-895800.0
                               d:-1935.8
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.04
                              a0:5.76e-05
                               K:590000.0])) {
        [self setPhaseFormula:@"KAlSiO4"];
        [self setPhaseName:@"Kalsilite"];
    }
    return self;
}

@end

@implementation LeuciteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3029160.0
                               S:200.0
                               a:369.8
                               b:-0.016332
                               c:684700.0
                               d:-3683.1
                             Tc0:938.0
                            Smax:18.0
                            Vmax:0.482
                              v0:8.828
                              a0:3.67e-05
                               K:630000.0])) {
        [self setPhaseFormula:@"KAlSi2O6"];
        [self setPhaseName:@"Leucite"];
    }
    return self;
}

@end

@implementation MeioniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-13843570.0
                               S:752.0
                               a:1359.0
                               b:0.036442
                               c:-8594700.0
                               d:-9598.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:33.985
                              a0:3.16e-05
                               K:870000.0])) {
        [self setPhaseFormula:@"Ca4Al6Si6O27C"];
        [self setPhaseName:@"Meionite"];
    }
    return self;
}

@end

@implementation WairakiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-6666420.0
                               S:375.0
                               a:838.3000000000001
                               b:-0.02146
                               c:-2272000.0
                               d:-7292.3
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:19.04
                              a0:2.38e-05
                               K:1000000.0])) {
        [self setPhaseFormula:@"CaAl2Si4O14H4"];
        [self setPhaseName:@"Wairakite"];
    }
    return self;
}

@end

@implementation LaumontiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-7268470.0
                               S:457.0
                               a:1013.4000000000001
                               b:-0.021413
                               c:-2235800.0
                               d:-8806.699999999999
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:20.37
                              a0:2.38e-05
                               K:1000000.0])) {
        [self setPhaseFormula:@"CaAl2Si4O16H8"];
        [self setPhaseName:@"Laumontite"];
    }
    return self;
}

@end

@implementation HeulanditeHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-10599140.0
                               S:669.0
                               a:1504.8
                               b:-0.033224
                               c:-2959300.0
                               d:-13297.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:31.8
                              a0:2.38e-05
                               K:1000000.0])) {
        [self setPhaseFormula:@"CaAl2Si7O24H12"];
        [self setPhaseName:@"Heulandite"];
    }
    return self;
}

@end

@implementation StilbiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-10898830.0
                               S:710.0
                               a:1588.4
                               b:-0.032042999999999995
                               c:-3071600.0
                               d:-13966.900000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:32.87
                              a0:2.38e-05
                               K:1000000.0])) {
        [self setPhaseFormula:@"CaAl2Si7O25H14"];
        [self setPhaseName:@"Stilbite"];
    }
    return self;
}

@end

@implementation AnalciteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3309900.0
                               S:232.0
                               a:643.5
                               b:-0.016067
                               c:9302300.0
                               d:-9179.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:9.74
                              a0:5e-05
                               K:400000.0])) {
        [self setPhaseFormula:@"NaAlSi2O7H2"];
        [self setPhaseName:@"Analcite"];
    }
    return self;
}

@end

@implementation LimeHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-634950.0
                               S:38.1
                               a:52.400000000000006
                               b:0.0036729999999999996
                               c:-750700.0
                               d:-51.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:1.6764
                              a0:6.65e-05
                               K:1160000.0])) {
        [self setPhaseFormula:@"CaO"];
        [self setPhaseName:@"Lime"];
    }
    return self;
}

@end

@implementation RutileHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-944180.0
                               S:50.6
                               a:90.39999999999999
                               b:0.0029000000000000002
                               c:0.0
                               d:-623.8000000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:1.882
                              a0:4.43e-05
                               K:2225000.0])) {
        [self setPhaseFormula:@"TiO2"];
        [self setPhaseName:@"Rutile"];
    }
    return self;
}

@end

@implementation PericlaseHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-601600.0
                               S:26.9
                               a:60.5
                               b:0.00036199999999999996
                               c:-535800.0
                               d:-299.20000000000005
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:1.125
                              a0:6.2e-05
                               K:1650000.0])) {
        [self setPhaseFormula:@"MgO"];
        [self setPhaseName:@"Periclase"];
    }
    return self;
}

@end

@implementation ManganositeHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-385150.0
                               S:59.7
                               a:59.8
                               b:0.0036
                               c:-31400.0
                               d:-282.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:1.3221
                              a0:6.3e-05
                               K:1640000.0])) {
        [self setPhaseFormula:@"MnO"];
        [self setPhaseName:@"Manganosite"];
    }
    return self;
}

@end

@implementation CorundumHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1675250.0
                               S:50.9
                               a:139.5
                               b:0.00589
                               c:-2460600.0
                               d:-589.1999999999999
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.558
                              a0:4.19e-05
                               K:2520000.0])) {
        [self setPhaseFormula:@"Al2O3"];
        [self setPhaseName:@"Corundum"];
    }
    return self;
}

@end

@implementation HematiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-825710.0
                               S:87.4
                               a:163.89999999999998
                               b:0.0
                               c:-2257200.0
                               d:-657.5999999999999
                             Tc0:955.0
                            Smax:15.6
                            Vmax:0.0
                              v0:3.0274
                              a0:5.99e-05
                               K:1996000.0])) {
        [self setPhaseFormula:@"Fe2O3"];
        [self setPhaseName:@"Hematite"];
    }
    return self;
}

@end

@implementation Nickel_OxideHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-239440.0
                               S:38.0
                               a:47.7
                               b:0.007824
                               c:-392500.0
                               d:0.0
                             Tc0:520.0
                            Smax:5.7
                            Vmax:0.0
                              v0:1.097
                              a0:6.2e-05
                               K:1650000.0])) {
        [self setPhaseFormula:@"NiO"];
        [self setPhaseName:@"Nickel_Oxide"];
    }
    return self;
}

@end

@implementation PyrophaniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1359240.0
                               S:104.89999999999999
                               a:141.9
                               b:0.003373
                               c:-1940700.0
                               d:-407.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.288
                              a0:4.95e-05
                               K:1770000.0])) {
        [self setPhaseFormula:@"MnTiO3"];
        [self setPhaseName:@"Pyrophanite"];
    }
    return self;
}

@end

@implementation GeikieliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1567490.0
                               S:74.6
                               a:151.0
                               b:0.0
                               c:-1890400.0
                               d:-652.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.086
                              a0:4.95e-05
                               K:1770000.0])) {
        [self setPhaseFormula:@"MgTiO3"];
        [self setPhaseName:@"Geikielite"];
    }
    return self;
}

@end

@implementation IlmeniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1231300.0
                               S:108.89999999999999
                               a:138.9
                               b:0.0050810000000000004
                               c:-1288800.0
                               d:-463.7
                             Tc0:1900.0
                            Smax:11.0
                            Vmax:0.02
                              v0:3.169
                              a0:4.95e-05
                               K:1770000.0])) {
        [self setPhaseFormula:@"FeTiO3"];
        [self setPhaseName:@"Ilmenite"];
    }
    return self;
}

@end

@implementation BaddeleyiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1101030.0
                               S:50.4
                               a:90.7
                               b:0.0
                               c:-813300.0
                               d:-438.8
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.115
                              a0:3.76e-05
                               K:2225000.0])) {
        [self setPhaseFormula:@"ZrO2"];
        [self setPhaseName:@"Baddeleyite"];
    }
    return self;
}

@end

@implementation SpinelHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2300720.0
                               S:81.5
                               a:242.7
                               b:-0.006037
                               c:-2315100.0
                               d:-1678.1
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.978
                              a0:4.31e-05
                               K:1945000.0])) {
        [self setPhaseFormula:@"MgAl2O4"];
        [self setPhaseName:@"Spinel"];
    }
    return self;
}

@end

@implementation HercyniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1959150.0
                               S:107.5
                               a:283.3
                               b:-0.005376
                               c:609800.0
                               d:-2713.6
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.075
                              a0:3.95e-05
                               K:2120000.0])) {
        [self setPhaseFormula:@"FeAl2O4"];
        [self setPhaseName:@"Hercynite"];
    }
    return self;
}

@end

@implementation MagnetiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1115510.0
                               S:146.1
                               a:262.5
                               b:-0.007204
                               c:-1926200.0
                               d:-1655.7
                             Tc0:848.0
                            Smax:35.0
                            Vmax:0.0
                              v0:4.452
                              a0:6.96e-05
                               K:1850000.0])) {
        [self setPhaseFormula:@"Fe3O4"];
        [self setPhaseName:@"Magnetite"];
    }
    return self;
}

@end

@implementation MagnesioferriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1440660.0
                               S:126.5
                               a:217.9
                               b:0.000355
                               c:-3108000.0
                               d:-745.9
                             Tc0:665.0
                            Smax:12.9
                            Vmax:0.0
                              v0:4.457
                              a0:6.96e-05
                               K:1850000.0])) {
        [self setPhaseFormula:@"MgFe2O4"];
        [self setPhaseName:@"Magnesioferrite"];
    }
    return self;
}

@end

@implementation UlvospinelHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1497490.0
                               S:175.0
                               a:-102.6
                               b:0.14252
                               c:-9144500.0
                               d:5270.7
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.682
                              a0:6.9e-05
                               K:1850000.0])) {
        [self setPhaseFormula:@"Fe2TiO4"];
        [self setPhaseName:@"Ulvospinel"];
    }
    return self;
}

@end

@implementation BruciteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-924920.0
                               S:64.5
                               a:158.4
                               b:-0.004076
                               c:-1052300.0
                               d:-1171.3
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.463
                              a0:0.00013
                               K:485000.0])) {
        [self setPhaseFormula:@"MgO2H2"];
        [self setPhaseName:@"Brucite"];
    }
    return self;
}

@end

@implementation DiasporeHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-999470.0
                               S:35.0
                               a:145.1
                               b:0.008709000000000001
                               c:584400.0
                               d:-1741.1000000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:1.776
                              a0:7.97e-05
                               K:2300000.0])) {
        [self setPhaseFormula:@"AlO2H"];
        [self setPhaseName:@"Diaspore"];
    }
    return self;
}

@end

@implementation GoethiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-561650.0
                               S:60.400000000000006
                               a:139.3
                               b:0.000147
                               c:-212700.0
                               d:-1077.8000000000002
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.082
                              a0:7.97e-05
                               K:2300000.0])) {
        [self setPhaseFormula:@"FeO2H"];
        [self setPhaseName:@"Goethite"];
    }
    return self;
}

@end

@implementation CalciteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1207470.0
                               S:92.5
                               a:140.9
                               b:0.005029
                               c:-950700.0
                               d:-858.4000000000001
                             Tc0:1240.0
                            Smax:10.0
                            Vmax:0.04
                              v0:3.689
                              a0:4.4e-05
                               K:760000.0])) {
        [self setPhaseFormula:@"CaCO3"];
        [self setPhaseName:@"Calcite"];
    }
    return self;
}

@end

@implementation AragoniteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1207580.0
                               S:89.5
                               a:192.3
                               b:-0.003052
                               c:1149700.0
                               d:-2118.3
                             Tc0:1240.0
                            Smax:9.0
                            Vmax:0.04
                              v0:3.415
                              a0:0.000115
                               K:650000.0])) {
        [self setPhaseFormula:@"CaO3C"];
        [self setPhaseName:@"Aragonite"];
    }
    return self;
}

@end

@implementation MagnesiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1111360.0
                               S:65.10000000000001
                               a:186.4
                               b:-0.003772
                               c:0.0
                               d:-1886.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.803
                              a0:6.48e-05
                               K:1460000.0])) {
        [self setPhaseFormula:@"MgO3C"];
        [self setPhaseName:@"Magnesite"];
    }
    return self;
}

@end

@implementation SideriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-761440.0
                               S:95.0
                               a:168.4
                               b:0.0
                               c:0.0
                               d:-1483.6000000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.938
                              a0:0.00011
                               K:1200000.0])) {
        [self setPhaseFormula:@"FeCO3"];
        [self setPhaseName:@"Siderite"];
    }
    return self;
}

@end

@implementation RhodochrositeHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-891090.0
                               S:98.0
                               a:169.5
                               b:0.0
                               c:0.0
                               d:-1534.3
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.107
                              a0:6.5e-05
                               K:800000.0])) {
        [self setPhaseFormula:@"MnCO3"];
        [self setPhaseName:@"Rhodochrosite"];
    }
    return self;
}

@end

@implementation DolomiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2324430.0
                               S:156.0
                               a:358.9
                               b:-0.004905
                               c:0.0
                               d:-3456.2
                             Tc0:1373.0
                            Smax:13.0
                            Vmax:0.015
                              v0:6.434
                              a0:6.35e-05
                               K:900000.0])) {
        [self setPhaseFormula:@"CaMgO6C2"];
        [self setPhaseName:@"Dolomite"];
    }
    return self;
}

@end

@implementation AnkeriteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1971270.0
                               S:187.0
                               a:341.0
                               b:-0.001161
                               c:0.0
                               d:-3054.8
                             Tc0:1273.0
                            Smax:9.0
                            Vmax:0.01
                              v0:6.606
                              a0:6.35e-05
                               K:900000.0])) {
        [self setPhaseFormula:@"CaFeO6C2"];
        [self setPhaseName:@"Ankerite"];
    }
    return self;
}

@end

@implementation SylviteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-436500.0
                               S:82.60000000000001
                               a:46.199999999999996
                               b:0.01797
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.752
                              a0:0.000247
                               K:170000.0])) {
        [self setPhaseFormula:@"KCl"];
        [self setPhaseName:@"Sylvite"];
    }
    return self;
}

@end

@implementation HaliteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-411300.0
                               S:72.1
                               a:45.199999999999996
                               b:0.01797
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.702
                              a0:0.000269
                               K:240000.0])) {
        [self setPhaseFormula:@"NaCl"];
        [self setPhaseName:@"Halite"];
    }
    return self;
}

@end

@implementation IronHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-0.0
                               S:27.32
                               a:46.199999999999996
                               b:0.005158
                               c:723100.0
                               d:-556.2
                             Tc0:1042.0
                            Smax:8.3
                            Vmax:0.0
                              v0:0.7092
                              a0:7.46e-05
                               K:1680000.0])) {
        [self setPhaseFormula:@"Fe"];
        [self setPhaseName:@"Iron"];
    }
    return self;
}

@end

@implementation NickelHollandAndPowell

-(id)init {
    if ((self = [super initWithH:0.0
                               S:29.87
                               a:49.8
                               b:0.0
                               c:585900.0
                               d:-533.9000000000001
                             Tc0:631.0
                            Smax:3.0
                            Vmax:0.0
                              v0:0.6588
                              a0:8.86e-05
                               K:1870000.0])) {
        [self setPhaseFormula:@"Ni"];
        [self setPhaseName:@"Nickel"];
    }
    return self;
}

@end

@implementation GraphiteHollandAndPowell

-(id)init {
    if ((self = [super initWithH:0.0
                               S:5.8500000000000005
                               a:51.0
                               b:-0.004428
                               c:488600.0
                               d:-805.5
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:0.5298
                              a0:4.84e-05
                               K:390000.0])) {
        [self setPhaseFormula:@"C"];
        [self setPhaseName:@"Graphite"];
    }
    return self;
}

@end

@implementation DiamondHollandAndPowell

-(id)init {
    if ((self = [super initWithH:2070.0
                               S:2.3
                               a:24.299999999999997
                               b:0.006272000000000001
                               c:-377400.0
                               d:-273.4
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:0.3417
                              a0:1.65e-05
                               K:5800000.0])) {
        [self setPhaseFormula:@"C"];
        [self setPhaseName:@"Diamond"];
    }
    return self;
}

@end

@implementation Water_FluidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-241810.0
                               S:188.79999999999998
                               a:40.099999999999994
                               b:0.008656
                               c:487500.0
                               d:-251.2
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:0.0
                              a0:0.0
                               K:0.0])) {
        [self setPhaseFormula:@"H2O"];
        [self setPhaseName:@"Water_Fluid"];
    }
    return self;
}

@end

@implementation Carbon_DioxideHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-393510.0
                               S:213.7
                               a:87.8
                               b:-0.002644
                               c:706400.0
                               d:-998.9
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:0.0
                              a0:0.0
                               K:0.0])) {
        [self setPhaseFormula:@"CO2"];
        [self setPhaseName:@"Carbon_Dioxide"];
    }
    return self;
}

@end

@implementation Carbon_MonoxideHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-110530.0
                               S:197.67000000000002
                               a:45.699999999999996
                               b:-9.7e-05
                               c:662700.0
                               d:-414.7
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:0.0
                              a0:0.0
                               K:0.0])) {
        [self setPhaseFormula:@"CO"];
        [self setPhaseName:@"Carbon_Monoxide"];
    }
    return self;
}

@end

@implementation MethaneHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-74810.0
                               S:186.26000000000002
                               a:150.10000000000002
                               b:0.002062
                               c:3427700.0
                               d:-2650.4
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:0.0
                              a0:0.0
                               K:0.0])) {
        [self setPhaseFormula:@"CH4"];
        [self setPhaseName:@"Methane"];
    }
    return self;
}

@end

@implementation OxygenHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-0.0
                               S:205.2
                               a:48.300000000000004
                               b:-0.000691
                               c:499200.0
                               d:-420.70000000000005
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:0.0
                              a0:0.0
                               K:0.0])) {
        [self setPhaseFormula:@"O2"];
        [self setPhaseName:@"Oxygen"];
    }
    return self;
}

@end

@implementation HydrogenHollandAndPowell

-(id)init {
    if ((self = [super initWithH:0.0
                               S:130.70000000000002
                               a:23.3
                               b:0.004627
                               c:0.0
                               d:76.30000000000001
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:0.0
                              a0:0.0
                               K:0.0])) {
        [self setPhaseFormula:@"H2"];
        [self setPhaseName:@"Hydrogen"];
    }
    return self;
}

@end

@implementation Sylvite_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-416490.0
                               S:95.3
                               a:66.9
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:3.857
                              a0:0.0005
                               K:59000.0])) {
        [self setPhaseFormula:@"KCl"];
        [self setPhaseName:@"Sylvite_Liquid"];
    }
    return self;
}

@end

@implementation Halite_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-393420.0
                               S:79.69999999999999
                               a:72.0
                               b:-0.0032229999999999997
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.965
                              a0:0.0005
                               K:66000.0])) {
        [self setPhaseFormula:@"NaCl"];
        [self setPhaseName:@"Halite_Liquid"];
    }
    return self;
}

@end

@implementation Forsterite_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2225160.0
                               S:-55.0
                               a:267.90000000000003
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.243
                              a0:0.000145
                               K:730000.0])) {
        [self setPhaseFormula:@"Mg2SiO4"];
        [self setPhaseName:@"Forsterite_Liquid"];
    }
    return self;
}

@end

@implementation Fayalite_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-1459210.0
                               S:102.5
                               a:239.7
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:4.695
                              a0:0.000169
                               K:410000.0])) {
        [self setPhaseFormula:@"Fe2SiO4"];
        [self setPhaseName:@"Fayalite_Liquid"];
    }
    return self;
}

@end

@implementation Sillimanite_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-2583100.0
                               S:39.0
                               a:237.6
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.419
                              a0:1e-05
                               K:300000.0])) {
        [self setPhaseFormula:@"Al2SiO5"];
        [self setPhaseName:@"Sillimanite_Liquid"];
    }
    return self;
}

@end

@implementation Anorthite_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-4257750.0
                               S:52.0
                               a:417.5
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.206
                              a0:4.9e-05
                               K:200000.0])) {
        [self setPhaseFormula:@"CaAl2Si2O8"];
        [self setPhaseName:@"Anorthite_Liquid"];
    }
    return self;
}

@end

@implementation H2O_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-295710.0
                               S:45.5
                               a:80.0
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:1.414
                              a0:0.001079
                               K:40000.0])) {
        [self setPhaseFormula:@"H2O"];
        [self setPhaseName:@"H2O_Liquid"];
    }
    return self;
}

@end

@implementation Enstatite_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3091760.0
                               S:-2.0
                               a:354.9
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:6.9
                              a0:0.000129
                               K:260000.0])) {
        [self setPhaseFormula:@"Mg2Si2O6"];
        [self setPhaseName:@"Enstatite_Liquid"];
    }
    return self;
}

@end

@implementation K_Feldspar_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3992670.0
                               S:129.5
                               a:367.3
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:11.468
                              a0:6e-05
                               K:260000.0])) {
        [self setPhaseFormula:@"KAlSi3O8"];
        [self setPhaseName:@"K_Feldspar_Liquid"];
    }
    return self;
}

@end

@implementation Silica_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-920850.0
                               S:16.5
                               a:82.5
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:2.64
                              a0:-5e-06
                               K:470000.0])) {
        [self setPhaseFormula:@"SiO2"];
        [self setPhaseName:@"Silica_Liquid"];
    }
    return self;
}

@end

@implementation Diopside_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3208180.0
                               S:23.8
                               a:345.3
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:7.34
                              a0:0.000129
                               K:300000.0])) {
        [self setPhaseFormula:@"CaMgSi2O6"];
        [self setPhaseName:@"Diopside_Liquid"];
    }
    return self;
}

@end

@implementation Albite_LiquidHollandAndPowell

-(id)init {
    if ((self = [super initWithH:-3934370.0
                               S:145.0
                               a:358.5
                               b:0.0
                               c:0.0
                               d:0.0
                             Tc0:0.0
                            Smax:0.0
                            Vmax:0.0
                              v0:10.71
                              a0:4.5e-05
                               K:390000.0])) {
        [self setPhaseFormula:@"NaAlSi3O8"];
        [self setPhaseName:@"Albite_Liquid"];
    }
    return self;
}

@end
