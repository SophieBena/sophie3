//
//  StixrudeStoichiometricPhases.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 3/10/11.
//  Higher precision parameters added by Bob Myhill on 07/16/18.
//  Copyright 2018 OFM Research Inc. All rights reserved.
//
//  Stixrude and Lithgow-Bertelloni (2011) database
//

#import "StixrudeStoichiometricPhases.h"

@implementation AnorthiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-4014.619],
										  [NSNumber numberWithDouble:13.0],
										  [NSNumber numberWithDouble:100.61],
										  [NSNumber numberWithDouble:84.08915],
										  [NSNumber numberWithDouble:4.0],
										  [NSNumber numberWithDouble:752.3911],
										  [NSNumber numberWithDouble:0.39241],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"CaAl2Si2O8"];
		[self setPhaseName:@"Anorthite"];
	}
	return self;
}
@end

@implementation AlbiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-3718.799],
										  [NSNumber numberWithDouble:13.0],
										  [NSNumber numberWithDouble:100.452],
										  [NSNumber numberWithDouble:59.76162],
										  [NSNumber numberWithDouble:4.0],
										  [NSNumber numberWithDouble:713.7824],
										  [NSNumber numberWithDouble:0.56704],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"NaAlSi3O8"];
		[self setPhaseName:@"Albite"];
	}
	return self;
}
@end

@implementation MgSpinelStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-8667.568],
										  [NSNumber numberWithDouble:28.0],
										  [NSNumber numberWithDouble:159.048],
										  [NSNumber numberWithDouble:196.9428],
										  [NSNumber numberWithDouble:5.68282],
										  [NSNumber numberWithDouble:842.8104],
										  [NSNumber numberWithDouble:1.02283],
										  [NSNumber numberWithDouble:2.71208],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg4Al8O16"];
		[self setPhaseName:@"MgSpinel"];
	}
	return self;
}
@end

@implementation HercyniteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-7324.009],
										  [NSNumber numberWithDouble:28.0],
										  [NSNumber numberWithDouble:163.372],
										  [NSNumber numberWithDouble:208.8965],
										  [NSNumber numberWithDouble:5.68282],
										  [NSNumber numberWithDouble:763.231],
										  [NSNumber numberWithDouble:1.21719],
										  [NSNumber numberWithDouble:2.71208],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Fe4Al8O16"];
		[self setPhaseName:@"Hercynite"];
	}
	return self;
}
@end

@implementation ForsteriteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2055.403],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:43.603],
										  [NSNumber numberWithDouble:127.9555],
										  [NSNumber numberWithDouble:4.21796],
										  [NSNumber numberWithDouble:809.1703],
										  [NSNumber numberWithDouble:0.99282],
										  [NSNumber numberWithDouble:2.10672],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg2SiO4"];
		[self setPhaseName:@"Forsterite"];
	}
	return self;
}
@end

@implementation FayaliteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1370.519],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:46.29],
										  [NSNumber numberWithDouble:134.9622],
										  [NSNumber numberWithDouble:4.21796],
										  [NSNumber numberWithDouble:618.7007],
										  [NSNumber numberWithDouble:1.06023],
										  [NSNumber numberWithDouble:3.6466],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Fe2SiO4"];
		[self setPhaseName:@"Fayalite"];
	}
	return self;
}
@end

@implementation MgWadsleyiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2027.837],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:40.515],
										  [NSNumber numberWithDouble:168.6948],
										  [NSNumber numberWithDouble:4.3229],
										  [NSNumber numberWithDouble:843.4973],
										  [NSNumber numberWithDouble:1.2061],
										  [NSNumber numberWithDouble:2.0188],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg2SiO4"];
		[self setPhaseName:@"MgWadsleyite"];
	}
	return self;
}
@end

@implementation FeWadsleyiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1364.668],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:42.8],
										  [NSNumber numberWithDouble:168.591],
										  [NSNumber numberWithDouble:4.3229],
										  [NSNumber numberWithDouble:665.4492],
										  [NSNumber numberWithDouble:1.2061],
										  [NSNumber numberWithDouble:2.0188],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Fe2SiO4"];
		[self setPhaseName:@"FeWadsleyite"];
	}
	return self;
}
@end

@implementation MgRingwooditeStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2017.557],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:39.493],
										  [NSNumber numberWithDouble:184.9009],
										  [NSNumber numberWithDouble:4.22035],
										  [NSNumber numberWithDouble:877.7094],
										  [NSNumber numberWithDouble:1.10791],
										  [NSNumber numberWithDouble:2.3914],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg2SiO4"];
		[self setPhaseName:@"MgRingwoodite"];
	}
	return self;
}
@end

@implementation FeRingwooditeStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1362.772],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:41.86],
										  [NSNumber numberWithDouble:213.412],
										  [NSNumber numberWithDouble:4.22035],
										  [NSNumber numberWithDouble:677.7177],
										  [NSNumber numberWithDouble:1.27193],
										  [NSNumber numberWithDouble:2.3914],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Fe2SiO4"];
		[self setPhaseName:@"FeRingwoodite"];
	}
	return self;
}
@end

@implementation EnstatiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2913.596],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:62.676],
										  [NSNumber numberWithDouble:107.0768],
										  [NSNumber numberWithDouble:7.02751],
										  [NSNumber numberWithDouble:812.1848],
										  [NSNumber numberWithDouble:0.78479],
										  [NSNumber numberWithDouble:3.43846],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg2Si2O6"];
		[self setPhaseName:@"Enstatite"];
	}
	return self;
}
@end

@implementation FerrosiliteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2225.718],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:65.941],
										  [NSNumber numberWithDouble:100.5386],
										  [NSNumber numberWithDouble:7.02751],
										  [NSNumber numberWithDouble:674.4769],
										  [NSNumber numberWithDouble:0.71889],
										  [NSNumber numberWithDouble:3.43846],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Fe2Si2O6"];
		[self setPhaseName:@"Ferrosilite"];
	}
	return self;
}
@end

@implementation MgTschermaksStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-3002.47],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:59.14],
										  [NSNumber numberWithDouble:107.0768],
										  [NSNumber numberWithDouble:7.02751],
										  [NSNumber numberWithDouble:783.8404],
										  [NSNumber numberWithDouble:0.78479],
										  [NSNumber numberWithDouble:3.43846],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"MgAl2SiO6"];
		[self setPhaseName:@"MgTschermaks"];
	}
	return self;
}
@end

@implementation OrthoDiopsideStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-3015.827],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:68.054],
										  [NSNumber numberWithDouble:107.0768],
										  [NSNumber numberWithDouble:7.02751],
										  [NSNumber numberWithDouble:744.6988],
										  [NSNumber numberWithDouble:0.78479],
										  [NSNumber numberWithDouble:3.43846],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"CaMgSi2O6"];
		[self setPhaseName:@"OrthoDiopside"];
	}
	return self;
}
@end

@implementation DiopsideStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-3029.531],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:66.039],
										  [NSNumber numberWithDouble:112.2413],
										  [NSNumber numberWithDouble:5.23885],
										  [NSNumber numberWithDouble:781.6146],
										  [NSNumber numberWithDouble:0.95873],
										  [NSNumber numberWithDouble:1.52852],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"CaMgSi2O6"];
		[self setPhaseName:@"Diopside"];
	}
	return self;
}
@end

@implementation HedenbergiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2677.33],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:67.867],
										  [NSNumber numberWithDouble:119.2555],
										  [NSNumber numberWithDouble:5.23885],
										  [NSNumber numberWithDouble:701.5851],
										  [NSNumber numberWithDouble:0.93516],
										  [NSNumber numberWithDouble:1.52852],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"CaFeSi2O6"];
		[self setPhaseName:@"Hedenbergite"];
	}
	return self;
}
@end

@implementation ClinoenstatiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2905.918],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:62.5],
										  [NSNumber numberWithDouble:112.2413],
										  [NSNumber numberWithDouble:5.23885],
										  [NSNumber numberWithDouble:805.0547],
										  [NSNumber numberWithDouble:0.95873],
										  [NSNumber numberWithDouble:1.52852],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg2Si2O6"];
		[self setPhaseName:@"Clinoenstatite"];
	}
	return self;
}
@end

@implementation CaTschermaksStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-3120.253],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:63.574],
										  [NSNumber numberWithDouble:112.2413],
										  [NSNumber numberWithDouble:5.23885],
										  [NSNumber numberWithDouble:803.6626],
										  [NSNumber numberWithDouble:0.78126],
										  [NSNumber numberWithDouble:1.52852],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"CaAl2SiO6"];
		[self setPhaseName:@"CaTschermaks"];
	}
	return self;
}
@end

@implementation JadeiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2855.192],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:60.508],
										  [NSNumber numberWithDouble:142.2873],
										  [NSNumber numberWithDouble:5.23885],
										  [NSNumber numberWithDouble:820.7623],
										  [NSNumber numberWithDouble:0.903],
										  [NSNumber numberWithDouble:0.39234],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"NaAlSi2O6"];
		[self setPhaseName:@"Jadeite"];
	}
	return self;
}
@end

@implementation MgHpCpxStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2905.788],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:60.76],
										  [NSNumber numberWithDouble:116.0254],
										  [NSNumber numberWithDouble:6.23685],
										  [NSNumber numberWithDouble:824.4439],
										  [NSNumber numberWithDouble:1.12473],
										  [NSNumber numberWithDouble:0.20401],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg2Si2O6"];
		[self setPhaseName:@"MgHpCpx"];
	}
	return self;
}
@end

@implementation FeHpCpxStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2222.183],
										  [NSNumber numberWithDouble:10.0],
										  [NSNumber numberWithDouble:63.85413],
										  [NSNumber numberWithDouble:116.0254],
										  [NSNumber numberWithDouble:6.23685],
										  [NSNumber numberWithDouble:691.564],
										  [NSNumber numberWithDouble:1.12473],
										  [NSNumber numberWithDouble:0.20401],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Fe2Si2O6"];
		[self setPhaseName:@"FeHpCpx"];
	}
	return self;
}
@end

@implementation CaPerovskiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1463.358],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:27.45],
										  [NSNumber numberWithDouble:236.0],
										  [NSNumber numberWithDouble:3.9],
										  [NSNumber numberWithDouble:795.779],
										  [NSNumber numberWithDouble:1.88839],
										  [NSNumber numberWithDouble:0.89769],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"CaSiO3"];
		[self setPhaseName:@"CaPerovskite"];
	}
	return self;
}
@end

@implementation MgAkimotoiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1410.85],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:26.354],
										  [NSNumber numberWithDouble:210.706],
										  [NSNumber numberWithDouble:5.62088],
										  [NSNumber numberWithDouble:935.9778],
										  [NSNumber numberWithDouble:1.18984],
										  [NSNumber numberWithDouble:2.34514],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"MgSiO3"];
		[self setPhaseName:@"MgAkimotoite"];
	}
	return self;
}
@end

@implementation FeAkimotoiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1067.598],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:26.854],
										  [NSNumber numberWithDouble:210.706],
										  [NSNumber numberWithDouble:5.62088],
										  [NSNumber numberWithDouble:887.8709],
										  [NSNumber numberWithDouble:1.18984],
										  [NSNumber numberWithDouble:2.34514],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"FeSiO3"];
		[self setPhaseName:@"FeAkimotoite"];
	}
	return self;
}
@end

@implementation AlAkimotoiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1582.454],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:25.577],
										  [NSNumber numberWithDouble:252.5457],
										  [NSNumber numberWithDouble:4.33728],
										  [NSNumber numberWithDouble:932.5696],
										  [NSNumber numberWithDouble:1.32442],
										  [NSNumber numberWithDouble:1.30316],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Al2O3"];
		[self setPhaseName:@"AlAkimotoite"];
	}
	return self;
}
@end

@implementation PyropeStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-5936.538],
										  [NSNumber numberWithDouble:20.0],
										  [NSNumber numberWithDouble:113.08],
										  [NSNumber numberWithDouble:170.2396],
										  [NSNumber numberWithDouble:4.11067],
										  [NSNumber numberWithDouble:823.2102],
										  [NSNumber numberWithDouble:1.01424],
										  [NSNumber numberWithDouble:1.42169],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg3Al2Si3O12"];
		[self setPhaseName:@"Pyrope"];
	}
	return self;
}
@end

@implementation AlmandineStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-4935.516],
										  [NSNumber numberWithDouble:20.0],
										  [NSNumber numberWithDouble:115.43],
										  [NSNumber numberWithDouble:173.8963],
										  [NSNumber numberWithDouble:4.91341],
										  [NSNumber numberWithDouble:741.356],
										  [NSNumber numberWithDouble:1.06495],
										  [NSNumber numberWithDouble:1.42169],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Fe3Al2Si3O12"];
		[self setPhaseName:@"Almandine"];
	}
	return self;
}
@end

@implementation GrossularStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-6277.935],
										  [NSNumber numberWithDouble:20.0],
										  [NSNumber numberWithDouble:125.12],
										  [NSNumber numberWithDouble:167.0622],
										  [NSNumber numberWithDouble:3.91544],
										  [NSNumber numberWithDouble:822.743],
										  [NSNumber numberWithDouble:1.05404],
										  [NSNumber numberWithDouble:1.88887],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Ca3Al2Si3O12"];
		[self setPhaseName:@"Grossular"];
	}
	return self;
}
@end

@implementation MajoriteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-5691.614],
										  [NSNumber numberWithDouble:20.0],
										  [NSNumber numberWithDouble:114.324],
										  [NSNumber numberWithDouble:165.1183],
										  [NSNumber numberWithDouble:4.21183],
										  [NSNumber numberWithDouble:822.458],
										  [NSNumber numberWithDouble:0.97682],
										  [NSNumber numberWithDouble:1.53581],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Mg4Si4O12"];
		[self setPhaseName:@"Majorite"];
	}
	return self;
}
@end

@implementation NaMajoriteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-5518.542],
										  [NSNumber numberWithDouble:20.0],
										  [NSNumber numberWithDouble:110.94],
										  [NSNumber numberWithDouble:177.0772],
										  [NSNumber numberWithDouble:4.11067],
										  [NSNumber numberWithDouble:895.914],
										  [NSNumber numberWithDouble:1.01424],
										  [NSNumber numberWithDouble:1.42169],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Na2Al2Si4O12"];
		[self setPhaseName:@"NaMajorite"];
	}
	return self;
}
@end

@implementation QuartzStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-858.8534],
										  [NSNumber numberWithDouble:3.0],
										  [NSNumber numberWithDouble:23.67003],
										  [NSNumber numberWithDouble:49.54743],
										  [NSNumber numberWithDouble:4.33155],
										  [NSNumber numberWithDouble:816.3307],
										  [NSNumber numberWithDouble:-0.00296],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]
						  andLandauTerms:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:847.0],
										  [NSNumber numberWithDouble:0.1222],
										  [NSNumber numberWithDouble:5.164],
										  nil]])) {
		[self setPhaseFormula:@"SiO2"];
		[self setPhaseName:@"Quartz"];
	}
	return self;
}
@end

@implementation CoesiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-855.0685],
										  [NSNumber numberWithDouble:3.0],
										  [NSNumber numberWithDouble:20.657],
										  [NSNumber numberWithDouble:113.5856],
										  [NSNumber numberWithDouble:4.0],
										  [NSNumber numberWithDouble:852.4267],
										  [NSNumber numberWithDouble:0.39157],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"SiO2"];
		[self setPhaseName:@"Coesite"];
	}
	return self;
}
@end

@implementation StishoviteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-818.9846],
										  [NSNumber numberWithDouble:3.0],
										  [NSNumber numberWithDouble:14.017],
										  [NSNumber numberWithDouble:314.3352],
										  [NSNumber numberWithDouble:3.75122],
										  [NSNumber numberWithDouble:1107.824],
										  [NSNumber numberWithDouble:1.37466],
										  [NSNumber numberWithDouble:2.83517],
										  [NSNumber numberWithDouble:0.0],
										  nil]
						  andLandauTerms:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-4250.0],
										  [NSNumber numberWithDouble:0.0001],
										  [NSNumber numberWithDouble:0.012],
										  nil]])) {
		[self setPhaseFormula:@"SiO2"];
		[self setPhaseName:@"Stishovite"];
	}
	return self;
}
@end

@implementation SeifertiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-794.3354],
										  [NSNumber numberWithDouble:3.0],
										  [NSNumber numberWithDouble:13.67],
										  [NSNumber numberWithDouble:327.5843],
										  [NSNumber numberWithDouble:4.01553],
										  [NSNumber numberWithDouble:1140.772],
										  [NSNumber numberWithDouble:1.37466],
										  [NSNumber numberWithDouble:2.83517],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"SiO2"];
		[self setPhaseName:@"Seifertite"];
	}
	return self;
}
@end

@implementation MgPerovskiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1368.283],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:24.445],
										  [NSNumber numberWithDouble:250.5264],
										  [NSNumber numberWithDouble:4.14],
										  [NSNumber numberWithDouble:905.9412],
										  [NSNumber numberWithDouble:1.56508],
										  [NSNumber numberWithDouble:1.10945],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"MgSiO3"];
		[self setPhaseName:@"MgPerovskite"];
	}
	return self;
}
@end

@implementation FePerovskiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1040.92],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:25.485],
										  [NSNumber numberWithDouble:272.1152],
										  [NSNumber numberWithDouble:4.14],
										  [NSNumber numberWithDouble:870.8122],
										  [NSNumber numberWithDouble:1.56508],
										  [NSNumber numberWithDouble:1.10945],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"FeSiO3"];
		[self setPhaseName:@"FePerovskite"];
	}
	return self;
}
@end

@implementation AlPerovskiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1533.878],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:24.944],
										  [NSNumber numberWithDouble:258.2],
										  [NSNumber numberWithDouble:4.14],
										  [NSNumber numberWithDouble:886.4601],
										  [NSNumber numberWithDouble:1.56508],
										  [NSNumber numberWithDouble:1.10945],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Al2O3"];
		[self setPhaseName:@"AlPerovskite"];
	}
	return self;
}
@end

@implementation MgPostPerovskiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1348.641],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:24.419],
										  [NSNumber numberWithDouble:231.2],
										  [NSNumber numberWithDouble:4.0],
										  [NSNumber numberWithDouble:855.8173],
										  [NSNumber numberWithDouble:1.89155],
										  [NSNumber numberWithDouble:1.09081],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"MgSiO3"];
		[self setPhaseName:@"MgPostPerovskite"];
	}
	return self;
}
@end

@implementation FePostPerovskiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-981.8069],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:25.459],
										  [NSNumber numberWithDouble:231.2],
										  [NSNumber numberWithDouble:4.0],
										  [NSNumber numberWithDouble:781.3465],
										  [NSNumber numberWithDouble:1.89155],
										  [NSNumber numberWithDouble:1.09081],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"FeSiO3"];
		[self setPhaseName:@"FePostPerovskite"];
	}
	return self;
}
@end

@implementation AlPostPerovskiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1377.582],
										  [NSNumber numberWithDouble:5.0],
										  [NSNumber numberWithDouble:23.847],
										  [NSNumber numberWithDouble:249.0],
										  [NSNumber numberWithDouble:4.0],
										  [NSNumber numberWithDouble:762.1951],
										  [NSNumber numberWithDouble:1.64573],
										  [NSNumber numberWithDouble:1.09081],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Al2O3"];
		[self setPhaseName:@"AlPostPerovskite"];
	}
	return self;
}
@end

@implementation PericlaseStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-569.4446],
										  [NSNumber numberWithDouble:2.0],
										  [NSNumber numberWithDouble:11.244],
										  [NSNumber numberWithDouble:161.3836],
										  [NSNumber numberWithDouble:3.84045],
										  [NSNumber numberWithDouble:767.0977],
										  [NSNumber numberWithDouble:1.36127],
										  [NSNumber numberWithDouble:1.7217],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"MgO"];
		[self setPhaseName:@"Periclase"];
	}
	return self;
}
@end

@implementation WuestiteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-242.146],
										  [NSNumber numberWithDouble:2.0],
										  [NSNumber numberWithDouble:12.264],
										  [NSNumber numberWithDouble:179.4442],
										  [NSNumber numberWithDouble:4.9376],
										  [NSNumber numberWithDouble:454.1592],
										  [NSNumber numberWithDouble:1.53047],
										  [NSNumber numberWithDouble:1.7217],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"FeO"];
		[self setPhaseName:@"Wuestite"];
	}
	return self;
}
@end

@implementation MgCFStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2122.169],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:36.177],
										  [NSNumber numberWithDouble:210.6663],
										  [NSNumber numberWithDouble:4.0528],
										  [NSNumber numberWithDouble:838.6291],
										  [NSNumber numberWithDouble:1.31156],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"MgAl2O4"];
		[self setPhaseName:@"MgCF"];
	}
	return self;
}
@end

@implementation FeCFStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1790.284],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:37.258],
										  [NSNumber numberWithDouble:210.6663],
										  [NSNumber numberWithDouble:4.0528],
										  [NSNumber numberWithDouble:804.1986],
										  [NSNumber numberWithDouble:1.31156],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"FeAl2O4"];
		[self setPhaseName:@"FeCF"];
	}
	return self;
}
@end

@implementation NaCFStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1844.129],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:36.27],
										  [NSNumber numberWithDouble:161.3385],
										  [NSNumber numberWithDouble:4.32479],
										  [NSNumber numberWithDouble:812.4769],
										  [NSNumber numberWithDouble:0.69428],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"NaAlSiO4"];
		[self setPhaseName:@"NaCF"];
	}
	return self;
}
@end

@implementation KyaniteStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-2446.058],
										  [NSNumber numberWithDouble:8.0],
										  [NSNumber numberWithDouble:44.227],
										  [NSNumber numberWithDouble:160.0],
										  [NSNumber numberWithDouble:4.0],
										  [NSNumber numberWithDouble:943.1665],
										  [NSNumber numberWithDouble:0.9255],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"Al2SiO5"];
		[self setPhaseName:@"Kyanite"];
	}
	return self;
}
@end

@implementation NephelineStixrude
-(id)init {
	if ((self = [super initWithParameters:[NSArray arrayWithObjects:
										  [NSNumber numberWithDouble:-1992.104],
										  [NSNumber numberWithDouble:7.0],
										  [NSNumber numberWithDouble:54.6684],
										  [NSNumber numberWithDouble:53.07799],
										  [NSNumber numberWithDouble:4.0],
										  [NSNumber numberWithDouble:700.9422],
										  [NSNumber numberWithDouble:0.69428],
										  [NSNumber numberWithDouble:1.0],
										  [NSNumber numberWithDouble:0.0],
										  nil]])) {
		[self setPhaseFormula:@"NaAlSiO4"];
		[self setPhaseName:@"Nepheline"];
	}
	return self;
}
@end
