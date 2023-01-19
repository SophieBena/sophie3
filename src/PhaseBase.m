//
//  PhaseBase.m
//  PhasePlot
//
//  Created by Mark Ghiorso on 6/14/10.
//  Copyright 2010 OFM Research Inc.. All rights reserved.
//

#import "PhaseBase.h"
#import "DoubleVector.h"

@implementation PhaseBase

typedef struct _elementsData {
	const char *name;
	const char *sym;
	double     atomicWt;
} ElementsData;

static const ElementsData elements[] = {
	{ "dummy"           , "!" ,   0.0     },
	{ "hydrogen"        , "H" ,   1.0079  },
	{ "helium"          , "He",   4.00260 },
	{ "lithium"         , "Li",   6.94    },
	{ "beryllium"       , "Be",   9.01218 },
	{ "boron"           , "B" ,  10.81    },
	{ "carbon"          , "C" ,  12.011   },
	{ "nitrogen"        , "N" ,  14.0067  },
	{ "oxygen"          , "O" ,  15.9994  },
	{ "fluorine"        , "F" ,  18.998403},
	{ "neon"            , "Ne",  20.179   },
	{ "sodium"          , "Na",  22.98977 },
	{ "magnesium"       , "Mg",  24.305   },
	{ "aluminum"        , "Al",  26.98154 },
	{ "silicon"         , "Si",  28.0855  },
	{ "phosphorous"     , "P" ,  30.97376 },
	{ "sulfur"          , "S" ,  32.06    },
	{ "chlorine"        , "Cl",  35.453   },
	{ "argon"           , "Ar",  39.948   },
	{ "potassium"       , "K" ,  39.102   },
	{ "calcium"         , "Ca",  40.08    },
	{ "scandium"        , "Sc",  44.9559  },
	{ "titanium"        , "Ti",  47.90    },
	{ "vanadium"        , "V" ,  50.9415  },
	{ "chromium"        , "Cr",  51.996   },
	{ "manganese"       , "Mn",  54.9380  },
	{ "iron"            , "Fe",  55.847   },
	{ "cobalt"          , "Co",  58.9332  },
	{ "nickel"          , "Ni",  58.71    },
	{ "copper"          , "Cu",  63.546   },
	{ "zinc"            , "Zn",  65.38    },
	{ "gallium"         , "Ga",  69.735   },
	{ "germanium"       , "Ge",  72.59    },
	{ "arsenic"         , "As",  74.9216  },
	{ "selenium"        , "Se",  78.96    },
	{ "bromine"         , "Br",  79.904   },
	{ "krypton"         , "Kr",  83.80    },
	{ "rubidium"        , "Rb",  85.4678  },
	{ "strontium"       , "Sr",  87.62    },
	{ "yttrium"         , "Y" ,  88.9059  },
	{ "zirconium"       , "Zr",  91.22    },
	{ "niobium"         , "Nb",  92.9064  },
	{ "molybdenum"      , "Mo",  95.94    },
	{ "technetium"      , "Tc",  98.9062  },
	{ "ruthenium"       , "Ru", 101.07    },
	{ "rhodium"         , "Rh", 102.9055  },
	{ "palladium"       , "Pd", 106.4     },
	{ "silver"          , "Ag", 107.868   },
	{ "cadmium"         , "Cd", 112.41    },
	{ "indium"          , "In", 114.82    },
	{ "tin"             , "Sn", 118.69    },
	{ "antimony"        , "Sb", 121.75    },
	{ "tellurium"       , "Te", 127.60    },
	{ "iodine"          , "I" , 126.9045  },
	{ "xenon"           , "Xe", 131.30    },
	{ "cesium"          , "Cs", 132.9054  },
	{ "barium"          , "Ba", 137.33    },
	{ "lantahnum"       , "La", 138.9055  },
	{ "cerium"          , "Ce", 140.12    },
	{ "praseodymium"    , "Pr", 140.9077  },
	{ "neodymium"       , "Nd", 144.24    },
	{ "promethium"      , "Pm", 145.      },
	{ "samarium"        , "Sm", 150.4     },
	{ "europium"        , "Eu", 151.96    },
	{ "gadolinium"      , "Gd", 157.25    },
	{ "terbium"         , "Tb", 158.9254  },
	{ "dysprosium"      , "Dy", 162.50    },
	{ "holmium"         , "Ho", 164.9304  },
	{ "erbium"          , "Er", 167.26    },
	{ "thulium"         , "Tm", 168.9342  },
	{ "ytterbium"       , "Yb", 173.04    },
	{ "lutetium"        , "Lu", 174.967   },
	{ "hafnium"         , "Hf", 178.49    },
	{ "tantalum"        , "Ta", 180.9479  },
	{ "tungsten"        , "W" , 183.85    },
	{ "rhenium"         , "Re", 186.207   },
	{ "osmium"          , "Os", 190.2     },
	{ "iridium"         , "Ir", 192.22    },
	{ "platinum"        , "Pt", 195.09    },
	{ "gold"            , "Au", 196.9665  },
	{ "mercury"         , "Hg", 200.59    },
	{ "thallium"        , "Tl", 204.37    },
	{ "lead"            , "Pb", 207.2     },
	{ "bismuth"         , "Bi", 208.9804  },
	{ "polonium"        , "Po", 209.      },
	{ "astatine"        , "At", 210.      },
	{ "radon"           , "Rn", 222.      },
	{ "francium"        , "Fr", 223.      },
	{ "radium"          , "Ra", 226.0254  },
	{ "actinium"        , "Ac", 227.      },
	{ "thorium"         , "Th", 232.0381  },
	{ "protactinium"    , "Pa", 231.0359  },
	{ "uranium"         , "U" , 238.029   },
	{ "neptunium"       , "Np", 237.0482  },
	{ "plutonium"       , "Pu", 244.      },
	{ "americium"       , "Am", 243.      },
	{ "curium"          , "Cm", 247.      },
	{ "berkelium"       , "Bk", 247.      },
	{ "californium"     , "Cf", 251.      },
	{ "einsteinium"     , "Es", 254.      },
	{ "fermium"         , "Fm", 257.      },
	{ "mendelevium"     , "Md", 258.      },
	{ "nobelium"        , "No", 259.      },
	{ "lawrencium"      , "Lw", 260.      },
	{ "ruferfordium"    , "Rf", 260.      },
	{ "hahnium"         , "Ha", 260.      },
	{ "106ium"          , "ZZ", 263.      }
};
static const NSUInteger ne = (sizeof elements / sizeof (struct _elementsData));

typedef struct _entropyData {
    const char *name;
    const char *sym;
    double     entropy;
} EntropyData;

// These entropy values are from Robie, Hemingway and Fisher (1979) USGS Bull 1452
//    as stipulated by Berman (1988).  They are NOT the most recent values (e.g.NIST)
static const EntropyData RobieEtAl1979[] = {
    { "dummy"           , "!" ,     0.0      },
    { "hydrogen"        , "H" ,   130.68/2.0 },
    { "helium"          , "He",   126.15     },
    { "lithium"         , "Li",    29.12     },
    { "beryllium"       , "Be",     9.54     },
    { "boron"           , "B" ,     5.90     },
    { "carbon"          , "C" ,     5.74     },
    { "nitrogen"        , "N" ,   191.61/2.0 },
    { "oxygen"          , "O" ,   205.15/2.0 },
    { "fluorine"        , "F" ,   202.79/2.0 },
    { "neon"            , "Ne",   146.32     },
    { "sodium"          , "Na",    51.30     },
    { "magnesium"       , "Mg",    32.68     },
    { "aluminum"        , "Al",    28.35     },
    { "silicon"         , "Si",    18.81     },
    { "phosphorous"     , "P" ,    22.85     },
    { "sulfur"          , "S" ,    31.80     },
    { "chlorine"        , "Cl",   223.08/2.0 },
    { "argon"           , "Ar",   154.84     },
    { "potassium"       , "K" ,    64.68     },
    { "calcium"         , "Ca",    41.63     },
    { "scandium"        , "Sc",    34.64     },
    { "titanium"        , "Ti",    30.63     },
    { "vanadium"        , "V" ,    28.91     },
    { "chromium"        , "Cr",    23.64     },
    { "manganese"       , "Mn",    32.01     },
    { "iron"            , "Fe",    27.28     },
    { "cobalt"          , "Co",    30.04     },
    { "nickel"          , "Ni",    29.87     },
    { "copper"          , "Cu",    33.15     },
    { "zinc"            , "Zn",    41.63     },
    { "gallium"         , "Ga",    40.83     },
    { "germanium"       , "Ge",    31.09     },
    { "arsenic"         , "As",    35.69     },
    { "selenium"        , "Se",    42.27     },
    { "bromine"         , "Br",   245.46/2.0 },
    { "krypton"         , "Kr",   164.08     },
    { "rubidium"        , "Rb",    76.78     },
    { "strontium"       , "Sr",    55.40     },
    { "yttrium"         , "Y" ,    44.43     },
    { "zirconium"       , "Zr",    38.99     },
    { "niobium"         , "Nb",    36.40     },
    { "molybdenum"      , "Mo",    28.66     },
    { "technetium"      , "Tc",      DBL_MAX },
    { "ruthenium"       , "Ru",    28.53     },
    { "rhodium"         , "Rh",    31.54     },
    { "palladium"       , "Pd",    37.82     },
    { "silver"          , "Ag",    42.55     },
    { "cadmium"         , "Cd",    51.80     },
    { "indium"          , "In",    57.84     },
    { "tin"             , "Sn",    51.20     },
    { "antimony"        , "Sb",    45.52     },
    { "tellurium"       , "Te",    49.50     },
    { "iodine"          , "I" ,   116.15/2.0 },
    { "xenon"           , "Xe",   169.68     },
    { "cesium"          , "Cs",    85.23     },
    { "barium"          , "Ba",    62.42     },
    { "lantahnum"       , "La",    56.90     },
    { "cerium"          , "Ce",    69.46     },
    { "praseodymium"    , "Pr",    73.93     },
    { "neodymium"       , "Nd",    71.09     },
    { "promethium"      , "Pm",      DBL_MAX },
    { "samarium"        , "Sm",    69.50     },
    { "europium"        , "Eu",    80.79     },
    { "gadolinium"      , "Gd",    68.45     },
    { "terbium"         , "Tb",    73.30     },
    { "dysprosium"      , "Dy",    74.89     },
    { "holmium"         , "Ho",    75.02     },
    { "erbium"          , "Er",    73.18     },
    { "thulium"         , "Tm",    74.01     },
    { "ytterbium"       , "Yb",    59.83     },
    { "lutetium"        , "Lu",    50.96     },
    { "hafnium"         , "Hf",    43.56     },
    { "tantalum"        , "Ta",    41.51     },
    { "tungsten"        , "W" ,    32.64     },
    { "rhenium"         , "Re",    36.53     },
    { "osmium"          , "Os",    32.64     },
    { "iridium"         , "Ir",    35.48     },
    { "platinum"        , "Pt",    41.63     },
    { "gold"            , "Au",    47.49     },
    { "mercury"         , "Hg",    75.90     },
    { "thallium"        , "Tl",    64.18     },
    { "lead"            , "Pb",    65.06     },
    { "bismuth"         , "Bi",    56.74     },
    { "polonium"        , "Po",      DBL_MAX },
    { "astatine"        , "At",      DBL_MAX },
    { "radon"           , "Rn",   176.23     },
    { "francium"        , "Fr",      DBL_MAX },
    { "radium"          , "Ra",      DBL_MAX },
    { "actinium"        , "Ac",      DBL_MAX },
    { "thorium"         , "Th",    53.39     },
    { "protactinium"    , "Pa",      DBL_MAX },
    { "uranium"         , "U" ,    50.29     },
    { "neptunium"       , "Np",      DBL_MAX },
    { "plutonium"       , "Pu",    51.46     },
    { "americium"       , "Am",      DBL_MAX },
    { "curium"          , "Cm",      DBL_MAX },
    { "berkelium"       , "Bk",      DBL_MAX },
    { "californium"     , "Cf",      DBL_MAX },
    { "einsteinium"     , "Es",      DBL_MAX },
    { "fermium"         , "Fm",      DBL_MAX },
    { "mendelevium"     , "Md",      DBL_MAX },
    { "nobelium"        , "No",      DBL_MAX },
    { "lawrencium"      , "Lw",      DBL_MAX },
    { "ruferfordium"    , "Rf",      DBL_MAX },
    { "hahnium"         , "Ha",      DBL_MAX },
    { "106ium"          , "ZZ",      DBL_MAX }
};

@synthesize phaseFormula, phaseName, mw, formulaAsElementArray, entropyFromRobieEtAl1979;

+(NSString *)elementNameFromAtomicNumber:(NSUInteger)atomicNumber {
	if      (atomicNumber ==   0) return nil;
	else if (atomicNumber ==   1) return @"H";
	else if (atomicNumber ==   2) return @"He";
	else if (atomicNumber ==   3) return @"Li";
	else if (atomicNumber ==   4) return @"Be";
	else if (atomicNumber ==   5) return @"B";
	else if (atomicNumber ==   6) return @"C";
	else if (atomicNumber ==   7) return @"N";
	else if (atomicNumber ==   8) return @"O";
	else if (atomicNumber ==   9) return @"F";
	else if (atomicNumber ==  10) return @"Ne";
	else if (atomicNumber ==  11) return @"Na";
	else if (atomicNumber ==  12) return @"Mg";
	else if (atomicNumber ==  13) return @"Al";
	else if (atomicNumber ==  14) return @"Si";
	else if (atomicNumber ==  15) return @"P";
	else if (atomicNumber ==  16) return @"S";
	else if (atomicNumber ==  17) return @"Cl";
	else if (atomicNumber ==  18) return @"Ar";
	else if (atomicNumber ==  19) return @"K";
	else if (atomicNumber ==  20) return @"Ca";
	else if (atomicNumber ==  21) return @"Sc";
	else if (atomicNumber ==  22) return @"Ti";
	else if (atomicNumber ==  23) return @"V";
	else if (atomicNumber ==  24) return @"Cr";
	else if (atomicNumber ==  25) return @"Mn";
	else if (atomicNumber ==  26) return @"Fe";
	else if (atomicNumber ==  27) return @"Co";
	else if (atomicNumber ==  28) return @"Ni";
	else if (atomicNumber ==  29) return @"Cu";
	else if (atomicNumber ==  30) return @"Zn";
	else if (atomicNumber ==  31) return @"Ga";
	else if (atomicNumber ==  32) return @"Ge";
	else if (atomicNumber ==  33) return @"As";
	else if (atomicNumber ==  34) return @"Se";
	else if (atomicNumber ==  35) return @"Br";
	else if (atomicNumber ==  36) return @"Kr";

	else if (atomicNumber ==  37) return @"Rb";
	else if (atomicNumber ==  38) return @"Sr";
	else if (atomicNumber ==  39) return @"Y";
	else if (atomicNumber ==  40) return @"Zr";
	else if (atomicNumber ==  41) return @"Nb";
	else if (atomicNumber ==  42) return @"Mo";
	else if (atomicNumber ==  43) return @"Tc";
	else if (atomicNumber ==  44) return @"Ru";
	else if (atomicNumber ==  45) return @"Rh";
	else if (atomicNumber ==  46) return @"Pd";
	else if (atomicNumber ==  47) return @"Ag";
	else if (atomicNumber ==  48) return @"Cd";
	else if (atomicNumber ==  49) return @"In";
	else if (atomicNumber ==  50) return @"Sn";
	else if (atomicNumber ==  51) return @"Sb";
	else if (atomicNumber ==  52) return @"Te";
	else if (atomicNumber ==  53) return @"I";
	else if (atomicNumber ==  54) return @"Xe";
	else if (atomicNumber ==  55) return @"Cs";
	else if (atomicNumber ==  56) return @"Ba";
	else if (atomicNumber ==  57) return @"La";
	else if (atomicNumber ==  58) return @"Ce";
	else if (atomicNumber ==  59) return @"Pr";
	else if (atomicNumber ==  60) return @"Nd";
	else if (atomicNumber ==  61) return @"Pm";
	else if (atomicNumber ==  62) return @"Sm";
	else if (atomicNumber ==  63) return @"Eu";
	else if (atomicNumber ==  64) return @"Gd";
	else if (atomicNumber ==  65) return @"Tb";
	else if (atomicNumber ==  66) return @"Dy";
	else if (atomicNumber ==  67) return @"Ho";
	else if (atomicNumber ==  68) return @"Er";
	else if (atomicNumber ==  69) return @"Tm";
	else if (atomicNumber ==  70) return @"Yb";
	else if (atomicNumber ==  71) return @"Lu";
	else if (atomicNumber ==  72) return @"Hf";
	else if (atomicNumber ==  73) return @"Ta";
	else if (atomicNumber ==  74) return @"W";
	else if (atomicNumber ==  75) return @"Re";
	else if (atomicNumber ==  76) return @"Os";
	else if (atomicNumber ==  77) return @"Ir";
	else if (atomicNumber ==  78) return @"Pt";
	else if (atomicNumber ==  79) return @"Au";
	else if (atomicNumber ==  80) return @"Hg";
	else if (atomicNumber ==  81) return @"Tl";
	else if (atomicNumber ==  82) return @"Pb";
	else if (atomicNumber ==  83) return @"Bi";
	else if (atomicNumber ==  84) return @"Po";
	else if (atomicNumber ==  85) return @"At";
	else if (atomicNumber ==  86) return @"Rn";
	else if (atomicNumber ==  87) return @"Fr";
	else if (atomicNumber ==  88) return @"Ra";
	else if (atomicNumber ==  89) return @"Ac";
	else if (atomicNumber ==  90) return @"Th";
	else if (atomicNumber ==  91) return @"Pa";
	else if (atomicNumber ==  92) return @"U";
	else if (atomicNumber ==  93) return @"Np";
	else if (atomicNumber ==  94) return @"Pu";
	else if (atomicNumber ==  95) return @"Am";
	else if (atomicNumber ==  96) return @"Cm";
	else if (atomicNumber ==  97) return @"Bk";
	else if (atomicNumber ==  98) return @"Cf";
	else if (atomicNumber ==  99) return @"Es";
	else if (atomicNumber == 100) return @"Fm";
	else if (atomicNumber == 101) return @"Md";
	else if (atomicNumber == 102) return @"No";
	else if (atomicNumber == 103) return @"Lr";
	else if (atomicNumber == 104) return @"Rf";
	else if (atomicNumber == 105) return @"Db";
	else if (atomicNumber == 106) return @"Sg";
	else                         return nil;
}

-(id)init {
	if ((self = [super init])) {
		phaseName = [NSString string];
		phaseFormula = [NSString string];
		mw = 0.0;
        formulaAsElementArray = [[DoubleVector alloc]initWithSize:ne andInitialValue:0.0];
	}
	return self;
}

-(double)convertElementsToMolesOfPhase:(double *)e {
	double sumE = 0.0, sumRef = 0.0;
	double *ref = [[self formulaAsElementArray] pointerToDouble];
	for (NSUInteger i=1; i<107; i++) {
		sumE += e[i];
		sumRef += ref[i];
	}
	return (sumRef != 0.0) ? sumE/sumRef : 0.0;
}

-(double)convertElementsToMassOfPhase:(double *)e {
	return [self mw]*[self convertElementsToMolesOfPhase:e];
}

/**
 "C" helper function to parse formula string and calculate a conversion vector and a molecular weight
 */
static double formulaToMwStoich(const char *formula, double *stoich)
{
	NSUInteger len, i, j;
	char c, *sym, *num, *temp;
	double mw = 0.0, mult = 1.0;

	if((len = strlen(formula)) == 0) return 0.0;

	sym  = (char *) malloc((unsigned)       3*sizeof(char));
	temp = (char *) malloc((unsigned) (len+2)*sizeof(char));
	num  = (char *) malloc((unsigned) (len+2)*sizeof(char));
	sym[0] = '\0';
	num[0] = '\0';
	j = 0;

	temp = strcpy(temp, formula);
	temp = strcat(temp," ");

	for (i=0; i<=len; i++) {
		c = temp[i];
		if(isalpha(c) || isspace(c) || c == '(' || c == ')') {
			if( !islower(c) && strlen(sym) > 0) {
				for (j=1; j<ne; j++) {
					if(strcmp(sym, elements[j].sym) == 0) break;
				}
				if (j == ne) {free(temp); free(num); return 0.0;}
				sym[0] = '\0';
			}
			if (j > 0) {
				if (num[0] == '\0') {
					stoich[j] += mult;
					mw += elements[j].atomicWt*mult;
				} else {
					stoich[j] += atof(num)*mult;
					mw += atof(num) * elements[j].atomicWt * mult;
				}
				j = 0;
				num[0] = '\0';
			}
			if (c == '(') {
				NSInteger close, finish, k;

				for(close=i+1; temp[close] != ')' && close < len; close++);
				if(close == len) {free(temp); free(num); free(sym); return 0.0;}
				close++;
				for(finish=close; !isupper(temp[finish]) &&
					(temp[finish] != '(') && (finish < len); finish++);

				if (finish > close) {
					for(k=0; k < (finish-close); k++) {
						num[k] = temp[close+k];
					}
					num[k] = '\0';
					mult = atof(num);
					for(k=0; k <= (len+1-finish); k++) {
						temp[close+k] = temp[finish+k];
					}
					len -= strlen(num);
					num[0] = '\0';
				}

			} else if (c == ')') {
				mult = 1.0;
			} else {
				sym = strncat(sym, &temp[i], 1);
			}
		} else if (isdigit(c) || c == '.' || c == '-') { /* decimal amount of element */
			if(num[0] == '\0') {
				for (j=1; j<ne; j++) {
					if(strcmp(sym, elements[j].sym) == 0) break;
				}
				if (j == ne) {free(temp); free(num); return 0.0;}
				sym[0] = '\0';
			}
			num = strncat(num, &temp[i], 1);
		} else {                                  /* Illegal character         */
			free(temp); free(num); return 0.0;
		}
	}

	free(sym);
	free(temp);
	free(num);

	return mw;
}

-(void)setPhaseFormula:(NSString *)aFormula {
	phaseFormula = [NSString stringWithString:aFormula];
    double *pointerToFormulaAsElementArray = [[self formulaAsElementArray] pointerToDouble];
	[self setMw:formulaToMwStoich([aFormula cStringUsingEncoding:[NSString defaultCStringEncoding]],
                                  pointerToFormulaAsElementArray)];
    double sum = 0.0;
    for (NSUInteger i=1; i<107; i++) if (pointerToFormulaAsElementArray[i] != 0.0) {
        sum += pointerToFormulaAsElementArray[i] * RobieEtAl1979[i].entropy;
    }
    [self setEntropyFromRobieEtAl1979:sum];
}

#pragma mark -
#pragma mark NSSecureCoding protocol

static NSString *kphaseFormula             = @"phaseFormula";
static NSString *kphaseName                = @"phaseName";
static NSString *kformulaAsElementArray    = @"formulaAsElementArray";
static NSString *kmw                       = @"mw";
static NSString *kentropyFromRobieEtAl1979 = @"entropyFromRobieEtAl1979";

- (id)initWithCoder:(NSCoder *)aDecoder {
    if ((self = [super init])) {
#ifdef __APPLE__
        phaseFormula = (NSString *) [aDecoder decodeObjectOfClass:[NSString class] forKey:kphaseFormula];
        phaseName    = (NSString *) [aDecoder decodeObjectOfClass:[NSString class] forKey:kphaseName];
        formulaAsElementArray = (DoubleVector *) [aDecoder decodeObjectOfClass:[DoubleVector class] forKey:kformulaAsElementArray];
#else
        phaseFormula = (NSString *) [aDecoder decodeObjectForKey:kphaseFormula];
        phaseName    = (NSString *) [aDecoder decodeObjectForKey:kphaseName];
        formulaAsElementArray = (DoubleVector *) [aDecoder decodeObjectForKey:kformulaAsElementArray];
#endif
        mw = (double) [aDecoder decodeDoubleForKey:kmw];
        entropyFromRobieEtAl1979 = (double) [aDecoder decodeDoubleForKey:kentropyFromRobieEtAl1979];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)aCoder {
    if ([aCoder isKindOfClass:[NSKeyedArchiver class]]) {
        [aCoder encodeObject:phaseFormula forKey:kphaseFormula];
        [aCoder encodeObject:phaseName forKey:kphaseName];
        [aCoder encodeObject:formulaAsElementArray forKey:kformulaAsElementArray];
        [aCoder encodeDouble:mw forKey:kmw];
        [aCoder encodeDouble:entropyFromRobieEtAl1979 forKey:kentropyFromRobieEtAl1979];
    } else [NSException raise:NSInvalidArchiveOperationException format:@"Class DoubleVector only supports NSKeyedArchiver coders."];
}

+ (BOOL)supportsSecureCoding {
    return YES;
}


@end
