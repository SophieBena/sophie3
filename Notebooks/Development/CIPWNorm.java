package client;

import java.text.*;
import java.util.*;

/**CIPW Norm class
 *
 * @author Mark S. Ghiorso
 * @version 1.0.0, 4/27/98 
 */
public class CIPWNorm extends Object {

  private static final int Al2O3 =  0;
  private static final int BaO   =  1;
  private static final int CO2   =  2;
  private static final int CaO   =  3;
  private static final int Cl2O  =  4;
  private static final int CoO   =  5;
  private static final int Cr2O3 =  6;
  private static final int F2O   =  7;
  private static final int Fe2O3 =  8;
  private static final int FeO   =  9;
  private static final int H2O   = 10;
  private static final int K2O   = 11;
  private static final int MgO   = 12;
  private static final int MnO   = 13;
  private static final int Na2O  = 14;
  private static final int NiO   = 15;
  private static final int P2O5  = 16;
  private static final int S     = 17;
  private static final int SO3   = 18;
  private static final int SiO2  = 19;
  private static final int SrO   = 20;
  private static final int TiO2  = 21;
  private static final int ZrO2  = 22;
  
  private static final int nO    = 23;
  private static double oxideMW[];

  private static final int Acmite		 =  0;
  private static final int Albite		 =  1;
  private static final int Anorthite		 =  2;
  private static final int Apatite		 =  3;
  private static final int Calcite		 =  4;
  private static final int CalciumOrthosilicate  =  5;
  private static final int Chromite		 =  6;
  private static final int Corundum		 =  7;
  private static final int Diopside		 =  8;
  private static final int Enstatite		 =  9;
  private static final int Fayalite		 = 10;
  private static final int Ferrosilite  	 = 11;
  private static final int Fluorite		 = 12;
  private static final int Forsterite		 = 13;
  private static final int Halite		 = 14;
  private static final int Hedenbergite 	 = 15;
  private static final int Hematite		 = 16;
  private static final int Hypersthene  	 = 17;
  private static final int Ilmenite		 = 18;
  private static final int Kaliophilite 	 = 19;
  private static final int Leucite		 = 20;
  private static final int Magnetite		 = 21;
  private static final int Nepheline		 = 22;
  private static final int Olivine		 = 23;
  private static final int Orthoclase		 = 24;
  private static final int Perofskite		 = 25;
  private static final int PotassiumMetasilicate = 26;
  private static final int Pyrite		 = 27;
  private static final int Quartz		 = 28;
  private static final int Rutile		 = 29;
  private static final int SodiumCarbonate	 = 30;
  private static final int SodiumMetasilicate	 = 31;
  private static final int Thenardite		 = 32;
  private static final int Titanite		 = 33;
  private static final int Wollastonite 	 = 34;
  private static final int Zircon		 = 35;
  
  public  static final int nM                    = 36;
  public  static String mineralLabel[];

  private static double mineralMW[];
  private static NumberFormat nf;

  static {
    oxideMW = new double[nO];
    oxideMW[Al2O3] = 101.96128;
    oxideMW[BaO  ] = 153.3394;
    oxideMW[CO2  ] =  44.0098;
    oxideMW[CaO  ] =  56.0794;
    oxideMW[Cl2O ] =  86.9054;
    oxideMW[CoO  ] =  74.9326;
    oxideMW[Cr2O3] = 151.9902;
    oxideMW[F2O  ] =  53.9962;
    oxideMW[Fe2O3] = 159.6922;
    oxideMW[FeO  ] =  71.8464;
    oxideMW[H2O  ] =  18.0152;
    oxideMW[K2O  ] =  94.196;
    oxideMW[MgO  ] =  40.3044;
    oxideMW[MnO  ] =  70.9374;
    oxideMW[Na2O ] =  61.97894;
    oxideMW[NiO  ] =  74.6994;
    oxideMW[P2O5 ] = 141.94452;
    oxideMW[S    ] =  32.06;
    oxideMW[SO3  ] =  80.0582;
    oxideMW[SiO2 ] =  60.0843;
    oxideMW[SrO  ] = 103.6194;
    oxideMW[TiO2 ] =  79.8988;
    oxideMW[ZrO2 ] = 123.2188;   
  
    mineralLabel = new String[nM];  
    mineralLabel[Acmite               ] = "ac";
    mineralLabel[Albite               ] = "ab";
    mineralLabel[Anorthite            ] = "an";
    mineralLabel[Apatite              ] = "ap";
    mineralLabel[Calcite              ] = "cc";
    mineralLabel[CalciumOrthosilicate ] = "cs";
    mineralLabel[Chromite             ] = "cm";
    mineralLabel[Corundum             ] = "C ";
    mineralLabel[Diopside             ] = "di";
    mineralLabel[Enstatite            ] = "  ";
    mineralLabel[Fayalite             ] = "  ";
    mineralLabel[Ferrosilite          ] = "  ";
    mineralLabel[Fluorite             ] = "fr";
    mineralLabel[Forsterite           ] = "  ";
    mineralLabel[Halite               ] = "hl";
    mineralLabel[Hedenbergite         ] = "  ";
    mineralLabel[Hematite             ] = "hm";
    mineralLabel[Hypersthene          ] = "hy";
    mineralLabel[Ilmenite             ] = "il";
    mineralLabel[Kaliophilite         ] = "kp";
    mineralLabel[Leucite              ] = "lc";
    mineralLabel[Magnetite            ] = "mt";
    mineralLabel[Nepheline            ] = "ne";
    mineralLabel[Olivine              ] = "ol";
    mineralLabel[Orthoclase           ] = "or";
    mineralLabel[Perofskite           ] = "pf";
    mineralLabel[PotassiumMetasilicate] = "ks";
    mineralLabel[Pyrite               ] = "pr";
    mineralLabel[Quartz               ] = "Q ";
    mineralLabel[Rutile               ] = "ru";
    mineralLabel[SodiumCarbonate      ] = "nc";
    mineralLabel[SodiumMetasilicate   ] = "ns";
    mineralLabel[Thenardite           ] = "th";
    mineralLabel[Titanite             ] = "tn";
    mineralLabel[Wollastonite         ] = "wo";
    mineralLabel[Zircon               ] = "Z ";

    mineralMW = new double[nM];
    mineralMW[Acmite               ] = 231.00417*2.0; // Na2O Fe2O3 4SiO2
    mineralMW[Albite               ] = 262.22301*2.0; // Na2O Al2O3 6SiO2
    mineralMW[Anorthite            ] = 278.20928;     // CaO Al2O3 2SiO2
    mineralMW[Apatite              ] = 930.54816;     // 3(3CaO P2O5)
    mineralMW[Calcite              ] = 100.0892;      // CaO CO2
    mineralMW[CalciumOrthosilicate ] = 172.2431;      // 2CaO SiO2
    mineralMW[Chromite             ] = 223.8366;      // FeO Cr2O3
    mineralMW[Corundum             ] = 101.9613;      // Al2O3
    mineralMW[Diopside             ] = 216.55240;     // CaO MgO 2SiO2
    mineralMW[Enstatite            ] = 100.3887;      // MgO SiO2
    mineralMW[Fayalite             ] = 203.7771;      // 2FeO SiO2
    mineralMW[Ferrosilite          ] = 131.9307;      // FeO SiO2
    mineralMW[Fluorite             ] =  78.076806;    // CaF2
    mineralMW[Forsterite           ] = 140.6931;      // 2MgO SiO4
    mineralMW[Halite               ] =  58.44277;     // NaCl
    mineralMW[Hedenbergite         ] = 248.0944;      // CaO FeO 2SiO2
    mineralMW[Hematite             ] = 159.6922;      // Fe2O3
    mineralMW[Hypersthene          ] =   0.0;         // (Mg, Fe)O SiO3
    mineralMW[Ilmenite             ] = 151.7452;      // FeO TiO2
    mineralMW[Kaliophilite         ] = 158.16294*2.0; // K2O Al2O3 2SiO2
    mineralMW[Leucite              ] = 218.24724*2.0; // K2O Al2O3 4SiO2
    mineralMW[Magnetite            ] = 231.5386;      // FeO Fe2O3
    mineralMW[Nepheline            ] = 142.05441*2.0; // Na2O Al2O3 2SiO2
    mineralMW[Olivine              ] =   0.0;         // 2(Mg, Fe)O SiO2
    mineralMW[Orthoclase           ] = 278.33154*2.0; // K2O Al2O3 6SiO2
    mineralMW[Perofskite           ] = 135.9782;      // CaO TiO2
    mineralMW[PotassiumMetasilicate] = 154.2803;      // K2O SiO2
    mineralMW[Pyrite               ] = 119.967;       // FeS2
    mineralMW[Quartz               ] =  60.0843;      // SiO2
    mineralMW[Rutile               ] =  79.8988;      // TiO2
    mineralMW[SodiumCarbonate      ] = 105.98874;     // Na2O CO2
    mineralMW[SodiumMetasilicate   ] = 122.06324;     // Na2O SiO2
    mineralMW[Thenardite           ] = 142.03714;     // Na2O SO3
    mineralMW[Titanite             ] = 196.0625;      // CaO TiO2 SiO2
    mineralMW[Wollastonite         ] = 116.1637;      // CaO SiO3
    mineralMW[Zircon               ] = 183.3031;      // ZrO2 SiO2

    nf = NumberFormat.getInstance();
    nf.setMaximumFractionDigits(3);
    nf.setMinimumFractionDigits(3);
    nf.setGroupingUsed(false);
  }
  
  private int index[];
  private int nI;

  public CIPWNorm(String[] oxLabels) {
    if ((nI = oxLabels.length) == 0) return;

    index   = new int[nI];    
    for (int i=0; i<nI; i++) {
      index[i] = -1;
      if     (oxLabels[i].equals("Al2O3" )) index[i] = Al2O3;
      else if(oxLabels[i].equals("CO2"   )) index[i] = CO2  ;
      else if(oxLabels[i].equals("CaO"   )) index[i] = CaO  ;
      else if(oxLabels[i].equals("Cl2O-1")) index[i] = Cl2O ;
      else if(oxLabels[i].equals("CoO"   )) index[i] = CoO  ;
      else if(oxLabels[i].equals("Cr2O3" )) index[i] = Cr2O3;
      else if(oxLabels[i].equals("F2O-1" )) index[i] = F2O  ;
      else if(oxLabels[i].equals("Fe2O3" )) index[i] = Fe2O3;
      else if(oxLabels[i].equals("FeO"   )) index[i] = FeO  ;
      else if(oxLabels[i].equals("H2O"   )) index[i] = H2O  ;
      else if(oxLabels[i].equals("K2O"   )) index[i] = K2O  ;
      else if(oxLabels[i].equals("MgO"   )) index[i] = MgO  ;
      else if(oxLabels[i].equals("MnO"   )) index[i] = MnO  ;
      else if(oxLabels[i].equals("Na2O"  )) index[i] = Na2O ;
      else if(oxLabels[i].equals("NiO"   )) index[i] = NiO  ;
      else if(oxLabels[i].equals("P2O5"  )) index[i] = P2O5 ;
      else if(oxLabels[i].equals("S"     )) index[i] = S    ;
      else if(oxLabels[i].equals("SO3"   )) index[i] = SO3  ;
      else if(oxLabels[i].equals("SiO2"  )) index[i] = SiO2 ;
      else if(oxLabels[i].equals("TiO2"  )) index[i] = TiO2 ;
      else if(oxLabels[i].equals("ZrO2"  )) index[i] = ZrO2 ;
    }
  }
  
  public double[] norm(double[] oxWt) {
    double oxMole[]  = new double[nO];
    double minMole[] = new double[nM];
    
    if (nI == 0) return null;
    
    // Step 1
    for (int i=0; i<nO; i++) oxMole[i] = 0.0;
    for (int i=0; i<nI; i++) if (index[i] != -1) oxMole[index[i]] = oxWt[i]/oxideMW[index[i]];
    for (int i=0; i<nM; i++) minMole[i] = 0.0;

    // Step 2
    oxMole[FeO] += oxMole[MnO] + oxMole[NiO];
    oxMole[CaO] += oxMole[BaO] + oxMole[SrO];
    
    // Step 3a (fluorine is ignored in apatite)
    minMole[Apatite]  = oxMole[P2O5]/3.0; 
    oxMole[CaO]      -= 3.33*oxMole[P2O5];
    oxMole[P2O5]      = 0.0;
    
    // Step 3b
    minMole[Halite]  = 2.0*oxMole[Cl2O];
    oxMole[Na2O]    -= oxMole[Cl2O];
    oxMole[Cl2O]     = 0.0;
    
    // Step 3c (SO3 is oxidized sulfur)
    minMole[Thenardite]  = oxMole[SO3];
    oxMole[Na2O]        -= oxMole[SO3];
    oxMole[SO3]          = 0.0;
    
    // Step 3d (S is reduced sulfur)
    minMole[Pyrite]  = oxMole[S]/2.0;
    oxMole[FeO]     -= oxMole[S]/2.0;
    oxMole[S]        = 0.0;
    
    // Step 3e 
    minMole[Chromite]  = oxMole[Cr2O3];
    oxMole[FeO]       -= oxMole[Cr2O3];
    oxMole[Cr2O3]      = 0.0;
    
    // Step 3f
    if(oxMole[FeO] > oxMole[TiO2]) {
      minMole[Ilmenite]  = oxMole[TiO2];
      oxMole[FeO]       -= oxMole[TiO2];
      oxMole[TiO2]       = 0.0;
    } else {
      minMole[Ilmenite]  = oxMole[FeO];
      oxMole[TiO2]      -= oxMole[FeO];
      oxMole[FeO]        = 0.0;
    }
    
    // Step 3g
    minMole[Fluorite]  = oxMole[F2O]/2.0;
    oxMole[CaO]       -= oxMole[F2O]/2.0;
    oxMole[F2O]        = 0.0;
    
    // Step 3h (assume CO2 is 1/2 calcite, 1/2 Cancrinite)
    minMole[SodiumCarbonate]  = oxMole[CO2]/2.0;
    minMole[Calcite]          = oxMole[CO2]/2.0;
    oxMole[Na2O]             -= oxMole[CO2]/2.0;
    oxMole[CaO]              -= oxMole[CO2]/2.0;
    oxMole[CO2]               = 0.0;
    
    // Step 3i
    minMole[Zircon] = oxMole[ZrO2];
    oxMole[ZrO2]    = 0.0;
    
    // Step 4a
    if(oxMole[Al2O3] >= oxMole[K2O]) {
      minMole[Orthoclase]  = oxMole[K2O];
      oxMole[Al2O3]       -= oxMole[K2O];
      oxMole[K2O]          = 0.0;
      // Step 4c
      if(oxMole[Al2O3] >= oxMole[Na2O]) {
  	minMole[Albite]  = oxMole[Na2O];
  	oxMole[Al2O3]   -= oxMole[Na2O];
  	oxMole[Na2O]     = 0.0;
  	// Step 4d/f
  	if(oxMole[Al2O3] <= oxMole[CaO]) {
  	  minMole[Anorthite]  = oxMole[Al2O3];
  	  oxMole[CaO]        -= oxMole[Al2O3];
  	  oxMole[Al2O3]       = 0.0;
  	// Step 4e
  	} else {
  	  minMole[Anorthite]  = oxMole[CaO];
  	  oxMole[Al2O3]      -= oxMole[CaO];
  	  oxMole[CaO]         = 0.0;
  	  minMole[Corundum]   = oxMole[Al2O3];
  	  oxMole[Al2O3]       = 0.0;
  	}
      // Step 4g
      } else {
  	minMole[Albite]  = oxMole[Al2O3];
  	oxMole[Na2O]    -= oxMole[Al2O3];
  	oxMole[Al2O3]    = 0.0;
      }
    // Step 4b
    } else {
      minMole[Orthoclase]            = oxMole[Al2O3];
      minMole[PotassiumMetasilicate] = oxMole[K2O] - oxMole[Al2O3];
      oxMole[K2O]                    = 0.0;
      oxMole[Al2O3]                  = 0.0;
    }
    
    if(oxMole[Na2O] > 0.0) {
      // Step 5a
      if(oxMole[Na2O] <= oxMole[Fe2O3]) {
  	minMole[Acmite]  = oxMole[Na2O];
  	oxMole[Fe2O3]   -= oxMole[Na2O];
  	oxMole[Na2O]     = 0.0;
      // Step 5b
      } else {
  	minMole[Acmite]             = oxMole[Fe2O3];
  	minMole[SodiumMetasilicate] = oxMole[Na2O] - oxMole[Fe2O3];
  	oxMole[Na2O]                = 0.0;
  	oxMole[Fe2O3]               = 0.0;
      }
    }
    
    if(oxMole[Fe2O3] > 0.0) {
      // Step 5c
      if(oxMole[FeO] >= oxMole[Fe2O3]) {
  	minMole[Magnetite]  = oxMole[Fe2O3];
  	oxMole[FeO]        -= oxMole[Fe2O3];
  	oxMole[Fe2O3]       = 0.0;
      // Step 5d
      } else {
  	minMole[Magnetite] = oxMole[FeO];
  	minMole[Hematite]  = oxMole[Fe2O3] - oxMole[FeO];
  	oxMole[FeO]        = 0.0;
  	oxMole[Fe2O3]      = 0.0;
      }
    }

    // Step 6
    double ratio = ((oxMole[MgO]+oxMole[FeO])) > 0.0 ? oxMole[MgO]/(oxMole[MgO]+oxMole[FeO]) : 0.0;

    // Step 3f (continued)
    if(oxMole[TiO2] > 0.0) {
      if(oxMole[CaO] >= oxMole[TiO2]) {
  	minMole[Titanite]  = oxMole[TiO2];
  	oxMole[CaO]       -= oxMole[TiO2];
  	oxMole[TiO2]       = 0.0;
      } else {
  	minMole[Titanite] = oxMole[CaO];
  	minMole[Rutile]   = oxMole[TiO2] - oxMole[CaO];
  	oxMole[CaO]       = 0.0;
  	oxMole[TiO2]      = 0.0;
      }
    }
    
    // Step 7a/b
    if((oxMole[MgO]+oxMole[FeO]) < oxMole[CaO]) {
      minMole[Diopside]     = oxMole[MgO] + oxMole[FeO];
      minMole[Wollastonite] = oxMole[CaO] - oxMole[MgO] - oxMole[FeO];
      oxMole[MgO]           = 0.0;
      oxMole[FeO]           = 0.0;
      oxMole[CaO]           = 0.0;
    // Step 7c
    } else {
      minMole[Diopside]    = oxMole[CaO];
      minMole[Hypersthene] = oxMole[MgO] + oxMole[FeO] - oxMole[CaO];
      oxMole[MgO]          = 0.0;
      oxMole[FeO]          = 0.0;
      oxMole[CaO]          = 0.0;
    }

    // Step 8a
    double sum = minMole[Titanite] + 4.0*minMole[Acmite] + minMole[PotassiumMetasilicate] 
      + minMole[SodiumMetasilicate] + 6.0*minMole[Orthoclase] + 6.0*minMole[Albite] 
      + 2.0*minMole[Anorthite] + 2.0*minMole[Diopside] + minMole[Wollastonite]
      + minMole[Hypersthene] + minMole[Zircon];
    oxMole[SiO2] = oxMole[SiO2] - sum;

    // Step 8b
    if(oxMole[SiO2] >= 0.0) {
      minMole[Quartz] = oxMole[SiO2];
    // Step 8c
    } else {
      oxMole[SiO2] = oxMole[SiO2] + minMole[Hypersthene];
      if(oxMole[SiO2] >= 0.5*minMole[Hypersthene]) {
  	sum                  = minMole[Hypersthene];
  	minMole[Hypersthene] = 2.0*oxMole[SiO2] - sum;
  	minMole[Olivine]     = sum - oxMole[SiO2];
  	oxMole[SiO2]         = 0.0;
      // Step 8d/e
      } else {
  	minMole[Olivine] = 0.5*minMole[Hypersthene];
  	oxMole[SiO2]         = oxMole[SiO2] - minMole[Olivine] + minMole[Titanite] + 6.0*minMole[Albite];
  	minMole[Hypersthene] = 0.0;
  	minMole[Perofskite]  = minMole[Titanite];
  	minMole[Titanite]    = 0.0;
  	if(oxMole[SiO2] >= 2.0*minMole[Albite]) {
  	  oxMole[Na2O]       = minMole[Albite];
  	  minMole[Albite]    = (oxMole[SiO2]-2.0*oxMole[Na2O])/4.0;
  	  minMole[Nepheline] = oxMole[Na2O] - minMole[Albite];
  	  oxMole[SiO2]       = 0.0;
  	// Step 8f
  	} else {
  	  minMole[Nepheline] = minMole[Albite];
  	  minMole[Albite]    = 0.0;
  	  oxMole[SiO2] = oxMole[SiO2] - 2.0*minMole[Nepheline] + 6.0*minMole[Orthoclase];
  	  if(oxMole[SiO2] >= 4.0*minMole[Orthoclase]) {
  	    oxMole[K2O]         = minMole[Orthoclase];
  	    minMole[Orthoclase] = (oxMole[SiO2]-4.0*oxMole[K2O])/2.0;
  	    minMole[Leucite]    = oxMole[K2O] - minMole[Orthoclase];
  	    oxMole[SiO2]        = 0.0;
  	  // Step 8g
  	  } else {
  	    if((oxMole[SiO2]+minMole[Wollastonite]) < 4.0*minMole[Orthoclase]) {
  	      minMole[Leucite] = minMole[Orthoclase];
  	      minMole[Orthoclase] = 0.0;
  	      oxMole[SiO2] = oxMole[SiO2] - 4.0*minMole[Leucite] + minMole[Wollastonite] + 2.0*minMole[Diopside];
  	      // Step 8h
  	      if((oxMole[SiO2]+minMole[Wollastonite]+minMole[Olivine]) < (minMole[Diopside]+minMole[Olivine])) {
  		oxMole[SiO2] = oxMole[SiO2] - minMole[Diopside] - minMole[Wollastonite]*0.5 + minMole[Leucite]*4.0;
  		minMole[CalciumOrthosilicate] = 0.5*(minMole[Diopside]+minMole[Wollastonite]);
  		minMole[Olivine]      = minMole[Olivine] + 0.5*minMole[Diopside];
  		minMole[Wollastonite] = 0.0;
  		minMole[Diopside]     = 0.0;
  		oxMole[K2O]           = minMole[Leucite];
  		minMole[Leucite]      = (oxMole[SiO2]-2.0*oxMole[K2O])/2.0;
  		minMole[Kaliophilite] = (4.0*oxMole[K2O]-oxMole[SiO2])/2.0;
  		oxMole[SiO2]          = 0.0;
  	      } else {
  		double x = (2.0*oxMole[SiO2]-minMole[Diopside]*2.0-minMole[Wollastonite])/2.0;
  		double y = (minMole[Diopside]-x)/2.0;
  		double z = (minMole[Diopside]+minMole[Wollastonite]-x)/2.0;
  		minMole[Diopside]              = x;
  		minMole[Olivine]              += y;
  		minMole[CalciumOrthosilicate]  = z;
  		minMole[Wollastonite]          = 0.0;
  		oxMole[SiO2]                   = 0.0;
  	      }
  	    } else {
  	      minMole[Leucite]               = minMole[Orthoclase];
  	      minMole[Orthoclase]            = 0.0;
  	      oxMole[SiO2]                   = oxMole[SiO2] - 4.0*minMole[Leucite];
  	      minMole[CalciumOrthosilicate]  = -oxMole[SiO2];
  	      minMole[Wollastonite]         -= 2.0*minMole[CalciumOrthosilicate];
  	      oxMole[SiO2]                   = 0.0;
  	    }
  	  }
  	}
      }
    }

    minMole[Hedenbergite] = (1.0-ratio)*minMole[Diopside];
    minMole[Diopside]     = ratio*minMole[Diopside];
    minMole[Enstatite]    = ratio*minMole[Hypersthene];
    minMole[Ferrosilite]  = (1.0-ratio)*minMole[Hypersthene];
    minMole[Forsterite]   = ratio*minMole[Olivine];
    minMole[Fayalite]     = (1.0-ratio)*minMole[Olivine];
    
    for (int i=0; i<nM; i++) minMole[i] *= mineralMW[i];
    minMole[Diopside]     = minMole[Diopside] + minMole[Hedenbergite];
    minMole[Hypersthene]  = minMole[Enstatite] + minMole[Ferrosilite];
    minMole[Olivine]      = minMole[Forsterite] + minMole[Fayalite];
    minMole[Hedenbergite] = 0.0;
    minMole[Enstatite]    = 0.0;
    minMole[Ferrosilite]  = 0.0;
    minMole[Forsterite]   = 0.0;
    minMole[Fayalite]     = 0.0;
        
    return minMole;
  }
}
