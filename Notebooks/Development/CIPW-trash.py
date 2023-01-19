


    Al2O3 =  0
    BaO   =  1
    CO2   =  2
    CaO   =  3
    Cl2O  =  4
    CoO   =  5
    Cr2O3 =  6
    F2O   =  7
    Fe2O3 =  8
    FeO   =  9
    H2O   = 10
    K2O   = 11
    MgO   = 12
    MnO   = 13
    Na2O  = 14
    NiO   = 15
    P2O5  = 16
    S     = 17
    SO3   = 18
    SiO2  = 19
    SrO   = 20
    TiO2  = 21
    ZrO2  = 22

    n0=23


    Acmite		 =  0
    Albite		 =  1
    Anorthite		 =  2
    Apatite		 =  3
    Calcite		 =  4
    CalciumOrthosilicate  =  5
    Chromite		 =  6
    Corundum		 =  7
    Diopside		 =  8
    Enstatite		 =  9
    Fayalite		 = 10
    Ferrosilite  	 = 11
    Fluorite		 = 12
    Forsterite		 = 13
    Halite		 = 14
    Hedenbergite 	 = 15
    Hematite		 = 16
    Hypersthene  	 = 17
    Ilmenite		 = 18
    Kaliophilite 	 = 19
    Leucite		 = 20
    Magnetite		 = 21
    Nepheline		 = 22
    Olivine		 = 23
    Orthoclase		 = 24
    Perofskite		 = 25
    PotassiumMetasilicate = 26
    Pyrite		 = 27
    Quartz		 = 28
    Rutile		 = 29
    SodiumCarbonate	 = 30
    SodiumMetasilicate	 = 31
    Thenardite		 = 32
    Titanite		 = 33
    Wollastonite 	 = 34
    Zircon		 = 35


    nM                    = 36


    oxideMW = new double[nO];



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



    mineral_names = [
        'Acmite', 'Albite', 'Anorthite', 'Apatite', 'Calcite',
        'CalciumOrthosilicate', 'Chromite', 'Corundum', 'Diopside',
        'Enstatite', 'Fayalite', 'Ferrosilite', 'Fluorite', 'Forsterite',
        'Halite', 'Hedenbergite', 'Hematite', 'Hypersthene', 'Ilmenite',
        'Kaliophilite', 'Leucite', 'Magnetite', 'Nepheline', 'Olivine',
        'Orthoclase', 'Perofskite', 'PotassiumMetasilicate', 'Pyrite',
        'Quartz', 'Rutile', 'SodiumCarbonate', 'SodiumMetasilicate',
        'Thenardite', 'Titanite', 'Wollastonite', 'Zircon' ]


    mineralMW = new double[nM];


    # 'Acmite??': 'Acm'
    # 'CalciumOrthosilicate?': CaOrsi
    # 'Fluorite?': 'Fl'
    # 'Halite?': 'Hl'
    # 'Hypersthene?': 'Hyp'
    # 'Kaliophilite?': 'Klpl'
    # 'PotassiumMetasilicate?': 'KMtsi'
    # 'Pyrite?': 'Py'
    # 'SodiumCarbonate?': 'NaCb'
    # 'SodiumMetasilicate?': 'NaMtsi'
    # 'Thenardite?': 'Thnd'
    # 'Titanite?': 'Ttn'
    # 'Zircon?': 'Zrn'

    # nM  = len(mineral_names)


def CIPW_norm(oxLabels):
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
