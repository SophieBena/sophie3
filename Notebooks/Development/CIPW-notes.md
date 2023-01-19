
    'Acm'    = 'Na2O Fe2O3 4SiO2'
    'Ab'     = 'Na2O Al2O3 6SiO2'
    'An'     = 'CaO Al2O3 2SiO2'
    'Ap'     = '3(3CaO P2O5)'
    'Cal'    = 'CaO CO2'
    'CaOrsi' = '2CaO SiO2'
    'Chr'    = 'FeO Cr2O3'
    'Crn'    = 'Al2O3'
    'Di'     = 'CaO MgO 2SiO2'
    'En'     = 'MgO SiO2'
    'Fa'     = '2FeO SiO2'
    'Fs'     = 'FeO SiO2'
    'Fl'     = 'CaF2'
    'Fo'     = '2MgO SiO2'
    'Hl'     = 'NaCl'
    'Hd'     = 'CaO FeO 2SiO2'
    'Hem'    = 'Fe2O3'
    'Hyp'    = '(Mg, Fe)O SiO3'
    'Ilm'    = 'FeO TiO2'
    'Klpl'   = 'K2O Al2O3 2SiO2'
    'Lct'    = 'K2O Al2O3 4SiO2'
    'Mag'    = 'FeO Fe2O3'
    'Nph'    = 'Na2O Al2O3 2SiO2'
    'Ol'     = '2(Mg, Fe)O SiO2'
    'Or'     = 'K2O Al2O3 6SiO2'
    'Prv'    = 'CaO TiO2'
    'KMtsi'  = 'K2O SiO2'
    'Py'     = 'FeS2'
    'Qz'     = 'SiO2'
    'Rt'     = 'TiO2'
    'NaCb'   = 'Na2O CO2'
    'NaMtsi' = 'Na2O SiO2'
    'Thnd'   = 'Na2O SO3'
    'Ttn'    = 'CaO TiO2 SiO2'
    'Wo'     = 'CaO SiO3'
    'Zrn'    = 'ZrO2 SiO2'



# Group like oxides
* FeO = FeO + MnO + NiO
* CaO = CaO + BaO + SrO

# mineral order
* Ap: 3(3CaO P2O5)
* Hl: NaCl
* Thnd: Na2O SO3
* Py: FeS2
* Chr: FeO Cr2O3
* Ilm: FeO TiO2
* Fl: CaF2
* 1/2 NaCb + 1/2 Cal: Na2O CO2 + CaO CO2
* Zrn: ZrO2 SiO2
* Or: K2O Al2O3 6SiO2
* KMtsi: K2O SiO2
* Ab: Na2O Al2O3 6SiO2
* An: CaO Al2O3 2SiO2
* Crn: Al2O3
* Acm: Na2O Fe2O3 4SiO2
* NaMtsi: Na2O SiO2
* Mag: FeO Fe2O3
* Hem: Fe2O3

* **remaining Mg_num = Mg/(Mg+Fe)**
* Ttn: CaO TiO2 SiO2
* Rt: TiO2
* Di: CaO MgO 2SiO2
      * Hd: CaO FeO 2SiO2
* Wo: CaO SiO3
* Hyp: (Mg, Fe)O SiO3
      * En: MgO SiO2
      * Fs: FeO SiO2
* Qz: SiO2

* **If Qz present, we're done**

* **If SiO2 undersaturated:**
* -Hyp: (Mg, Fe)O SiO3
      * +Hyp: (Mg, Fe)O SiO3
      * +Ol :2(Mg, Fe)O SiO2

* -Ttn -Hyp +Ol + Prv
      * Ol = 1/2 Hyp, Hyp=0, Prv = Ttn, Ttn=0

* -Ab + Nph
* -Or + Lct + Or?
* -Di -Wo + CaOrsi + Ol
* -Lct + Klpl + Lct




-------------------
# Non-negative scheme
## SiO2-free steps
* +Ap: 3(3CaO P2O5)
* +Hl: NaCl
* +Thnd: Na2O SO3
* +Py: FeS2
* +Chr: FeO Cr2O3
* +Ilm: FeO TiO2
* +Fl: CaF2
* +1/2 NaCb +1/2 Cal: Na2O CO2 + CaO CO2

## SiO2-full steps
* +Zrn: ZrO2 SiO2


* +Or: K2O Al2O3 6SiO2
*     +Lct: K2O Al2O3 4SiO2
* +KMtsi: K2O SiO2

* +Ab: Na2O Al2O3 6SiO2
*     +Nph: Na2O Al2O3 2SiO2
* +An: CaO Al2O3 2SiO2
* +Crn: Al2O3
* +Acm: Na2O Fe2O3 4SiO2
* +NaMtsi: Na2O SiO2
* +Mag: FeO Fe2O3
* +Hem: Fe2O3


* **remaining Mg_num = Mg/(Mg+Fe)**
* +Ttn: CaO TiO2 SiO2
*      +Prv: CaO TiO2
* +Rt: TiO2
* +Di: CaO MgO 2SiO2
      * +Hd: CaO FeO 2SiO2
* +Wo: CaO SiO3
* +Hyp: (Mg, Fe)O SiO3
      * +En: MgO SiO2
      * +Fs: FeO SiO2
* +Qz: SiO2

* **step 8c/d**
* +Ol-Hyp: (Mg, Fe)O
      * +Fo-En: MgO
            * 2MgO SiO2 -(MgO SiO2)
      * +Fa-Fs: FeO
            * 2FeO SiO2 -(FeO SiO2)
* **step 8d**
* +Prv-Ttn: -SiO2
      * CaO TiO2-(CaO TiO2 SiO2)

* **step 8e**
* +Ab-Nph:  +4SiO2
      * Na2O Al2O3 6SiO2 - (Na2O Al2O3 2SiO2)

* **step 8f**
* +Nph-Ab:  -4SiO2
      * Na2O Al2O3 2SiO2 - (Na2O Al2O3 6SiO2)
* +Or-Lct: 2SiO2
      * K2O Al2O3 6SiO2 - (K2O Al2O3 4SiO2)

* **step 8g**
* +Lct-Or: -2SiO2
      * K2O Al2O3 4SiO2 - (K2O Al2O3 6SiO2)

* **step 8h**

K2O = Lct
Lct  = (SiO2- 2*K2O)/2
Klpl = (4*K2O -SiO2)/2

Klpl: K2O Al2O3 2SiO2


# C

* +Klpl-Lct: -2SiO2
      * K2O Al2O3 2SiO2 - (K2O Al2O3 4SiO2)
* use all Di and Wo
* +CaOrsi+Ol-Di-Wo: CaO
      * 2CaO SiO2 + 2(Mg,Fe)O SiO2 -(CaO (Mg,Fe)O 2SiO2 + CaO SiO3)
      * -SiO2 + (Mg,Fe)O
    CaOrsi = 2CaO SiO2
    Ol = 2(Mg, Fe)O SiO2
    Di = CaO MgO 2SiO2
    Wo = CaO SiO3


* +CaOrsi-Wo: CaO
      * 2CaO SiO2 -(CaO SiO2)


* +CaOrsi+Ol-Di: CaO + (Mg, Fe)O
      * 2CaO SiO2 + 2(Mg, Fe)O SiO2 - (CaO (Mg, Fe)O 2SiO2) = CaO + (Mg, Fe)O
      * +CaOrsi+Fo-En: CaO + MgO
      * +CaOrsi+Fa-Fs: CaO + FeO







## K-Al2O3 step
* **Al2O3 > K2O**
* Or: K2O Al2O3 6SiO2
* Ab: Na2O Al2O3 6SiO2
      * An: CaO Al2O3 2SiO2
      * Crn: Al2O3

* **Al2O3 < K2O**
* Or: K2O Al2O3 6SiO2
* KMtsi: K2O SiO2
