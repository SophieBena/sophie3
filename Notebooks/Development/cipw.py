
# CIPW Norm class
#
# @authors Mark S. Ghiorso, Aaron S. Wolf
# @version 1.0.0, 4/12/19 (4/27/98)

import pandas as pd
import numpy as np


def get_oxide_info():
    oxides = [
        'Al2O3', 'BaO', 'CO2', 'CaO', 'Cl2O', 'CoO', 'Cr2O3', 'F2O',
        'Fe2O3','FeO','H2O','K2O','MgO','MnO','Na2O','NiO','P2O5','S',
        'SO3','SiO2', 'SrO','TiO2','ZrO2']

    # Missing: 'CoO', 'H2O',

    oxide_molwts = pd.Series()
    oxide_molwts['Al2O3'] = 101.96128
    oxide_molwts['BaO'  ] = 153.3394
    oxide_molwts['CO2'  ] =  44.0098
    oxide_molwts['CaO'  ] =  56.0794
    oxide_molwts['Cl2O' ] =  86.9054
    oxide_molwts['CoO'  ] =  74.9326
    oxide_molwts['Cr2O3'] = 151.9902
    oxide_molwts['F2O'  ] =  53.9962
    oxide_molwts['Fe2O3'] = 159.6922
    oxide_molwts['FeO'  ] =  71.8464
    oxide_molwts['H2O'  ] =  18.0152
    oxide_molwts['K2O'  ] =  94.196
    oxide_molwts['MgO'  ] =  40.3044
    oxide_molwts['MnO'  ] =  70.9374
    oxide_molwts['Na2O' ] =  61.97894
    oxide_molwts['NiO'  ] =  74.6994
    oxide_molwts['P2O5' ] = 141.94452
    oxide_molwts['S'    ] =  32.06
    oxide_molwts['SO3'  ] =  80.0582
    oxide_molwts['SiO2' ] =  60.0843
    oxide_molwts['SrO'  ] = 103.6194
    oxide_molwts['TiO2' ] =  79.8988
    oxide_molwts['ZrO2' ] = 123.2188

    oxide_info = {}
    oxide_info['symbols'] = oxides
    oxide_info['molwts'] = oxide_molwts
    oxide_info['num'] = len(oxides)

    return oxide_info

def get_mineral_info():
    mineral_symbols = [
        'Acm', 'Ab', 'An', 'Ap', 'Cal',
        'CaOrsi', 'Chr', 'Crn', 'Di',
        'En', 'Fa', 'Fs', 'Fl', 'Fo',
        'Hl', 'Hd', 'Hem', 'Hyp', 'Ilm',
        'Klpl', 'Lct', 'Mag', 'Nph', 'Ol',
        'Or', 'Prv', 'KMtsi', 'Py',
        'Qz', 'Rt', 'NaCb', 'NaMtsi',
        'Thnd', 'Ttn', 'Wo', 'Zrn' ]

    mineral_names = pd.Series()
    mineral_names['Acm']    = 'Acmite'
    mineral_names['Ab']     = 'Albite'
    mineral_names['An']     = 'Anorthite'
    mineral_names['Ap']     = 'Apatite'
    mineral_names['Cal']    = 'Calcite'
    mineral_names['CaOrsi'] = 'CalciumOrthosilicate'
    mineral_names['Chr']    = 'Chromite'
    mineral_names['Crn']    = 'Corundum'
    mineral_names['Di']     = 'Diopside'
    mineral_names['En']     = 'Enstatite'
    mineral_names['Fa']     = 'Fayalite'
    mineral_names['Fs']     = 'Ferrosilite'
    mineral_names['Fl']     = 'Fluorite'
    mineral_names['Fo']     = 'Forsterite'
    mineral_names['Hl']     = 'Halite'
    mineral_names['Hd']     = 'Hedenbergite'
    mineral_names['Hem']    = 'Hematite'
    mineral_names['Hyp']    = 'Hypersthene'
    mineral_names['Ilm']    = 'Ilmenite'
    mineral_names['Klpl']   = 'Kaliophilite'
    mineral_names['Lct']    = 'Leucite'
    mineral_names['Mag']    = 'Magnetite'
    mineral_names['Nph']    = 'Nepheline'
    mineral_names['Ol']     = 'Olivine'
    mineral_names['Or']     = 'Orthoclase'
    mineral_names['Prv']    = 'Perovskite'
    mineral_names['KMtsi']  = 'PotassiumMetasilicate'
    mineral_names['Py']     = 'Pyrite'
    mineral_names['Qz']     = 'Quartz'
    mineral_names['Rt']     = 'Rutile'
    mineral_names['NaCb']   = 'SodiumCarbonate'
    mineral_names['NaMtsi'] = 'SodiumMetasilicate'
    mineral_names['Thnd']   = 'Thenardite'
    mineral_names['Ttn']    = 'Titanite'
    mineral_names['Wo']     = 'Wollastonite'
    mineral_names['Zrn']    = 'Zircon'


    mineral_molwts = pd.Series()
    mineral_molwts['Acm'   ] = 231.00417*2.0
    mineral_molwts['Ab'    ] = 262.22301*2.0
    mineral_molwts['An'    ] = 278.20928
    mineral_molwts['Ap'    ] = 930.54816
    mineral_molwts['Cal'   ] = 100.0892
    mineral_molwts['CaOrsi'] = 172.2431
    mineral_molwts['Chr'   ] = 223.8366
    mineral_molwts['Crn'   ] = 101.9613
    mineral_molwts['Di'    ] = 216.55240
    mineral_molwts['En'    ] = 100.3887
    mineral_molwts['Fa'    ] = 203.7771
    mineral_molwts['Fs'    ] = 131.9307
    mineral_molwts['Fl'    ] =  78.076806
    mineral_molwts['Fo'    ] = 140.6931
    mineral_molwts['Hl'    ] =  58.44277
    mineral_molwts['Hd'    ] = 248.0944
    mineral_molwts['Hem'   ] = 159.6922
    mineral_molwts['Hyp'   ] =   0.0
    mineral_molwts['Ilm'   ] = 151.7452
    mineral_molwts['Klpl'  ] = 158.16294*2.0
    mineral_molwts['Lct'   ] = 218.24724*2.0
    mineral_molwts['Mag'   ] = 231.5386
    mineral_molwts['Nph'   ] = 142.05441*2.0
    mineral_molwts['Ol'    ] =   0.0
    mineral_molwts['Or'    ] = 278.33154*2.0
    mineral_molwts['Prv'   ] = 135.9782
    mineral_molwts['KMtsi' ] = 154.2803
    mineral_molwts['Py'    ] = 119.967
    mineral_molwts['Qz'    ] =  60.0843
    mineral_molwts['Rt'    ] =  79.8988
    mineral_molwts['NaCb'  ] = 105.98874
    mineral_molwts['NaMtsi'] = 122.06324
    mineral_molwts['Thnd'  ] = 142.03714
    mineral_molwts['Ttn'   ] = 196.0625
    mineral_molwts['Wo'    ] = 116.1637
    mineral_molwts['Zrn'   ] = 183.3031

    mineral_formulas = pd.Series()
    mineral_formulas['Acm'   ] = 'Na2O Fe2O3 4SiO2'
    mineral_formulas['Ab'    ] = 'Na2O Al2O3 6SiO2'
    mineral_formulas['An'    ] = 'CaO Al2O3 2SiO2'
    mineral_formulas['Ap'    ] = '3(3CaO P2O5)'
    mineral_formulas['Cal'   ] = 'CaO CO2'
    mineral_formulas['CaOrsi'] = '2CaO SiO2'
    mineral_formulas['Chr'   ] = 'FeO Cr2O3'
    mineral_formulas['Crn'   ] = 'Al2O3'
    mineral_formulas['Di'    ] = 'CaO MgO 2SiO2'
    mineral_formulas['En'    ] = 'MgO SiO2'
    mineral_formulas['Fa'    ] = '2FeO SiO2'
    mineral_formulas['Fs'    ] = 'FeO SiO2'
    mineral_formulas['Fl'    ] = 'CaF2'
    mineral_formulas['Fo'    ] = '2MgO SiO4'
    mineral_formulas['Hl'    ] = 'NaCl'
    mineral_formulas['Hd'    ] = 'CaO FeO 2SiO2'
    mineral_formulas['Hem'   ] = 'Fe2O3'
    mineral_formulas['Hyp'   ] = '(Mg, Fe)O SiO3'
    mineral_formulas['Ilm'   ] = 'FeO TiO2'
    mineral_formulas['Klpl'  ] = 'K2O Al2O3 2SiO2'
    mineral_formulas['Lct'   ] = 'K2O Al2O3 4SiO2'
    mineral_formulas['Mag'   ] = 'FeO Fe2O3'
    mineral_formulas['Nph'   ] = 'Na2O Al2O3 2SiO2'
    mineral_formulas['Ol'    ] = '2(Mg, Fe)O SiO2'
    mineral_formulas['Or'    ] = 'K2O Al2O3 6SiO2'
    mineral_formulas['Prv'   ] = 'CaO TiO2'
    mineral_formulas['KMtsi' ] = 'K2O SiO2'
    mineral_formulas['Py'    ] = 'FeS2'
    mineral_formulas['Qz'    ] = 'SiO2'
    mineral_formulas['Rt'    ] = 'TiO2'
    mineral_formulas['NaCb'  ] = 'Na2O CO2'
    mineral_formulas['NaMtsi'] = 'Na2O SiO2'
    mineral_formulas['Thnd'  ] = 'Na2O SO3'
    mineral_formulas['Ttn'   ] = 'CaO TiO2 SiO2'
    mineral_formulas['Wo'    ] = 'CaO SiO3'
    mineral_formulas['Zrn'   ] = 'ZrO2 SiO2'

    mineral_info = {}
    mineral_info['symbols'] = mineral_symbols
    mineral_info['names'] = mineral_names
    mineral_info['molwts'] = mineral_molwts
    mineral_info['formulas'] = mineral_formulas
    mineral_info['num'] = len(mineral_symbols)

    return mineral_info




def norm(oxide_wts):
    oxide_info = get_oxide_info()
    mineral_info = get_mineral_info()
    

    # Step 1
    oxide_mols = pd.Series(np.zeros(oxide_info['num']),
                           index=oxide_info['symbols'])
    # print(oxide_mols)
    for (iname, iwt) in oxide_wts.iteritems():
        oxide_mols[iname] = iwt/oxide_info['molwts'][iname]
        
    # oxide_mols *= oxide_wts/oxide_info['molwts']
    # print(oxide_mols)
    
    # oxide_mols = oxide_wts/oxide_info['molwts']
    min_mols = pd.Series(np.zeros(mineral_info['num']),
                         index=mineral_info['symbols'])
    # print(min_mols)


    # Step 2
    oxide_mols['FeO'] += oxide_mols['MnO'] + oxide_mols['NiO']
    oxide_mols['CaO'] += oxide_mols['BaO'] + oxide_mols['SrO']

    # Step 3a (fluorine is ignored in apatite)
    min_mols['Ap']  = oxide_mols['P2O5']/3.0
    oxide_mols['CaO']  -= 3.33*oxide_mols['P2O5']
    oxide_mols['P2O5']  = 0.0

    # Step 3b
    min_mols['Hl']  = 2.0*oxide_mols['Cl2O']
    oxide_mols['Na2O'] -= oxide_mols['Cl2O']
    oxide_mols['Cl2O']  = 0.0

    # Step 3c ('SO3' is oxidized sulfur)
    min_mols['Thnd']  = oxide_mols['SO3']
    oxide_mols['Na2O'] -= oxide_mols['SO3']
    oxide_mols['SO3']   = 0.0

    # Step 3d (S is reduced sulfur)
    min_mols['Py']  = oxide_mols['S']/2.0
    oxide_mols['FeO'] -= oxide_mols['S']/2.0
    oxide_mols['S']    = 0.0

    # Step 3e
    min_mols['Chr']  = oxide_mols['Cr2O3']
    oxide_mols['FeO'] -= oxide_mols['Cr2O3']
    oxide_mols['Cr2O3'] = 0.0

    # Step 3f
    if oxide_mols['FeO'] > oxide_mols['TiO2']:
        min_mols['Ilm'] = oxide_mols['TiO2']
        oxide_mols['FeO'] -= oxide_mols['TiO2']
        oxide_mols['TiO2'] = 0.0
    else:
        min_mols['Ilm'] = oxide_mols['FeO']
        oxide_mols['TiO2'] -= oxide_mols['FeO']
        oxide_mols['FeO']  = 0.0

    # Step 3g
    min_mols['Fl']  = oxide_mols['F2O']/2.0
    oxide_mols['CaO'] -= oxide_mols['F2O']/2.0
    oxide_mols['F2O']  = 0.0

    # Step 3h (assume 'CO2' is 1/2 calcite, 1/2 Cancrinite)
    min_mols['NaCb']  = oxide_mols['CO2']/2.0
    min_mols['Cal']   = oxide_mols['CO2']/2.0
    oxide_mols['Na2O'] -= oxide_mols['CO2']/2.0
    oxide_mols['CaO']  -= oxide_mols['CO2']/2.0
    oxide_mols['CO2']   = 0.0

    # Step 3i
    min_mols['Zrn'] = oxide_mols['ZrO2']
    oxide_mols['ZrO2'] = 0.0

    # Step 4a
    if oxide_mols['Al2O3'] >= oxide_mols['K2O']:
        min_mols['Or']  = oxide_mols['K2O']
        oxide_mols['Al2O3'] -= oxide_mols['K2O']
        oxide_mols['K2O'] = 0.0

        # Step 4c
        if oxide_mols['Al2O3'] >= oxide_mols['Na2O']:
            min_mols['Ab']  = oxide_mols['Na2O']
            oxide_mols['Al2O3'] -= oxide_mols['Na2O']
            oxide_mols['Na2O']  = 0.0

            # Step 4d/f
            if oxide_mols['Al2O3'] <= oxide_mols['CaO']:
                min_mols['An']  = oxide_mols['Al2O3']
                oxide_mols['CaO'] -= oxide_mols['Al2O3']
                oxide_mols['Al2O3'] = 0.0

            # Step 4e
            else:
                min_mols['An']  = oxide_mols['CaO']
                oxide_mols['Al2O3'] -= oxide_mols['CaO']
                oxide_mols['CaO'] = 0.0
                min_mols['Crn'] = oxide_mols['Al2O3']
                oxide_mols['Al2O3'] = 0.0

        # Step 4g
        else:
            min_mols['Ab']  = oxide_mols['Al2O3']
            oxide_mols['Na2O'] -= oxide_mols['Al2O3']
            oxide_mols['Al2O3'] = 0.0

    # Step 4b
    else:
        min_mols['Or'] = oxide_mols['Al2O3']
        min_mols['KMtsi'] = oxide_mols['K2O'] - oxide_mols['Al2O3']
        oxide_mols['K2O'] = 0.0
        oxide_mols['Al2O3'] = 0.0

    if oxide_mols['Na2O'] > 0.0:
        # Step 5a
        if oxide_mols['Na2O'] <= oxide_mols['Fe2O3']:
            min_mols['Acm']  = oxide_mols['Na2O']
            oxide_mols['Fe2O3'] -= oxide_mols['Na2O']
            oxide_mols['Na2O'] = 0.0

        # Step 5b
        else:
            min_mols['Acm'] = oxide_mols['Fe2O3']
            min_mols['NaMtsi'] = oxide_mols['Na2O'] - oxide_mols['Fe2O3']
            oxide_mols['Na2O'] = 0.0
            oxide_mols['Fe2O3'] = 0.0

    if oxide_mols['Fe2O3'] > 0.0:
        # Step 5c
        if oxide_mols['FeO'] >= oxide_mols['Fe2O3']:
            min_mols['Mag'] = oxide_mols['Fe2O3']
            oxide_mols['FeO'] -= oxide_mols['Fe2O3']
            oxide_mols['Fe2O3'] = 0.0

        # Step 5d
        else:
            min_mols['Mag'] = oxide_mols['FeO']
            min_mols['Hem'] = oxide_mols['Fe2O3'] - oxide_mols['FeO']
            oxide_mols['FeO'] = 0.0
            oxide_mols['Fe2O3'] = 0.0

    # Step 6
    ratio = 0.0
    if (oxide_mols['MgO']+oxide_mols['FeO']) > 0.0:
        ratio = oxide_mols['MgO']/(oxide_mols['MgO']+oxide_mols['FeO'])

    # Step 3f (continued)
    if oxide_mols['TiO2'] > 0.0:
        if oxide_mols['CaO'] >= oxide_mols['TiO2']:
            min_mols['Ttn']  = oxide_mols['TiO2']
            oxide_mols['CaO'] -= oxide_mols['TiO2']
            oxide_mols['TiO2'] = 0.0
        else:
            min_mols['Ttn'] = oxide_mols['CaO']
            min_mols['Rt'] = oxide_mols['TiO2'] - oxide_mols['CaO']
            oxide_mols['CaO'] = 0.0
            oxide_mols['TiO2'] = 0.0

    # Step 7a/b
    if (oxide_mols['MgO']+oxide_mols['FeO']) < oxide_mols['CaO']:
        min_mols['Di'] = oxide_mols['MgO'] + oxide_mols['FeO']
        min_mols['Wo'] = (oxide_mols['CaO'] - oxide_mols['MgO']
                          - oxide_mols['FeO'])
        oxide_mols['MgO'] = 0.0
        oxide_mols['FeO'] = 0.0
        oxide_mols['CaO'] = 0.0

    # Step 7c
    else:
        min_mols['Di'] = oxide_mols['CaO']
        min_mols['Hyp'] = (oxide_mols['MgO'] + oxide_mols['FeO']
                           - oxide_mols['CaO'])
        oxide_mols['MgO'] = 0.0
        oxide_mols['FeO'] = 0.0
        oxide_mols['CaO'] = 0.0

    # Step 8a
    sum = (min_mols['Ttn'] + 4*min_mols['Acm'] + min_mols['KMtsi']
           + min_mols['NaMtsi'] + 6*min_mols['Or'] + 6*min_mols['Ab']
           + 2*min_mols['An'] + 2*min_mols['Di'] + min_mols['Wo']
           + min_mols['Hyp'] + min_mols['Zrn'])
    oxide_mols['SiO2'] = oxide_mols['SiO2'] - sum

    # Step 8b
    if oxide_mols['SiO2'] >= 0.0:
        min_mols['Qz'] = oxide_mols['SiO2']

    # Step 8c
    else:
        oxide_mols['SiO2'] = oxide_mols['SiO2'] + min_mols['Hyp']
        if oxide_mols['SiO2'] >= 0.5*min_mols['Hyp']:
            sum = min_mols['Hyp']
            min_mols['Hyp'] = 2.0*oxide_mols['SiO2'] - sum
            min_mols['Ol'] = sum - oxide_mols['SiO2']
            oxide_mols['SiO2'] = 0.0

        # Step 8d/e
        else:
            min_mols['Ol'] = 0.5*min_mols['Hyp']
            oxide_mols['SiO2'] = (oxide_mols['SiO2'] - min_mols['Ol'] +
                                  min_mols['Ttn'] + 6*min_mols['Ab'])
            min_mols['Hyp'] = 0.0
            min_mols['Prv'] = min_mols['Ttn']
            min_mols['Ttn'] = 0.0

            if oxide_mols['SiO2'] >= 2.0*min_mols['Ab']:
                oxide_mols['Na2O'] = min_mols['Ab']
                min_mols['Ab'] = (
                    (oxide_mols['SiO2']-2*oxide_mols['Na2O'])/4 )
                min_mols['Nph'] = oxide_mols['Na2O'] - min_mols['Ab']
                oxide_mols['SiO2'] = 0.0

            # Step 8f
            else:
                min_mols['Nph'] = min_mols['Ab']
                min_mols['Ab'] = 0.0
                oxide_mols['SiO2'] = (
                    oxide_mols['SiO2'] - 2.0*min_mols['Nph']
                    + 6.0*min_mols['Or'] )

                if oxide_mols['SiO2'] >= 4.0*min_mols['Or']:
                    oxide_mols['K2O'] = min_mols['Or']
                    min_mols['Or'] = (
                        (oxide_mols['SiO2']-4.0*oxide_mols['K2O'])/2.0 )
                    min_mols['Lct'] = oxide_mols['K2O'] - min_mols['Or']
                    oxide_mols['SiO2'] = 0.0

                # Step 8g
                else:
                    if ((oxide_mols['SiO2']+min_mols['Wo']) <
                        4*min_mols['Or']):
                        min_mols['Lct'] = min_mols['Or']
                        min_mols['Or'] = 0.0
                        oxide_mols['SiO2'] = (
                            oxide_mols['SiO2'] - 4.0*min_mols['Lct']
                            + min_mols['Wo'] + 2.0*min_mols['Di'] )

                        # Step 8h
                        if ((oxide_mols['SiO2']+min_mols['Wo']
                             +min_mols['Ol'])
                            < (min_mols['Di']+min_mols['Ol'])):

                            oxide_mols['SiO2'] = (
                                oxide_mols['SiO2'] - min_mols['Di']
                                - min_mols['Wo']*0.5
                                + min_mols['Lct']*4.0)
                            min_mols['CaOrsi'] = (
                                0.5*(min_mols['Di']+min_mols['Wo']))
                            min_mols['Ol'] = (min_mols['Ol']
                                              + 0.5*min_mols['Di'])
                            min_mols['Wo'] = 0.0
                            min_mols['Di'] = 0.0
                            oxide_mols['K2O'] = min_mols['Lct']
                            min_mols['Lct'] = (oxide_mols['SiO2']-
                                               2.0*oxide_mols['K2O'])/2
                            min_mols['Klpl'] = (4.0*oxide_mols['K2O']
                                                -oxide_mols['SiO2'])/2
                            oxide_mols['SiO2'] = 0.0

                        else:
                            x = (2*oxide_mols['SiO2']-min_mols['Di']*2
                                 -min_mols['Wo'])/2
                            y = (min_mols['Di']-x)/2.0
                            z = (min_mols['Di']+min_mols['Wo']-x)/2
                            min_mols['Di'] = x
                            min_mols['Ol'] += y
                            min_mols['CaOrsi'] = z
                            min_mols['Wo'] = 0.0
                            oxide_mols['SiO2'] = 0.0
                    else:
                        min_mols['Lct'] = min_mols['Or']
                        min_mols['Or'] = 0.0
                        oxide_mols['SiO2'] = (oxide_mols['SiO2']
                                              - 4.0*min_mols['Lct'])
                        min_mols['CaOrsi']  = -oxide_mols['SiO2']
                        min_mols['Wo'] -= 2.0*min_mols['CaOrsi']
                        oxide_mols['SiO2'] = 0.0

    min_mols['Hd'] = (1.0-ratio)*min_mols['Di']
    min_mols['Di'] = ratio*min_mols['Di']
    min_mols['En'] = ratio*min_mols['Hyp']
    min_mols['Fs'] = (1.0-ratio)*min_mols['Hyp']
    min_mols['Fo'] = ratio*min_mols['Ol']
    min_mols['Fa'] = (1.0-ratio)*min_mols['Ol']

    min_wts = min_mols*mineral_info['molwts']
    # for imin_sym, imolwt in mineral_info['molwts']:
    #     # print(mineral_sym)
    #     # molwt = mineral_info['molwts'][mineral_sym]
    #     min_mols[imin_sym] *= imolwt

    min
    min_wts['Di'] = min_mols['Di'] + min_mols['Hd']
    min_wts['Hyp'] = min_mols['En'] + min_mols['Fs']
    min_wts['Ol'] = min_mols['Fo'] + min_mols['Fa']
    min_wts['Hd'] = 0.0
    min_wts['En'] = 0.0
    min_wts['Fs'] = 0.0
    min_wts['Fo'] = 0.0
    min_wts['Fa'] = 0.0

    return min_mols
