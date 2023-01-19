MELTS v. 1.1.0
==============

| Versions of MELTS implemented are:
| - MELTS v. 1.0.2 ➞ (rhyolite-MELTS, Gualda et al., 2012)
| - MELTS v. 1.1.0 ➞ (rhyolite-MELTS + new CO2, works at the ternary
  minimum)
| - MELTS v. 1.2.0 ➞ (rhyolite-MELTS + new H2O + new CO2)
| - pMELTS v. 5.6.1

Initialize tools and packages that are required to execute this notebook.
-------------------------------------------------------------------------

.. code:: ipython3

    from thermoengine import equilibrate
    import matplotlib.pyplot as plt
    import numpy as np
    %matplotlib inline

Create a MELTS v 1.1.0 instance.
--------------------------------

Rhyolite-MELTS version 1.0.2 is the default model.

.. code:: ipython3

    melts = equilibrate.MELTSmodel('1.1.0')

Optional: Generate some information about the implemented model.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    oxides = melts.get_oxide_names()
    phases = melts.get_phase_names()
    #print (oxides)
    #print (phases)

Required: Input initial composition of the system (liquid), in wt% or grams of oxides.
--------------------------------------------------------------------------------------

Early Bishop Tuff average melt inlusion composition

.. code:: ipython3

    feasible = melts.set_bulk_composition({'SiO2':  77.5, 
                                           'TiO2':   0.08, 
                                           'Al2O3': 12.5, 
                                           'Fe2O3':  0.207,
                                           'Cr2O3':  0.0, 
                                           'FeO':    0.473, 
                                           'MnO':    0.0,
                                           'MgO':    0.03, 
                                           'NiO':    0.0, 
                                           'CoO':    0.0,
                                           'CaO':    0.43, 
                                           'Na2O':   3.98, 
                                           'K2O':    4.88, 
                                           'P2O5':   0.0, 
                                           'H2O':    5.5,
                                           'CO2':    0.05})

Optional: Suppress phases that are not required in the simulation.
------------------------------------------------------------------

.. code:: ipython3

    b = melts.get_phase_inclusion_status()
    melts.set_phase_inclusion_status({'Nepheline':False, 'OrthoOxide':False})
    a = melts.get_phase_inclusion_status()
    for phase in b.keys():
        if b[phase] != a[phase]:
            print ("{0:<15s} Before: {1:<5s} After: {2:<5s}".format(phase, repr(b[phase]), repr(a[phase])))


.. parsed-literal::

    Nepheline       Before: True  After: False
    OrthoOxide      Before: True  After: False


Compute the equilibrium state at some specified T (°C) and P (MPa).
-------------------------------------------------------------------

Print status of the calculation.

.. code:: ipython3

    output = melts.equilibrate_tp(760.0, 175.0, initialize=True)
    (status, t, p, xmlout) = output[0]
    print (status, t, p)


.. parsed-literal::

    success, Optimal residual norm. 760.0 175.0


Summary output of equilibrium state …
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    melts.output_summary(xmlout)


.. parsed-literal::

    T (°C)       760.00
    P (MPa)      175.00
    Quartz            5.4862 (g)  SiO2                                                        
    Plagioclase       1.1738 (g)  K0.17Na0.76Ca0.07Al1.07Si2.93O8                             
    Liquid           88.5151 (g)  wt%:SiO2 73.57 TiO2  0.09 Al2O3 11.85 Fe2O3  0.23 Cr2O3  0.00 FeO  0.53 MnO  0.00 MgO  0.03
                                      NiO  0.00 CoO  0.00 CaO  0.43 Na2O  3.78 K2O  4.60 P2O5  0.00 H2O  4.88 CO2  0.00
    Sanidine          9.2226 (g)  K0.48Na0.50Ca0.02Al1.02Si2.98O8                             
    Fluid             1.2323 (g)  X:H2O 0.984 CO2 0.016
                                     


Output thermodynamic properties of any phase present in the system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

… or the sum of all phases in the system

.. code:: ipython3

    props = melts.get_list_of_properties()
    for prop in props:
        print ("{0:<20s} {1:13.6e} {2:<10s}".format(prop, melts.get_property_of_phase(xmlout,'System', prop), \
                                                    melts.get_units_of_property(prop)))


.. parsed-literal::

    Mass                  1.056300e+02 g         
    GibbsFreeEnergy      -1.724198e+06 J         
    Enthalpy             -1.447163e+06 J         
    Entropy               2.681463e+02 J/K       
    HeatCapacity          1.420167e+02 J/K       
    DcpDt                 5.784923e-03 J/K^2     
    Volume                4.948637e+00 J/bar     
    DvDt                  1.067823e-03 J/bar-K   
    DvDp                 -1.967218e-04 J/bar^2   
    D2vDt2                2.456499e-07 J/bar-K^2 
    D2vDtDp              -3.977640e-07 J/bar^2-K 
    D2vDp2                1.873064e-07 J/bar^3   
    Density               2.134527e+00 g/cm^3    
    Alpha                 2.157813e-04 1/K       
    Beta                  3.975272e-05 1/bar     
    K                     2.515551e+00 GPa       
    K'                    2.295153e+01 none      
    Gamma                 3.990411e-02 none      


Output chemical affinities and potential compositions of undersaturated phases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    dict = melts.get_dictionary_of_affinities(xmlout, sort=True)
    for phase in dict:
        (affinity, formulae) = dict[phase]
        print ("{0:<20s} {1:10.2f} {2:<60s}".format(phase, affinity, formulae))


.. parsed-literal::

    Spinel                  1728.65 Fe''1.00Mg-0.00Fe'''-0.00Al1.00Cr1.00Ti-0.00O4              
    Leucite                 1877.28 K0.50Na0.50AlSi2O5.50(OH)1.00                               
    Tridymite               2141.04 SiO2                                                        
    Cristobalite            2204.71 SiO2                                                        
    Magnetite               4416.17 Fe3O4                                                       
    Ilmenite ss             5665.90 Mn-0.00Fe''-0.00Mg0.50Fe'''1.00Al-0.00Ti0.50O3              
    Orthopyroxene           6377.10 Na-0.00Ca0.50Fe''-0.00Mg1.50Fe'''-0.00Ti-0.00Al-0.00Si2.00O6
    Coesite                 6702.07 SiO2                                                        
    Pigeonite               6796.30 Na-0.00Ca0.50Fe''-0.00Mg1.50Fe'''-0.00Ti-0.00Al-0.00Si2.00O6
    Olivine                 7457.21 (Ca-0.00Mg0.00Fe''0.50Mn0.50Co0.00Ni0.00)2SiO4              
    Ilmenite                9211.89 FeTiO3                                                      
    Fayalite                9389.14 Fe2SiO4                                                     
    Rutile                 12105.23 TiO2                                                        
    Sillimanite            14226.37 Al2SiO5                                                     
    Andalusite             14585.72 Al2SiO5                                                     
    Corundum               14607.25 Al2O3                                                       
    Garnet                 17128.41 (Ca0.00Fe''0.00Mg1.00)3Al2Si3O12                            
    Hematite               17719.32 Fe2O3                                                       
    Biotite                18137.52 K(Fe''0.00Mg1.00)3AlSi3O10(OH)2                             
    Kyanite                18246.67 Al2SiO5                                                     
    Nepheline              20021.77 NaAlSiO4                                                    
    Sphene                 22664.99 CaTiSiO5                                                    
    Muscovite              25226.95 KAl2Si3AlO10(OH)2                                           
    OrthoOxide             26989.92 Fe''0.50Mg-0.00Fe'''1.00Ti1.50O5                            
    Phlogopite             27050.81 KMg3AlSi3O10(OH)2                                           
    Forsterite             28071.11 Mg2SiO4                                                     
    Cummingtonite          33263.63 Ca-0.00Fe3.50Mg3.50Si8O22(OH)2                              
    Kalsilite              33674.89 KAlSiO4                                                     
    Anthophyllite          35362.91 Ca-0.00Fe7.00Mg0.00Si8O22(OH)2                              
    Perovskite             37089.30 CaTiO3                                                      
    Periclase              42573.76 MgO                                                         
    Nepheline ss           43695.20 Na2.00K2.00Ca-0.00[]-0.00Al4.00Si4.00O16                    
    Melilite               47991.88 Na-0.00Ca2.00Al1.00Mg0.50Fe-0.00Si1.50O7                    
    Hornblende             51692.58 NaCa2Mg2.00Fe2+2.00Al1.00Fe3+-0.00Al2Si6O22(OH)2            
    Kalsilite ss           51807.73 Na2.00K2.00Ca-0.00[]-0.00Al4.00Si4.00O16                    
    Aenigmatite            54045.03 Na2Fe5TiSi6O20                                              
    Solid Alloy            68565.93 Fe1.00Ni0.00                                                
    Liquid Alloy           73162.93 Fe1.00Ni0.00                                                
    Akermanite             79473.52 Ca2MgSi2O7                                                  
    Gehlenite              83674.97 Ca2Al2SiO7                                                  
    Lime                  113223.36 CaO                                                         
    Aegirine              141031.17 NaFeSi2O6                                                   
    Actinolite            999999.00 Ca-0.00Fe3.50Mg3.50Si8O22(OH)2                              
    Whitlockite           999999.00 Ca3(PO4)2                                                   
    Apatite               999999.00 Ca5(PO4)3OH                                                 
    Chromite              999999.00 FeCr2O4                                                     
    Augite                999999.00 Na-0.00Ca0.50Fe''-0.00Mg1.50Fe'''-0.00Ti-0.00Al-0.00Si2.00O6
    Titanaugite           999999.00 Na-0.00Ca0.50Fe''-0.00Mg1.50Fe'''-0.00Ti-0.00Al-0.00Si2.00O6
    Panunzite             999999.00 Na2.00K2.00Ca-0.00[]-0.00Al4.00Si4.00O16                    


Run the sequence of calculations along a T, P=constant path:
------------------------------------------------------------

Output is sent to an Excel file and plotted in the notebook

.. code:: ipython3

    number_of_steps = 20
    t_increment_of_steps = -1.0
    p_increment_of_steps = 0.0
    
    plotPhases = ['Liquid', 'Sanidine', 'Plagioclase', 'Quartz', 'Fluid']
    # matplotlib colors b : blue, g : green, r : red, c : cyan, m : magenta, y : yellow, k : black, w : white.
    plotColors = [ 'ro', 'bo', 'go', 'co', 'mo']
    
    wb = melts.start_excel_workbook_with_sheet_name(sheetName="Summary")
    melts.update_excel_workbook(wb, xmlout)
    
    n = len(plotPhases)
    xPlot = np.zeros(number_of_steps+1)
    yPlot = np.zeros((n, number_of_steps+1))
    xPlot[0] = t
    for i in range (0, n):
        yPlot[i][0] = melts.get_property_of_phase(xmlout, plotPhases[i])
    
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim([min(t, t+t_increment_of_steps*number_of_steps), max(t, t+t_increment_of_steps*number_of_steps)])
    ax.set_ylim([0., 100.])
    graphs = []
    for i in range (0, n):
        graphs.append(ax.plot(xPlot, yPlot[i], plotColors[i]))
    handle = []
    for (graph,) in graphs:
        handle.append(graph)
    ax.legend(handle, plotPhases, loc='upper left')
    
    for i in range (1, number_of_steps):
        output = melts.equilibrate_tp(t+t_increment_of_steps, p+p_increment_of_steps)
        (status, t, p, xmlout) = output[0]
        print ("{0:<30s} {1:8.2f} {2:8.2f}".format(status, t, p))
        xPlot[i] = t
        for j in range (0, n):
            yPlot[j][i] = melts.get_property_of_phase(xmlout, plotPhases[j])
        j = 0
        for (graph,) in graphs:
            graph.set_xdata(xPlot)
            graph.set_ydata(yPlot[j])
            j = j + 1
        fig.canvas.draw()
        melts.update_excel_workbook(wb, xmlout)
    
    melts.write_excel_workbook(wb, "MELTSv102summary.xlsx")


.. parsed-literal::

    success, Optimal residual norm.   759.00   175.00
    success, Minimal energy computed.   758.00   175.00
    success, Optimal residual norm.   757.00   175.00
    success, Minimal energy computed.   756.00   175.00
    success, Minimal energy computed.   755.00   175.00
    success, Minimal energy computed.   754.00   175.00
    success, Minimal energy computed.   753.00   175.00
    success, Minimal energy computed.   752.00   175.00
    success, Minimal energy computed.   751.00   175.00
    success, Optimal residual norm.   750.00   175.00
    success, Minimal energy computed.   749.00   175.00
    success, Optimal residual norm.   748.00   175.00
    success, Minimal energy computed.   747.00   175.00
    success, Minimal energy computed.   746.00   175.00
    success, Optimal residual norm.   745.00   175.00
    success, Minimal energy computed.   744.00   175.00
    success, Minimal energy computed.   743.00   175.00
    success, Minimal energy computed.   742.00   175.00
    success, Minimal energy computed.   741.00   175.00



.. image:: MELTS-v1.1.0-equilibrium_files/MELTS-v1.1.0-equilibrium_20_1.png


.. code:: ipython3

    melts.output_summary(xmlout)


.. parsed-literal::

    T (°C)       741.00
    P (MPa)      175.00
    Quartz           34.1575 (g)  SiO2                                                        
    Spinel            0.4978 (g)  Fe''1.38Mg0.02Fe'''1.15Al0.04Cr0.00Ti0.40O4                 
    Liquid            1.6490 (g)  wt%:SiO2 74.49 TiO2  0.10 Al2O3  8.04 Fe2O3  0.17 Cr2O3  0.00 FeO  1.30 MnO  0.00 MgO 
                                      0.05 NiO  0.00 CoO  0.00 CaO  0.84 Na2O  5.84 K2O  3.68 P2O5  0.00 H2O  5.50 CO2 
    Olivine           0.1280 (g)  (Ca0.00Mg0.09Fe''0.91Mn0.00Co0.00Ni0.00)2SiO4               
    Ilmenite ss       0.0149 (g)  Mn0.00Fe''0.89Mg0.04Fe'''0.11Al0.04Ti0.92O3                 
    Orthopyroxene     0.3254 (g)  Na0.00Ca0.02Fe''1.52Mg0.43Fe'''0.02Ti0.00Al0.05Si1.97O6     
    Sanidine         47.5465 (g)  K0.54Na0.45Ca0.01Al1.01Si2.99O8                             
    Fluid             5.4593 (g)  X:H2O 0.996 CO2 0.004
                                     
    Plagioclase      15.8515 (g)  K0.14Na0.78Ca0.09Al1.09Si2.91O8                             


