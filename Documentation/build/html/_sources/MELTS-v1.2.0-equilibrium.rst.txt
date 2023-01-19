MELTS v. 1.2.0
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

Create a MELTS v 1.2.0 instance.
--------------------------------

Rhyolite-MELTS version 1.0.2 is the default model.

.. code:: ipython3

    melts = equilibrate.MELTSmodel('1.2.0')

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
    Plagioclase      11.1254 (g)  K0.09Na0.85Ca0.07Al1.07Si2.93O8                             
    Spinel            0.0752 (g)  Fe''1.11Mg0.02Fe'''1.64Al0.10Cr0.00Ti0.13O4                 
    Liquid           93.9072 (g)  wt%:SiO2 74.64 TiO2  0.08 Al2O3 10.88 Fe2O3  0.17 Cr2O3  0.00 FeO  0.48 MnO  0.00 MgO  0.03
                                      NiO  0.00 CoO  0.00 CaO  0.29 Na2O  3.06 K2O  5.01 P2O5  0.00 H2O  5.35 CO2  0.01
    Fluid             0.5223 (g)  X:H2O 0.964 CO2 0.036
                                     


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
    GibbsFreeEnergy      -1.724571e+06 J         
    Enthalpy             -1.445192e+06 J         
    Entropy               2.695581e+02 J/K       
    HeatCapacity          1.414445e+02 J/K       
    DcpDt                 6.430431e-03 J/K^2     
    Volume                4.881793e+00 J/bar     
    DvDt                  6.583345e-04 J/bar-K   
    DvDp                 -9.668888e-05 J/bar^2   
    D2vDt2               -3.343413e-08 J/bar-K^2 
    D2vDtDp              -1.352812e-07 J/bar^2-K 
    D2vDp2                7.450556e-08 J/bar^3   
    Density               2.163754e+00 g/cm^3    
    Alpha                 1.348551e-04 1/K       
    Beta                  1.980602e-05 1/bar     
    K                     5.048971e+00 GPa       
    K'                    3.790586e+01 none      
    Gamma                 4.976697e-02 none      


Output chemical affinities and potential compositions of undersaturated phases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    dict = melts.get_dictionary_of_affinities(xmlout, sort=True)
    for phase in dict:
        (affinity, formulae) = dict[phase]
        print ("{0:<20s} {1:10.2f} {2:<60s}".format(phase, affinity, formulae))


.. parsed-literal::

    Quartz                   230.43 SiO2                                                        
    Leucite                 1200.94 K0.00Na1.00AlSi2O6.00(OH)0.00                               
    Tridymite               2371.47 SiO2                                                        
    Cristobalite            2435.13 SiO2                                                        
    Magnetite               2517.18 Fe3O4                                                       
    Ilmenite ss             4450.46 Mn0.00Fe''0.96Mg0.00Fe'''0.00Al0.08Ti0.96O3                 
    Coesite                 6932.50 SiO2                                                        
    Orthopyroxene           7359.78 Na0.00Ca1.00Fe''0.93Mg0.00Fe'''0.01Ti0.00Al0.12Si1.94O6     
    Pigeonite               7679.63 Na0.01Ca0.99Fe''0.93Mg0.00Fe'''0.03Ti0.00Al0.11Si1.94O6     
    Olivine                 7959.14 (Ca0.01Mg1.00Fe''0.00Mn0.00Co0.00Ni0.00)2SiO4               
    Ilmenite                8630.75 FeTiO3                                                      
    Fayalite                9553.05 Fe2SiO4                                                     
    Rutile                 11557.34 TiO2                                                        
    Sillimanite            11688.61 Al2SiO5                                                     
    Corundum               11839.07 Al2O3                                                       
    Andalusite             12047.96 Al2SiO5                                                     
    Kyanite                15708.92 Al2SiO5                                                     
    Hematite               15853.59 Fe2O3                                                       
    Garnet                 16785.03 (Ca0.00Fe''0.00Mg1.00)3Al2Si3O12                            
    Nepheline              18951.56 NaAlSiO4                                                    
    Biotite                21198.06 K(Fe''0.00Mg1.00)3AlSi3O10(OH)2                             
    Muscovite              25037.67 KAl2Si3AlO10(OH)2                                           
    OrthoOxide             25263.60 Fe''0.00Mg1.00Fe'''0.00Ti2.00O5                             
    Sphene                 27674.83 CaTiSiO5                                                    
    Forsterite             31263.21 Mg2SiO4                                                     
    Phlogopite             34072.22 KMg3AlSi3O10(OH)2                                           
    Kalsilite              35620.86 KAlSiO4                                                     
    Cummingtonite          39205.76 Ca2.00Fe0.00Mg5.00Si8O22(OH)2                               
    Nepheline ss           41011.01 Na2.99K0.00Ca0.01[]1.00Al3.01Si4.99O16                      
    Anthophyllite          41250.21 Ca2.00Fe0.00Mg5.00Si8O22(OH)2                               
    Perovskite             41868.71 CaTiO3                                                      
    Periclase              44054.59 MgO                                                         
    Kalsilite ss           51092.31 Na2.99K0.00Ca0.01[]1.00Al3.01Si4.99O16                      
    Melilite               52369.61 Na1.53Ca0.47Al0.00Mg0.00Fe0.24Si2.76O7                      
    Aenigmatite            54880.30 Na2Fe5TiSi6O20                                              
    Hornblende             61274.22 NaCa2Mg4.00Fe2+0.00Al0.00Fe3+1.00Al2Si6O22(OH)2             
    Solid Alloy            70331.88 Fe1.00Ni0.00                                                
    Liquid Alloy           74928.88 Fe1.00Ni0.00                                                
    Gehlenite              91791.80 Ca2Al2SiO7                                                  
    Akermanite             92069.79 Ca2MgSi2O7                                                  
    Lime                  118550.65 CaO                                                         
    Aegirine              140642.61 NaFeSi2O6                                                   
    Actinolite            999999.00 Ca-0.00Fe0.00Mg7.00Si8O22(OH)2                              
    Whitlockite           999999.00 Ca3(PO4)2                                                   
    Apatite               999999.00 Ca5(PO4)3OH                                                 
    Chromite              999999.00 FeCr2O4                                                     
    Sanidine              999999.00 K1.00Na0.00Ca0.00Al1.00Si3.00O8                             
    Augite                999999.00 Na0.01Ca0.99Fe''0.93Mg0.00Fe'''0.03Ti0.00Al0.11Si1.94O6     
    Titanaugite           999999.00 Na0.01Ca0.99Fe''0.93Mg0.00Fe'''0.03Ti0.00Al0.11Si1.94O6     
    Panunzite             999999.00 Na2.99K0.00Ca0.01[]1.00Al3.01Si4.99O16                      


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
    success, Optimal residual norm.   758.00   175.00
    success, Optimal residual norm.   757.00   175.00
    success, Optimal residual norm.   756.00   175.00
    success, Optimal residual norm.   755.00   175.00
    success, Optimal residual norm.   754.00   175.00
    success, Optimal residual norm.   753.00   175.00
    success, Optimal residual norm.   752.00   175.00
    success, Optimal residual norm.   751.00   175.00
    success, Optimal residual norm.   750.00   175.00
    success, Optimal residual norm.   749.00   175.00
    success, Optimal residual norm.   748.00   175.00
    success, Optimal residual norm.   747.00   175.00
    success, Minimal energy computed.   746.00   175.00
    success, Optimal residual norm.   745.00   175.00
    success, Optimal residual norm.   744.00   175.00
    success, Minimal energy computed.   743.00   175.00
    success, Optimal residual norm.   742.00   175.00
    success, Optimal residual norm.   741.00   175.00



.. image:: MELTS-v1.2.0-equilibrium_files/MELTS-v1.2.0-equilibrium_20_1.png


.. code:: ipython3

    melts.output_summary(xmlout)


.. parsed-literal::

    T (°C)       741.00
    P (MPa)      175.00
    Plagioclase      17.9599 (g)  K0.11Na0.83Ca0.06Al1.06Si2.94O8                             
    Spinel            0.1876 (g)  Fe''1.15Mg0.02Fe'''1.57Al0.09Cr0.00Ti0.17O4                 
    Liquid           83.4833 (g)  wt%:SiO2 74.94 TiO2  0.08 Al2O3 10.60 Fe2O3  0.12 Cr2O3  0.00 FeO  0.48 MnO  0.00 MgO  0.04
                                      NiO  0.00 CoO  0.00 CaO  0.26 Na2O  2.67 K2O  5.42 P2O5  0.00 H2O  5.38 CO2  0.00
    Quartz            2.9467 (g)  SiO2                                                        
    Fluid             1.0524 (g)  X:H2O 0.981 CO2 0.019
                                     


