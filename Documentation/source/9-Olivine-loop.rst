Olivine Phase Loop
==================

| This notebook demonstrates calculation of the olivine liquid-solid
  phase loop under the assumption that both phases behave as ideal
  solutions.
| The workflow is: - Use the **coder** module to generate endmember
  properties of both solutions; the only thermodynamic properties that
  are specified are the enthalpy and entropy of fusion - Use the
  **coder** module to generate solid and liquid solution properties -
  Import the generated code using the **model** module - Use the
  **equilibrate** module to compute the liquid-solid phase loop - Plot
  results

.. code:: ipython3

    import numpy as np
    import scipy as sp
    import sympy as sym
    import matplotlib.pyplot as plt
    from thermoengine import model, equilibrate, coder
    %matplotlib inline

Endmember properties
--------------------

Write code into a working subdirectory

.. code:: ipython3

    working_dir = "working"
    !mkdir -p {working_dir}
    %cd {working_dir}


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Equilibrate/working


Model generation function

.. code:: ipython3

    def make_endmembers(module='none', name='none', formula='none', Hrefvalue=0.0, Srefvalue=0.0):
        mdl = coder.StdStateModel()
        T = mdl.get_symbol_for_t()
        GPr,Href,Sref = sym.symbols('GPr Href Sref')
        GPr = Href - T*Sref
        params = [('Href', 'J', Href), ('Sref', 'J/K', Sref)]
        mdl.add_expression_to_model(GPr, params)
        mdl.set_module_name(module)
        paramValues = {'Href':Hrefvalue, 'Sref':Srefvalue, 'T_r':298.15, 'P_r':1.0}
        mdl.create_code_module(phase=name, formula=formula, params=paramValues, 
                               module_type='calib', silent=True)

Forsterite Solid
~~~~~~~~~~~~~~~~

.. code:: ipython3

    make_endmembers(module='OlvSolid', name='Fo', formula='Mg(2)Si(1)O(4)', Hrefvalue=-100000.0, Srefvalue=0.0)
    %cp OlvSolid.pyx endmembersolids.pyx

Fayalite Solid
~~~~~~~~~~~~~~

.. code:: ipython3

    make_endmembers(module='OlvSolid', name='Fa', formula='Fe(2)Si(1)O(4)', Hrefvalue=-100000.0, Srefvalue=0.0)
    %cat OlvSolid.pyx >> endmembersolids.pyx

Forsterite Liquid
~~~~~~~~~~~~~~~~~

Fusion temperature is 2163 K, entropy is 57.2 J/K

.. code:: ipython3

    make_endmembers(module='OlvLiquid', name='Fo', formula='Mg(2)Si(1)O(4)', Hrefvalue=-100000.0+57.2*2163.0, Srefvalue=57.2)
    %cp OlvLiquid.pyx endmemberliquids.pyx

Fayalite Liquid
~~~~~~~~~~~~~~~

Fusion temperature is 1490 K, entropy is 59.9 J/K

.. code:: ipython3

    make_endmembers(module='OlvLiquid', name='Fa', formula='Fe(2)Si(1)O(4)', Hrefvalue=-100000.0+59.9*1490.0, Srefvalue=59.9)
    %cat OlvLiquid.pyx >> endmemberliquids.pyx

Solution Properties
-------------------

Model generation function

.. code:: ipython3

    def make_solution(module='none', name='none', endmembers=[]):
        c = 2
        mdl = coder.SimpleSolnModel(nc=c)
        n = mdl.n
        nT = mdl.nT
        X = n/nT
        T = mdl.get_symbol_for_t()
        mu = mdl.mu
        G_ss = (n.transpose()*mu)[0]
        S_config,R = sym.symbols('S_config R')
        S_config = 0
        for i in range(0,c):
            S_config += X[i]*sym.log(X[i])
        S_config *= -R*nT
        G_config = sym.simplify(-T*S_config)
        G = G_ss + G_config
        mdl.add_expression_to_model(G, [('dummy', 'none', sym.symbols('dummy'))])
        mdl.module = module
        mdl.formula_string = 'Mg[Mg]Fe[Fe]Si[Si]O4'
        mdl.conversion_string = ['[0]=[Mg]', '[1]=[Fe]']
        mdl.test_string = ['[0] > 0.0', '[1] > 0.0']
        mdl.create_code_module(phase=name, params={'dummy':0.0, 'T_r':298.15, 'P_r':1.0}, 
                               endmembers=endmembers, 
                               prefix="cy", module_type='calib', silent=True)

Solid solution
~~~~~~~~~~~~~~

.. code:: ipython3

    make_solution(module='OlvSolid', name='Olivine', endmembers=['Fo_OlvSolid', 'Fa_OlvSolid'])
    %cat endmembersolids.pyx >> OlvSolid.pyx

Liquid solution
~~~~~~~~~~~~~~~

.. code:: ipython3

    make_solution(module='OlvLiquid', name='Liquid', endmembers=['Fo_OlvLiquid', 'Fa_OlvLiquid'])
    %cat endmemberliquids.pyx >> OlvLiquid.pyx

.. code:: ipython3

    import OlvSolid
    import OlvLiquid
    %cd ..


.. parsed-literal::

    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/Cython/Compiler/Main.py:369: FutureWarning: Cython directive 'language_level' not set, using 2 for now (Py2). This will change in a later release! File: /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Equilibrate/working/OlvSolid.pyx
      tree = Parsing.p_module(s, pxd, full_module_name)
    /Users/ghiorso/anaconda3/lib/python3.7/site-packages/Cython/Compiler/Main.py:369: FutureWarning: Cython directive 'language_level' not set, using 2 for now (Py2). This will change in a later release! File: /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Equilibrate/working/OlvLiquid.pyx
      tree = Parsing.p_module(s, pxd, full_module_name)


.. parsed-literal::

    /Users/ghiorso/Documents/ARCHIVE_XCODE/ThermoEngine/Notebooks/Equilibrate


Set up phase loop calculation
-----------------------------

.. code:: ipython3

    modelDBsol = model.Database(database="CoderModule", calib='calib', 
                                phase_tuple=('OlvSolid', {
                                    'Ol':['Olivine','solution'],
                                    'Fo':['Fo','pure'],
                                    'Fa':['Fa','pure']
                                }))
    modelDBliq = model.Database(database="CoderModule", calib='calib', 
                                phase_tuple=('OlvLiquid', {
                                    'Liq':['Liquid','solution'],
                                    'Fo':['Fo','pure'],
                                    'Fa':['Fa','pure']
                                }))


.. parsed-literal::

    Solution phase code generated by the coder module does not yet provide information on solution species. Species are proxied by components.
    Solution phase code generated by the coder module does not yet provide information on species properties. Species are proxied by components.
    Solution phase code generated by the coder module does not yet provide information on solution species. Species are proxied by components.
    Solution phase code generated by the coder module does not yet provide information on species properties. Species are proxied by components.


.. code:: ipython3

    olivine = modelDBsol.get_phase("Ol")
    liquid = modelDBliq.get_phase("Liq")

.. code:: ipython3

    elm_sys = ['O','Mg','Si','Fe']
    phs_sys = [liquid, olivine]

Compute the loop
~~~~~~~~~~~~~~~~

.. code:: ipython3

    xFoSol = [1.0]
    xFoLiq = [1.0]
    tC = [2163.0-273.15]
    p = 1.0
    for i in range(1,20):
        XFo = 1.0 - i*0.05
        XFa = 1.0 - XFo
        blk_cmp = np.array([4.0*(XFo+XFa), 2.0*XFo, XFo+XFa, 2.0*XFa])
        equil = equilibrate.Equilibrate(elm_sys, phs_sys)
        t = 2163.0*XFo + 1490.0*XFa
        state = equil.execute(t, p, bulk_comp=blk_cmp, debug=0, stats=True)
        state.print_state()
        tC.append(t-273.15)
        xFoSol.append(state.compositions(phase_name='Olivine', units='mole_frac')[0])
        xFoLiq.append(state.compositions(phase_name='Liquid', units='mole_frac')[0])
    xFoSol.append(0.0)
    xFoLiq.append(0.0)
    tC.append(1490.0-273.15)


.. parsed-literal::

    Add: Olivine
    Quad (000) norm:  2.1104254762564e+00 Lin (016) step:  3.8779326997502e-01 func: -1.0402273326438e+05
    Quad (001) norm:  3.2418611320216e-02 Lin (017) step:  1.7951470449315e+00 func: -1.0405375729263e+05
    Quad (002) norm:  1.4763143793992e-02 Lin (010) step:  9.6932535179819e-01 func: -1.0405401325700e+05
    Quad (003) norm:  4.0623266711652e-05 Lin (029) step:  9.9979197987555e-01 func: -1.0405401326837e+05
    Quad (004) norm:  7.3677199507684e-09 Lin (039) step:  8.4138217960888e-01 func: -1.0405401326837e+05
    Quad (005) norm:  1.1686594034085e-09 Lin (040) step: -7.3123984754483e-01 func: -1.0405401326837e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1856.20 °C, P =        0.1 MPa
    Liquid          moles:   0.361955 grams:  53.547
                 Fo form:  Mg2SiO4        X:  0.8851  wt%    SiO2   40.61
                 Fa form:  Fe2SiO4        X:  0.1149  wt%     FeO   11.16
                                                      wt%     MgO   48.23
    Olivine         moles:   0.638045 grams:  90.300
                 Fo form:  Mg2SiO4        X:  0.9868  wt%    SiO2   42.45
                 Fa form:  Fe2SiO4        X:  0.0132  wt%     FeO    1.34
                                                      wt%     MgO   56.20
    Add: Olivine
    Quad (000) norm:  2.0114451413753e+00 Lin (015) step:  3.9174467592009e-01 func: -1.0663613100463e+05
    Quad (001) norm:  3.1142185393805e-02 Lin (032) step:  1.8317912801402e+00 func: -1.0670046985749e+05
    Quad (002) norm:  1.5416604556166e-02 Lin (019) step:  9.6691766198867e-01 func: -1.0670104487450e+05
    Quad (003) norm:  4.5064898567581e-05 Lin (028) step:  9.9975861375761e-01 func: -1.0670104490454e+05
    Quad (004) norm:  9.6146541687342e-09 Lin (039) step:  8.5998149639898e-01 func: -1.0670104490454e+05
    Quad (005) norm:  1.3462267133040e-09 Lin (039) step:  8.1909231079668e-01 func: -1.0670104490454e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1822.55 °C, P =        0.1 MPa
    Liquid          moles:   0.376308 grams:  58.172
                 Fo form:  Mg2SiO4        X:  0.7798  wt%    SiO2   38.87
                 Fa form:  Fe2SiO4        X:  0.2202  wt%     FeO   20.47
                                                      wt%     MgO   40.66
    Olivine         moles:   0.623692 grams:  88.829
                 Fo form:  Mg2SiO4        X:  0.9725  wt%    SiO2   42.19
                 Fa form:  Fe2SiO4        X:  0.0275  wt%     FeO    2.77
                                                      wt%     MgO   55.04
    Add: Olivine
    Quad (000) norm:  1.9106239644072e+00 Lin (016) step:  3.9582177358734e-01 func: -1.0863573839282e+05
    Quad (001) norm:  3.0573320893983e-02 Lin (014) step:  1.8719312484472e+00 func: -1.0873576578541e+05
    Quad (002) norm:  1.6092789244022e-02 Lin (014) step:  9.6415282496294e-01 func: -1.0873673638519e+05
    Quad (003) norm:  5.1487274533657e-05 Lin (029) step:  9.9945717817509e-01 func: -1.0873673644510e+05
    Quad (004) norm:  3.5575530590706e-09 Lin (038) step:  1.6519526824590e+00 func: -1.0873673644510e+05
    Quad (005) norm:  2.3193563531757e-09 Lin (042) step: -3.6785150611398e-01 func: -1.0873673644510e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1788.90 °C, P =        0.1 MPa
    Liquid          moles:   0.391316 grams:  62.871
                 Fo form:  Mg2SiO4        X:  0.6834  wt%    SiO2   37.40
                 Fa form:  Fe2SiO4        X:  0.3166  wt%     FeO   28.31
                                                      wt%     MgO   34.29
    Olivine         moles:   0.608684 grams:  87.285
                 Fo form:  Mg2SiO4        X:  0.9571  wt%    SiO2   41.90
                 Fa form:  Fe2SiO4        X:  0.0429  wt%     FeO    4.30
                                                      wt%     MgO   53.80
    Add: Olivine
    Quad (000) norm:  1.8079136182821e+00 Lin (024) step:  4.0003822981413e-01 func: -1.1019554484785e+05
    Quad (001) norm:  3.0787646794130e-02 Lin (016) step:  1.9160707731592e+00 func: -1.1033368599857e+05
    Quad (002) norm:  1.6791400729548e-02 Lin (020) step:  9.6095788326558e-01 func: -1.1033514496278e+05
    Quad (003) norm:  6.0794171941693e-05 Lin (024) step:  9.9954277934166e-01 func: -1.1033514506973e+05
    Quad (004) norm:  9.9720374388972e-09 Lin (040) step:  4.9983406065890e-01 func: -1.1033514506973e+05
    Quad (005) norm:  4.9876737609594e-09 Lin (047) step: -2.1334834772453e-02 func: -1.1033514506973e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1755.25 °C, P =        0.1 MPa
    Liquid          moles:   0.406994 grams:  67.643
                 Fo form:  Mg2SiO4        X:  0.5956  wt%    SiO2   36.15
                 Fa form:  Fe2SiO4        X:  0.4044  wt%     FeO   34.96
                                                      wt%     MgO   28.89
    Olivine         moles:   0.593006 grams:  85.667
                 Fo form:  Mg2SiO4        X:  0.9403  wt%    SiO2   41.59
                 Fa form:  Fe2SiO4        X:  0.0597  wt%     FeO    5.94
                                                      wt%     MgO   52.47
    Add: Olivine
    Quad (000) norm:  1.7032709372590e+00 Lin (014) step:  4.0441189414102e-01 func: -1.1139918109725e+05
    Quad (001) norm:  3.1794373696592e-02 Lin (029) step:  1.9648182746974e+00 func: -1.1157786603177e+05
    Quad (002) norm:  1.7512076248931e-02 Lin (021) step:  9.5724546055780e-01 func: -1.1157992569536e+05
    Quad (003) norm:  7.4166851734302e-05 Lin (022) step:  9.9954766440304e-01 func: -1.1157992587561e+05
    Quad (004) norm:  1.8561007700569e-08 Lin (038) step:  1.9941348560423e+00 func: -1.1157992587561e+05
    Quad (005) norm:  1.8452144971921e-08 Lin (039) step:  1.3449506850978e+00 func: -1.1157992587561e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1721.60 °C, P =        0.1 MPa
    Liquid          moles:   0.423353 grams:  72.490
                 Fo form:  Mg2SiO4        X:  0.5160  wt%    SiO2   35.09
                 Fa form:  Fe2SiO4        X:  0.4840  wt%     FeO   40.62
                                                      wt%     MgO   24.29
    Olivine         moles:   0.576647 grams:  83.975
                 Fo form:  Mg2SiO4        X:  0.9218  wt%    SiO2   41.26
                 Fa form:  Fe2SiO4        X:  0.0782  wt%     FeO    7.72
                                                      wt%     MgO   51.03
    Add: Olivine
    Quad (000) norm:  1.5966604996323e+00 Lin (021) step:  4.0896682151642e-01 func: -1.1229597492575e+05
    Quad (001) norm:  3.3531760786162e-02 Lin (038) step:  1.9999999695048e+00 func: -1.1251755736637e+05
    Quad (002) norm:  1.8730262736856e-02 Lin (021) step:  9.5145608560682e-01 func: -1.1252036721908e+05
    Quad (003) norm:  4.1909407472519e-05 Lin (026) step:  9.9978900144526e-01 func: -1.1252036728945e+05
    Quad (004) norm:  9.2042603294142e-09 Lin (040) step: -4.4760217220127e-01 func: -1.1252036728945e+05
    Quad (005) norm:  1.3324106387581e-08 Lin (040) step:  4.7007485184832e-01 func: -1.1252036728945e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1687.95 °C, P =        0.1 MPa
    Liquid          moles:   0.440402 grams:  77.409
                 Fo form:  Mg2SiO4        X:  0.4440  wt%    SiO2   34.18
                 Fa form:  Fe2SiO4        X:  0.5560  wt%     FeO   45.46
                                                      wt%     MgO   20.36
    Olivine         moles:   0.559598 grams:  82.209
                 Fo form:  Mg2SiO4        X:  0.9015  wt%    SiO2   40.90
                 Fa form:  Fe2SiO4        X:  0.0985  wt%     FeO    9.64
                                                      wt%     MgO   49.47
    Add: Olivine
    Quad (000) norm:  1.4880582357278e+00 Lin (012) step:  4.1373492082155e-01 func: -1.1291818273809e+05
    Quad (001) norm:  3.5885647757584e-02 Lin (037) step:  1.9999999451233e+00 func: -1.1318463428273e+05
    Quad (002) norm:  2.0995382584375e-02 Lin (013) step:  9.4673678156252e-01 func: -1.1318859383769e+05
    Quad (003) norm:  1.2095598641376e-04 Lin (009) step:  1.0006930291438e+00 func: -1.1318859460434e+05
    Quad (004) norm:  1.0234833107969e-07 Lin (032) step:  9.9640056598608e-01 func: -1.1318859460434e+05
    Quad (005) norm:  3.6836719171515e-10 Lin (027) step: -1.0557282694161e+00 func: -1.1318859460434e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1654.30 °C, P =        0.1 MPa
    Liquid          moles:   0.458145 grams:  82.401
                 Fo form:  Mg2SiO4        X:  0.3792  wt%    SiO2   33.41
                 Fa form:  Fe2SiO4        X:  0.6208  wt%     FeO   49.60
                                                      wt%     MgO   16.99
    Olivine         moles:   0.541855 grams:  80.372
                 Fo form:  Mg2SiO4        X:  0.8790  wt%    SiO2   40.51
                 Fa form:  Fe2SiO4        X:  0.1210  wt%     FeO   11.72
                                                      wt%     MgO   47.77
    Add: Olivine
    Quad (000) norm:  1.3774563841416e+00 Lin (016) step:  4.1875957633954e-01 func: -1.1328816860212e+05
    Quad (001) norm:  3.8714527790115e-02 Lin (037) step:  1.9999999507498e+00 func: -1.1360092303644e+05
    Quad (002) norm:  2.3554754466832e-02 Lin (021) step:  9.4853093477755e-01 func: -1.1360670263280e+05
    Quad (003) norm:  3.8209968150830e-04 Lin (017) step:  1.0016221307589e+00 func: -1.1360671021612e+05
    Quad (004) norm:  9.3148899886895e-07 Lin (038) step:  9.8610154707078e-01 func: -1.1360671021613e+05
    Quad (005) norm:  1.2944080631323e-08 Lin (040) step:  4.6997624564121e-01 func: -1.1360671021613e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1620.65 °C, P =        0.1 MPa
    Liquid          moles:   0.476580 grams:  87.461
                 Fo form:  Mg2SiO4        X:  0.3211  wt%    SiO2   32.74
                 Fa form:  Fe2SiO4        X:  0.6789  wt%     FeO   53.15
                                                      wt%     MgO   14.11
    Olivine         moles:   0.523420 grams:  78.466
                 Fo form:  Mg2SiO4        X:  0.8539  wt%    SiO2   40.08
                 Fa form:  Fe2SiO4        X:  0.1461  wt%     FeO   14.00
                                                      wt%     MgO   45.92
    Add: Olivine
    Quad (000) norm:  1.2648704862387e+00 Lin (014) step:  4.2409981494523e-01 func: -1.1342194473354e+05
    Quad (001) norm:  4.1869525150838e-02 Lin (037) step:  1.9999999596761e+00 func: -1.1378169604333e+05
    Quad (002) norm:  2.6613942137230e-02 Lin (019) step:  9.5579384524493e-01 func: -1.1379024996267e+05
    Quad (003) norm:  7.9668730517902e-04 Lin (014) step:  1.0024585337351e+00 func: -1.1379027996411e+05
    Quad (004) norm:  2.9607950072807e-06 Lin (032) step:  1.0007908725777e+00 func: -1.1379027996423e+05
    Quad (005) norm:  2.3600833745774e-09 Lin (038) step:  1.7283434109660e+00 func: -1.1379027996423e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1587.00 °C, P =        0.1 MPa
    Liquid          moles:   0.495698 grams:  92.587
                 Fo form:  Mg2SiO4        X:  0.2694  wt%    SiO2   32.17
                 Fa form:  Fe2SiO4        X:  0.7306  wt%     FeO   56.20
                                                      wt%     MgO   11.63
    Olivine         moles:   0.504302 grams:  76.494
                 Fo form:  Mg2SiO4        X:  0.8258  wt%    SiO2   39.61
                 Fa form:  Fe2SiO4        X:  0.1742  wt%     FeO   16.50
                                                      wt%     MgO   43.89
    Add: Olivine
    Quad (000) norm:  1.1503493938843e+00 Lin (010) step:  4.2983715680788e-01 func: -1.1333111264439e+05
    Quad (001) norm:  4.5202800591801e-02 Lin (037) step:  1.9999999484570e+00 func: -1.1373749443677e+05
    Quad (002) norm:  3.0406944185969e-02 Lin (015) step:  9.6703569319832e-01 func: -1.1375011298796e+05
    Quad (003) norm:  1.4006732752022e-03 Lin (015) step:  1.0033804526085e+00 func: -1.1375019792047e+05
    Quad (004) norm:  6.4058759280046e-06 Lin (018) step:  9.9998379375210e-01 func: -1.1375019792108e+05
    Quad (005) norm:  7.0408383870557e-11 Lin (028) step:  1.7770876763913e+00 func: -1.1375019792108e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1553.35 °C, P =        0.1 MPa
    Liquid          moles:   0.515485 grams:  97.774
                 Fo form:  Mg2SiO4        X:  0.2236  wt%    SiO2   31.68
                 Fa form:  Fe2SiO4        X:  0.7764  wt%     FeO   58.82
                                                      wt%     MgO    9.50
    Olivine         moles:   0.484515 grams:  74.461
                 Fo form:  Mg2SiO4        X:  0.7941  wt%    SiO2   39.10
                 Fa form:  Fe2SiO4        X:  0.2059  wt%     FeO   19.25
                                                      wt%     MgO   41.65
    Add: Olivine
    Quad (000) norm:  1.0339899082375e+00 Lin (016) step:  4.3608564888796e-01 func: -1.1302399689041e+05
    Quad (001) norm:  4.8567089825116e-02 Lin (037) step:  1.9999999486317e+00 func: -1.1347515653818e+05
    Quad (002) norm:  3.5170758113792e-02 Lin (016) step:  9.8102269403651e-01 func: -1.1349350833412e+05
    Quad (003) norm:  2.2408633540299e-03 Lin (011) step:  1.0046605966100e+00 func: -1.1349371088887e+05
    Quad (004) norm:  1.1680331732484e-05 Lin (019) step:  9.9995395647818e-01 func: -1.1349371089098e+05
    Quad (005) norm:  4.5315777914587e-10 Lin (039) step:  1.4164078897546e+00 func: -1.1349371089098e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1519.70 °C, P =        0.1 MPa
    Liquid          moles:   0.535917 grams: 103.015
                 Fo form:  Mg2SiO4        X:  0.1832  wt%    SiO2   31.26
                 Fa form:  Fe2SiO4        X:  0.8168  wt%     FeO   61.06
                                                      wt%     MgO    7.68
    Olivine         moles:   0.464083 grams:  72.375
                 Fo form:  Mg2SiO4        X:  0.7581  wt%    SiO2   38.53
                 Fa form:  Fe2SiO4        X:  0.2419  wt%     FeO   22.29
                                                      wt%     MgO   39.19
    Add: Olivine
    Quad (000) norm:  9.1595877115216e-01 Lin (016) step:  4.4300754818810e-01 func: -1.1250631657962e+05
    Quad (001) norm:  5.1809328381885e-02 Lin (037) step:  1.9999999486317e+00 func: -1.1299834763691e+05
    Quad (002) norm:  4.1117061534466e-02 Lin (024) step:  9.9680566612166e-01 func: -1.1302449337895e+05
    Quad (003) norm:  3.4076848385436e-03 Lin (020) step:  1.0066232763024e+00 func: -1.1302493065333e+05
    Quad (004) norm:  2.0123925274523e-05 Lin (019) step:  9.9996011819525e-01 func: -1.1302493065956e+05
    Quad (005) norm:  8.9683820119860e-10 Lin (040) step:  1.2909649445544e+00 func: -1.1302493065956e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1486.05 °C, P =        0.1 MPa
    Liquid          moles:   0.556960 grams: 108.302
                 Fo form:  Mg2SiO4        X:  0.1478  wt%    SiO2   30.90
                 Fa form:  Fe2SiO4        X:  0.8522  wt%     FeO   62.97
                                                      wt%     MgO    6.13
    Olivine         moles:   0.443040 grams:  70.241
                 Fo form:  Mg2SiO4        X:  0.7170  wt%    SiO2   37.90
                 Fa form:  Fe2SiO4        X:  0.2830  wt%     FeO   25.65
                                                      wt%     MgO   36.46
    Add: Olivine
    Quad (000) norm:  7.9652687717774e-01 Lin (013) step:  4.5083929279971e-01 func: -1.1178154576498e+05
    Quad (001) norm:  5.4759990379579e-02 Lin (037) step:  1.9999999486317e+00 func: -1.1230772601811e+05
    Quad (002) norm:  4.8408414712455e-02 Lin (019) step:  1.0134693868266e+00 func: -1.1234406869686e+05
    Quad (003) norm:  5.0764193451215e-03 Lin (012) step:  1.0097330712646e+00 func: -1.1234495654607e+05
    Quad (004) norm:  3.5796951091354e-05 Lin (027) step:  9.9996152538313e-01 func: -1.1234495656488e+05
    Quad (005) norm:  2.3817893499710e-09 Lin (040) step:  1.9257054217066e-01 func: -1.1234495656488e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1452.40 °C, P =        0.1 MPa
    Liquid          moles:   0.578571 grams: 113.626
                 Fo form:  Mg2SiO4        X:  0.1171  wt%    SiO2   30.59
                 Fa form:  Fe2SiO4        X:  0.8829  wt%     FeO   64.60
                                                      wt%     MgO    4.81
    Olivine         moles:   0.421429 grams:  68.071
                 Fo form:  Mg2SiO4        X:  0.6698  wt%    SiO2   37.20
                 Fa form:  Fe2SiO4        X:  0.3302  wt%     FeO   29.38
                                                      wt%     MgO   33.42
    Add: Olivine
    Quad (000) norm:  6.7612492383195e-01 Lin (021) step:  4.5993569447173e-01 func: -1.1085100202990e+05
    Quad (001) norm:  5.7217682254019e-02 Lin (037) step:  1.9999999486317e+00 func: -1.1140075998402e+05
    Quad (002) norm:  5.7138211438038e-02 Lin (027) step:  1.0297538888655e+00 func: -1.1144986398854e+05
    Quad (003) norm:  7.5762134765027e-03 Lin (013) step:  1.0147727341214e+00 func: -1.1145159773969e+05
    Quad (004) norm:  7.1041544163877e-05 Lin (031) step:  9.9977449540369e-01 func: -1.1145159780779e+05
    Quad (005) norm:  1.7571977345191e-08 Lin (039) step:  1.1827281053379e+00 func: -1.1145159780779e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1418.75 °C, P =        0.1 MPa
    Liquid          moles:   0.600695 grams: 118.976
                 Fo form:  Mg2SiO4        X:  0.0906  wt%    SiO2   30.34
                 Fa form:  Fe2SiO4        X:  0.9094  wt%     FeO   65.98
                                                      wt%     MgO    3.69
    Olivine         moles:   0.399305 grams:  65.876
                 Fo form:  Mg2SiO4        X:  0.6151  wt%    SiO2   36.42
                 Fa form:  Fe2SiO4        X:  0.3849  wt%     FeO   33.53
                                                      wt%     MgO   30.05
    Add: Olivine
    Quad (000) norm:  5.5543949202481e-01 Lin (021) step:  4.7085053128253e-01 func: -1.0971359206561e+05
    Quad (001) norm:  5.8925356540378e-02 Lin (037) step:  1.9999999486317e+00 func: -1.1027108471228e+05
    Quad (002) norm:  6.7305478732473e-02 Lin (013) step:  1.0435208574435e+00 func: -1.1033523428082e+05
    Quad (003) norm:  1.1505030474154e-02 Lin (016) step:  1.0232335975700e+00 func: -1.1033853169246e+05
    Quad (004) norm:  1.6369600248230e-04 Lin (021) step:  9.9973911783214e-01 func: -1.1033853201481e+05
    Quad (005) norm:  5.5089535171400e-08 Lin (038) step:  1.0838690740601e+00 func: -1.1033853201481e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1385.10 °C, P =        0.1 MPa
    Liquid          moles:   0.623264 grams: 124.337
                 Fo form:  Mg2SiO4        X:  0.0679  wt%    SiO2   30.12
                 Fa form:  Fe2SiO4        X:  0.9321  wt%     FeO   67.14
                                                      wt%     MgO    2.74
    Olivine         moles:   0.376736 grams:  63.669
                 Fo form:  Mg2SiO4        X:  0.5513  wt%    SiO2   35.55
                 Fa form:  Fe2SiO4        X:  0.4487  wt%     FeO   38.15
                                                      wt%     MgO   26.29
    Add: Olivine
    Quad (000) norm:  4.3559312804386e-01 Lin (014) step:  4.8449627414159e-01 func: -1.0836493609787e+05
    Quad (001) norm:  5.9529476050007e-02 Lin (037) step:  1.9999999486317e+00 func: -1.0890705815873e+05
    Quad (002) norm:  7.8766975062179e-02 Lin (019) step:  1.0509975357788e+00 func: -1.0898730390469e+05
    Quad (003) norm:  1.7891152136413e-02 Lin (017) step:  1.0382321492869e+00 func: -1.0899343676869e+05
    Quad (004) norm:  4.2822065197050e-04 Lin (025) step:  9.9908284219637e-01 func: -1.0899343867446e+05
    Quad (005) norm:  4.1366423993629e-07 Lin (035) step:  1.0540567761337e+00 func: -1.0899343867446e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1351.45 °C, P =        0.1 MPa
    Liquid          moles:   0.646198 grams: 129.694
                 Fo form:  Mg2SiO4        X:  0.0487  wt%    SiO2   29.94
                 Fa form:  Fe2SiO4        X:  0.9513  wt%     FeO   68.11
                                                      wt%     MgO    1.96
    Olivine         moles:   0.353802 grams:  61.466
                 Fo form:  Mg2SiO4        X:  0.4763  wt%    SiO2   34.59
                 Fa form:  Fe2SiO4        X:  0.5237  wt%     FeO   43.31
                                                      wt%     MgO   22.10
    Add: Olivine
    Quad (000) norm:  3.1852437204954e-01 Lin (013) step:  5.0249398225183e-01 func: -1.0679500499436e+05
    Quad (001) norm:  5.8501068165817e-02 Lin (037) step:  1.9999999486317e+00 func: -1.0728856655177e+05
    Quad (002) norm:  9.1131195096223e-02 Lin (019) step:  1.0458731597712e+00 func: -1.0738276874258e+05
    Quad (003) norm:  2.8316580022175e-02 Lin (018) step:  1.0666690482862e+00 func: -1.0739383048254e+05
    Quad (004) norm:  1.2011354132180e-03 Lin (020) step:  9.9659503150087e-01 func: -1.0739384283993e+05
    Quad (005) norm:  3.5716283239582e-06 Lin (026) step:  1.0005734051080e+00 func: -1.0739384284007e+05
    Quad (006) norm:  2.0085466993128e-09 Lin (040) step:  4.7216433831497e-01 func: -1.0739384284007e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1317.80 °C, P =        0.1 MPa
    Liquid          moles:   0.669403 grams: 135.030
                 Fo form:  Mg2SiO4        X:  0.0327  wt%    SiO2   29.79
                 Fa form:  Fe2SiO4        X:  0.9673  wt%     FeO   68.91
                                                      wt%     MgO    1.31
    Olivine         moles:   0.330597 grams:  59.285
                 Fo form:  Mg2SiO4        X:  0.3876  wt%    SiO2   33.51
                 Fa form:  Fe2SiO4        X:  0.6124  wt%     FeO   49.07
                                                      wt%     MgO   17.42
    Add: Olivine
    Quad (000) norm:  2.0794856357104e-01 Lin (021) step:  5.2806377702831e-01 func: -1.0498115451600e+05
    Quad (001) norm:  5.4960557582732e-02 Lin (037) step:  1.9999999486317e+00 func: -1.0537895573655e+05
    Quad (002) norm:  1.0349900326869e-01 Lin (016) step:  1.0188187654698e+00 func: -1.0547784639120e+05
    Quad (003) norm:  4.4619563918756e-02 Lin (013) step:  1.1247128399950e+00 func: -1.0549649928907e+05
    Quad (004) norm:  3.4487847296435e-03 Lin (022) step:  9.8752670798971e-01 func: -1.0549657516277e+05
    Quad (005) norm:  3.1391133892242e-05 Lin (016) step:  1.0000585534400e+00 func: -1.0549657516925e+05
    Quad (006) norm:  1.3416486172291e-09 Lin (035) step:  4.7207970565543e-01 func: -1.0549657516925e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1284.15 °C, P =        0.1 MPa
    Liquid          moles:   0.692772 grams: 140.323
                 Fo form:  Mg2SiO4        X:  0.0194  wt%    SiO2   29.66
                 Fa form:  Fe2SiO4        X:  0.9806  wt%     FeO   69.56
                                                      wt%     MgO    0.77
    Olivine         moles:   0.307228 grams:  57.145
                 Fo form:  Mg2SiO4        X:  0.2817  wt%    SiO2   32.30
                 Fa form:  Fe2SiO4        X:  0.7183  wt%     FeO   55.49
                                                      wt%     MgO   12.21
    Add: Olivine
    Quad (000) norm:  1.1254439960757e-01 Lin (013) step:  5.6886545119054e-01 func: -1.0286064747519e+05
    Quad (001) norm:  4.7212340551214e-02 Lin (037) step:  1.9999999470379e+00 func: -1.0309716844711e+05
    Quad (002) norm:  1.1367867550790e-01 Lin (023) step:  9.5893721444821e-01 func: -1.0317623263192e+05
    Quad (003) norm:  6.6929065305317e-02 Lin (014) step:  1.2552998124023e+00 func: -1.0320172077162e+05
    Quad (004) norm:  1.0133058640396e-02 Lin (011) step:  9.5627186609852e-01 func: -1.0320208266356e+05
    Quad (005) norm:  2.6970989681121e-04 Lin (022) step:  1.0004716759792e+00 func: -1.0320208286835e+05
    Quad (006) norm:  1.0937803239169e-07 Lin (039) step:  1.5470767854492e+00 func: -1.0320208286835e+05
    Minimal energy termination of quadratic loop.
    
     
    T =    1250.50 °C, P =        0.1 MPa
    Liquid          moles:   0.716183 grams: 145.553
                 Fo form:  Mg2SiO4        X:  0.0086  wt%    SiO2   29.56
                 Fa form:  Fe2SiO4        X:  0.9914  wt%     FeO   70.09
                                                      wt%     MgO    0.34
    Olivine         moles:   0.283817 grams:  55.070
                 Fo form:  Mg2SiO4        X:  0.1544  wt%    SiO2   30.97
                 Fa form:  Fe2SiO4        X:  0.8456  wt%     FeO   62.62
                                                      wt%     MgO    6.42


.. code:: ipython3

    plt.plot(xFoSol, tC, 'b-')
    plt.plot(xFoLiq, tC, 'r-')
    plt.ylabel('T °C')
    plt.xlabel('Mole fraction')
    plt.xlim(0.0, 1.0)




.. parsed-literal::

    (0.0, 1.0)




.. image:: 9-Olivine-loop_files/9-Olivine-loop_27_1.png


