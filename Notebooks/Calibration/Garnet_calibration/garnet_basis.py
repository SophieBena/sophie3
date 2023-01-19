import numpy as np
import pandas as pd
data_dir="data/"

def read_exp_data(oxides, deficiency_thresh=0.1):
    grt_exp=pd.read_excel(data_dir+'grt_bearing_expts.xls',sheetname=None)
    garnet_exp=grt_exp['Garnet']
    col_names = ['Wt: '+col for col in oxides]
    molar_data_exp = calc_mole_compos(garnet_exp, oxides, col_names)
    molar_data_exp = data_deficiency_filter(molar_data_exp, oxides, deficiency_thresh=deficiency_thresh)
    return molar_data_exp  

def read_nat_data(oxides, deficiency_thresh=0.1):
    grt_nat=pd.read_excel(data_dir+'natural_garnet_data.xls',sheetname=None)
    garnet_nat=grt_nat['Data']
    mask_null_SiO2 = ~pd.isnull(garnet_nat['SiO2'])
    # garnet_nat.loc[mask_null_SiO2][oxides]
    col_names = oxides
    molar_data_nat = calc_mole_compos(garnet_nat.loc[mask_null_SiO2], oxides, col_names)
    molar_data_nat = data_deficiency_filter(molar_data_nat, oxides, deficiency_thresh=deficiency_thresh)
    return molar_data_nat

def data_deficiency_filter(molar_data, oxides, deficiency_thresh=0.1):
    molar_data = molar_data.copy()
    site_occ_table, output = garnet_site_occupancy(molar_data['mol_frac_oxides'], oxides, full_output=True)
    high_deficiency = np.abs(output['total_deficiency'])>deficiency_thresh

    for key in molar_data:
        molar_data[key].mask(high_deficiency, inplace=True)
        molar_data[key].dropna(inplace=True)
        
    return molar_data
    
    

def get_oxide_info(oxide_cols):
    oxide_info = pd.read_csv(data_dir+'oxide_data.csv')
    oxide_info = oxide_info.set_index("name").T
    #oxide_data.index.name = None
    molec_weights = oxide_info.loc[['molecular_weight'],oxide_cols].values[0]
    cat_stoic = oxide_info.loc[['cation_stoic'],oxide_cols].values[0]
    oxy_stoic = oxide_info.loc[['oxygen_stoic'],oxide_cols].values[0]
    cat_names = oxide_info.loc[['cation'],oxide_cols].values[0]
    output = (molec_weights, cat_stoic, oxy_stoic, cat_names)
    return output

def lstsq_endmember_frac(rel_molar_data, rel_grt_stoic_mol_frac):
    """
    Least squares projection into component space
    """
    
    try: # Convert pandas dataframe into numpy array
        rel_molar_data = rel_molar_data.values
    except: # data is already numpy array
        pass
    
    Ndata = rel_molar_data.shape[0]
    Nendmem = rel_grt_stoic_mol_frac.shape[0]
    rel_molar_resid = np.zeros(rel_molar_data.shape)
    endmem_frac = np.zeros((Ndata, Nendmem))
    
    
    for i, idata in enumerate(rel_molar_data):
        # print(idata)
        iendmem_soln = np.linalg.lstsq(rel_grt_stoic_mol_frac.T, idata)
        iendmem_frac = iendmem_soln[0]
        
        idata_proj = np.dot(iendmem_frac, rel_grt_stoic_mol_frac)
        idata_resid = idata - idata_proj
        # print(idata_resid)
        endmem_frac[i,:] = iendmem_frac
        rel_molar_resid[i,:] = idata_resid
        
    return endmem_frac, rel_molar_resid

def calc_mole_compos(comp_data, oxides, col_names): 
    """
    Calculate molar composition from oxide weight % measurements
    
    Parameters
    ----------
    comp_data: pandas DataFrame
        Compositional data in wt%, one analysis per row
    oxides: string list
        list of oxide names in desired order
    col_names: string list
        list of column names from comp_data corresponding to desired oxides
        
    Outputs
    -------
    molar_data: dictionary containing molar composition data with keys:
        `mol_frac_oxides` - oxide mol fraction for each analysis
        `mol_cation_per_oxy` - moles of cations on per oxygen basis
        `mol_oxide_per_oxy` - moles of oxide on per oxygen basis (e.g. Al 2/3 per oxygen)
        `oxide_data` - nicely reformatted wt% oxides for each experiment
        
    """
    
    molec_weights, cation_stoic, oxy_stoic, cation_names = get_oxide_info(oxides) 
    
    #Extracting pertinent data and clean-up
    oxide_data=comp_data[col_names].copy() 
    oxide_data.columns=oxides
    oxide_data.dropna(axis=0, how='all', inplace=True)
    oxide_data.fillna(value=0, inplace=True)
    
    
    #Molar calculation of oxide, cation, and oxygen
    mol_oxide = oxide_data/molec_weights
    mol_cation = mol_oxide*cation_stoic
    mol_cation.columns=cation_names
    mol_oxy = mol_oxide*oxy_stoic
   
    mol_frac_oxides = mol_oxide.div(mol_oxide.sum(axis=1), axis=0)
    mol_cation_per_oxy = mol_cation.div(mol_oxy.sum(axis=1), axis=0)
    mol_oxide_per_oxy = mol_cation_per_oxy/(cation_stoic/oxy_stoic)
    mol_oxide_per_oxy.columns= [oxide_str +'/O' for oxide_str in oxides]

    
    #Dictionary of output
    molar_data = {} 
    molar_data['mol_frac_oxides'] = pd.DataFrame(data=mol_frac_oxides, columns=oxides)
    molar_data['mol_oxide_per_oxy'] = mol_oxide_per_oxy
    molar_data['mol_cation_per_oxy'] = mol_cation_per_oxy
    # molar_data['mol_oxy_sum'] = pd.DataFrame(data=mol_oxy_sum, columns=['O'])
    molar_data['oxide_data'] = pd.DataFrame(data=oxide_data, columns=oxides)
    
    return molar_data

def mol_frac_oxide_to_cat_per_oxy(mol_frac_oxide, oxides):
    Ndim = mol_frac_oxide.ndim
    molec_weights, cation_stoic, oxy_stoic, cation_names = get_oxide_info(oxides)
    
    try: # Convert pandas dataframe into numpy array
        mol_frac_oxide = mol_frac_oxide.values
    except: # data is already numpy array
        pass
    
    
    if Ndim==1:
        cat_per_oxy = (mol_frac_oxide*cation_stoic)/np.sum(mol_frac_oxide*oxy_stoic)
    else:
        cat_per_oxy = np.zeros(mol_frac_oxide.shape)
        for ind, imol_frac_oxide in enumerate(mol_frac_oxide):
            icat_per_oxy = (imol_frac_oxide*cation_stoic)/np.sum(imol_frac_oxide*oxy_stoic)
            cat_per_oxy[ind,:] = icat_per_oxy
    
    # Ideally, should return pandas dataframe instead of numpy array
    # cat_per_oxy = pd.DataFrame(data=cat_per_oxy, columns = oxides)
            
    return cat_per_oxy

def garnet_site_occupancy(mol_frac_oxide, oxides, full_output=False):
    oxy_basis = 12
    Y_site_num = 2
    X_site_num = 3
    T_site_num = 3
    
    molec_weights, cation_stoic, oxy_stoic, cation_names = get_oxide_info(oxides)
    
    mol_cation_per_oxy = pd.DataFrame(mol_frac_oxide_to_cat_per_oxy(mol_frac_oxide, oxides), columns=cation_names)
    
    Ndata = mol_cation_per_oxy.shape[0]
    X_occ_types = ['X:Mg','X:Ca','X:Mn','X:Fe','X:vacancy']
    Y_occ_types = ['Y:Al','Y:Ti','Y:Fe','Y:Cr']
    T_occ_types = ['T:Si','T:Al']
    occ_columns = X_occ_types+Y_occ_types+T_occ_types
    
    site_occ_table = pd.DataFrame(data=np.zeros((Ndata,len(occ_columns))), columns=occ_columns)
    
    total_site_num = Y_site_num+X_site_num+T_site_num
    
    cation_pfu = oxy_basis*mol_cation_per_oxy
    
    total_deficiency = total_site_num-cation_pfu.sum(axis=1)
    site_occ_table['X:vacancy'].where(total_deficiency<0, total_deficiency, inplace=True)
    
    site_occ_table['T:Si']=cation_pfu['Si']
    T_deficiency = T_site_num-cation_pfu['Si']
    #Si_ismissing = tet_deficiency > 0
    
    #Si_excess = np.zeros(Ndata)
    #Si_excess[~Si_ismissing] = np.abs(tet_deficiency[~Si_ismissing])
    
    site_occ_table['T:Al'].where(T_deficiency<0, T_deficiency, inplace=True)

    site_occ_table['Y:Al'] = cation_pfu['Al']-site_occ_table['T:Al']
    Y_deficiency = cation_pfu['Cr']-site_occ_table['Y:Al']
    
    site_occ_table['X:Mg']=cation_pfu['Mg']
    site_occ_table['X:Ca']=cation_pfu['Ca']
    site_occ_table['X:Mn']=cation_pfu['Mn']
    site_occ_table['Y:Cr']=cation_pfu['Cr']
    site_occ_table['Y:Ti']=cation_pfu['Ti']
    X_other = site_occ_table['X:Mg'] + site_occ_table['X:Mn'] + site_occ_table['X:Ca'] + site_occ_table['X:vacancy']
    X_Fe_ideal = X_site_num - X_other
    
    site_occ_table['X:Fe'] = np.minimum(X_Fe_ideal, cation_pfu['Fe'])
    site_occ_table['Y:Fe'] = np.maximum(cation_pfu['Fe'] - X_Fe_ideal, 0)
    
    output = {}
    output['total_deficiency']=total_deficiency
    output['Y_deficiency']=Y_deficiency
    output['T_deficiency']=T_deficiency
    output['mol_cation_per_oxy'] = mol_cation_per_oxy
    
    if full_output:
        return site_occ_table, output
    else:
        return site_occ_table
    
    