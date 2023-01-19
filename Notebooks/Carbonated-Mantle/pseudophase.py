class PseudoPhase:
    def __init__(self, T, P, database='Stixrude'):
        self._eval_system_props(database)

    def _eval_system_props(self, database):
        phases = get_subsolidus_phases(database=database)
        phs_sym, endmem_ids, mu, elem_comps, sys_elems = system_energy_landscape(
            T, P, phases)
        # display(phs_sym, endmem_ids, mu, elem_comps, sys_elems)
        phs_sym_uniq, endmem_ids_uniq, mu_uniq, elem_comps_uniq = (
            prune_polymorphs(phs_sym, endmem_ids, mu, elem_comps))
        Nelems = len(sys_elems)
        Npts = mu_uniq.size

        return phases,


    def get_subsolidus_phases(self, database='Berman'):
        remove_phases = ['Liq','H2O']

        modelDB = model.Database(database)
        phases = modelDB.phases
        [phases.pop(phs) for phs in remove_phases]

        return phases

    def system_energy_landscape(self, T, P, phases, TOL=1e-3):
        elem_comps = []
        phs_sym = []
        endmem_ids = []
        mu = []
        for phsnm in phases:
            phs = phases[phsnm]

            elem_comp = phs.props['element_comp']
            abbrev = phs.abbrev
            endmem_num = phs.endmember_num
            iendmem_ids = list(np.arange(endmem_num))

            if phs.phase_type=='pure':
                nelem = np.sum(elem_comp)
                mu += [phs.gibbs_energy(T, P)/nelem]
                # print(nelem)
            else:
                nelem = np.sum(elem_comp,axis=1)
                # print(nelem)
                for i in iendmem_ids:
                    imol = np.eye(phs.endmember_num)[i]
                    mu += [phs.gibbs_energy(T, P, mol=imol,deriv={"dmol":1})[0,i]/nelem[i]]
                    # print(nelem[i])

            endmem_ids.extend(iendmem_ids)
            phs_sym.extend(list(np.tile(abbrev,endmem_num)))
            # print(elem_comp)

            elem_comps.extend(elem_comp)
            # print(elem_comp)
            # print(phs)

        elem_comps = np.vstack(elem_comps)

        natoms = np.sum(elem_comps,axis=1)
        elem_comps = elem_comps/natoms[:,np.newaxis]

        elem_mask = ~np.all(elem_comps<TOL, axis=0)

        elem_comps = elem_comps[:, elem_mask]
        mu = np.array(mu)
        endmem_ids = np.array(endmem_ids)

        sys_elems = core.chem.PERIODIC_ORDER[elem_mask]
        return phs_sym, endmem_ids, mu, elem_comps, sys_elems

    def prune_polymorphs(self, phs_sym, endmem_ids, mu, elem_comps, decimals=4):
        elem_round_comps = np.round(elem_comps, decimals=decimals)
            # Drop identical comps
        elem_comps_uniq = np.unique(elem_round_comps, axis=0)

        # uniq_num = elem_comps_uniq.shape[0]
        mu_uniq = []
        phs_sym_uniq = []
        endmem_ids_uniq = []
        for elem_comp in elem_comps_uniq:
            is_equiv_comp = np.all(elem_round_comps == elem_comp[np.newaxis,:], axis=1)
            equiv_ind = np.where(is_equiv_comp)[0]
            min_ind = equiv_ind[np.argsort(mu[equiv_ind])[0]]
            min_mu = mu[min_ind]
            assert np.all(min_mu <= mu[equiv_ind]), 'fail'

            mu_uniq.append(min_mu)
            phs_sym_uniq.append(phs_sym[min_ind])
            endmem_ids_uniq.append(endmem_ids[min_ind])

        mu_uniq = np.array(mu_uniq)
        phs_sym_uniq = np.array(phs_sym_uniq)
        elem_comps_uniq = np.array(elem_comps_uniq)

        return phs_sym_uniq, endmem_ids_uniq, mu_uniq, elem_comps_uniq


    def min_energy_assemblage(self, bulk_comp, comp, mu, TOLmu=10, TOL=1e-5):
        xy = np.hstack((comp, mu[:,np.newaxis]))
        yavg = np.mean(mu)
        xy_bulk = np.hstack((bulk_comp, yavg))

        wt0, rnorm0 = opt.nnls(xy.T, xy_bulk)
        # print('rnorm',rnorm0)


        def fun(mu, shift=0):
            xy_bulk[-1] = mu
            wt, rnorm = opt.nnls(xy.T, xy_bulk)
            return rnorm-shift


        delmu = .1
        if rnorm0==0:
            shift_dir = -1
            soln_found = True
        else:
            output = opt.minimize_scalar(fun, bounds=[np.min(mu), np.max(mu)])
            xy_bulk[-1] = output['x']
            wt0, rnorm0 = opt.nnls(xy.T, xy_bulk)
            shift_dir = -1

        mu_prev=xy_bulk[-1]
        rnorm=rnorm0

        while True:
            mu_prev = xy_bulk[-1]
            rnorm_prev = rnorm

            xy_bulk[-1] += shift_dir*delmu
            wt, rnorm = opt.nnls(xy.T, xy_bulk)
            delmu *= 2

            # print(shift_dir, rnorm)
            if ((shift_dir==+1)&(rnorm>rnorm_prev)) or ((shift_dir==-1)&(rnorm>0)):
                break


        fun_fit = lambda mu, TOL=TOL: fun(mu, shift=TOL)
        if rnorm > TOL:
            mu_bulk = opt.brentq(fun_fit, mu_prev, xy_bulk[-1], xtol=TOLmu)
            xy_bulk[-1] = mu_bulk
            wt, rnorm = opt.nnls(xy.T, xy_bulk)

        mu_bulk = xy_bulk[-1]
        wt_bulk = wt


        ind_assem = np.where(wt_bulk>0)[0]
        return wt_bulk, mu_bulk, ind_assem


    def eval_curv(self, comps, method, cross_term_inds):
        single_pt = False
        if comps.ndim==1:
            single_pt = True
            comps = comps[np.newaxis,:]

        if method=='quad':
            XiXj = comps[:, cross_term_inds[0]]*comps[:, cross_term_inds[1]]
            X2_sum = np.sum(XiXj,axis=1)
            curv_term = X2_sum
        elif method=='quad-full':
            XiXj = comps[:, cross_term_inds[0]]*comps[:, cross_term_inds[1]]
            curv_term = XiXj
        elif method=='xlogx':
            logX = np.log(comps)
            logX[comps==0] = 0
            XlogX = comps*logX
            # XlogX[comps==0] = 0
            XlogX_sum = np.sum(XlogX,axis=1)
            curv_term = XlogX_sum
        elif method=='none':
            curv_term = np.zeros((comps.shape[0],0))
        else:
            assert False, method + ' is not a valid method for eval_curv.'

        if single_pt:
            curv_term = curv_term[0]

        return curv_term

    def init_lstsq(self, comps, mu, curv_method, cross_term_inds, yscl=None):
        curv_term = eval_curv(comps, curv_method, cross_term_inds)
        if curv_term.ndim==1:
            curv_term = curv_term[:,np.newaxis]

        print(curv_term)

        xobs = np.hstack((comps, curv_term))
        if yscl is None:
            yexp_scl = np.floor(np.log10(np.max(mu)-np.min(mu)))
            yscl = 10**yexp_scl

        yobs = mu/yscl

        return xobs, yobs, yscl
