import numpy as np
import scipy as sp
from scipy import spatial
from scipy import interpolate
import matplotlib.pyplot as plt


class Soln():
    _R = 8.3145   # J/mol/K

    def __init__(self, W, mult=1):
        self.W = W
        self.mult = mult

    def activity(self, X, T=300, P=1, endmem=1):
        # DEBUG FIX NEEDED, effect of multiplicity possibly wrong here
        # this is just a guess and needs verification
        R = self._R
        mult = self.mult

        RTloga_coef = self.RTloga_coef(X, T=T, P=P, endmem=endmem)
        loga_coef = 1/mult*RTloga_coef/(R*T)

        if endmem==1:
            activity = X*np.exp(loga_coef)
        elif endmem==2:
            activity = (1-X)*np.exp(loga_coef)
        else:
            assert False, 'endmem is not valid. endmem must be (1,2).'

        return activity

    def RTloga_coef(self, X, T=300, P=1, endmem=1):
        if endmem==1:
            RTloga_coef = self._calc_RTloga_coef(X,  T=T, P=P)
        elif endmem==2:
            RTloga_coef = self._calc_RTloga_coef(1-X, T=T, P=P)
        else:
            assert False, 'endmem is not valid. endmem must be (1,2).'

        return RTloga_coef

    def entropy_mix_ideal(self, X):
        S_mix_ideal = -(self._calc_XlogX(X)+self._calc_XlogX(1-X))
        return S_mix_ideal

    def gibbs_mix(self, X, T, P=1):
        R = self._R

        S_mix_ideal = self.entropy_mix_ideal(X)
        G_XS = self.gibbs_excess(X, T=T, P=P)
        G_mix_ideal = -R*T*S_mix_ideal
        G_mix = G_mix_ideal + G_XS
        return G_mix

    def solvus(self, T_grid, P=1, N=300, upsamp_fac=10, allow_interp=True):
        X = np.linspace(0,1,N)
        Nupsamp = int(np.ceil(N*upsamp_fac))

        X_solvus = []
        for T in T_grid:
            G_mix = self.gibbs_mix(X, T, P=P)
            X_gap_bnds = self._calc_solvus_bounds(X, G_mix, N=Nupsamp,
                                                  allow_interp=allow_interp)
            X_solvus.append(X_gap_bnds)

        X_solvus = np.vstack(X_solvus)
        return X_solvus

    def test_solvus_finder(self, T=300, P=1, N=300, upsamp_fac=10):
        W = self._calc_W(T=T, P=P)

        X = np.linspace(0,1,N)
        Nupsamp = int(np.ceil(N*upsamp_fac))

        G_mix = self.gibbs_mix(X, T, P=P)
        X_gap_bnds, hull = self._calc_solvus_bounds(
            X, G_mix, N=Nupsamp, output_hull=True)

        plt.figure()
        plt.plot(X, G_mix, 'k-')
        plt.plot(hull['X'], hull['G_mix'], 'r-')

    ###################

    def _calc_XlogX(self, X):
        mult = self.mult
        # could use function decorator here
        with np.errstate(divide='ignore', invalid='ignore'):
            XlogX = mult*X*np.log(X)

        XlogX[X==0] = 0
        return XlogX

    # def _calc_activity(self, X, T=300, P=1):
    #     R = self._R
    #     loga_coef = self.RTloga_coef(X, T=T, P=P)/(R*T)
    #     activity = X*np.exp(loga_coef)
    #     return activity

    def _get_hull_ind(X, G_mix, fac=None):
        if fac is None:
            fac = 100

        max_G = np.max(np.abs(G_mix))

        upper_hull = [0.5, max_G*fac]
        XG = np.vstack((X, G_mix)).T
        XG = np.vstack((XG, upper_hull))
        hull = spatial.ConvexHull(XG)

        ind_hull = list(hull.vertices)
        ind_hull.remove(len(X))
        ind_hull = np.sort(ind_hull)

        return ind_hull

    def _calc_solvus_bounds(self, X, G_mix, fac=None, N=1001,
                            output_hull=False, allow_interp=True):
        if allow_interp:
            f_hires = interpolate.interp1d(X, G_mix, kind='cubic')
            X_hires = np.linspace(X[0], X[-1], N)
            G_mix_hires = f_hires(X_hires)
        else:
            X_hires = X
            G_mix_hires = G_mix

        ind_hull = _get_hull_ind(X_hires, G_mix_hires, fac)


        if np.all(~(np.diff(ind_hull)>1)):
            X_gap_bnds = [np.nan, np.nan]
            X_hull = X_hires
            G_mix_hull = G_mix_hires

        else:
            # ind_gap = np.where(np.diff(ind_hull)>1)[0][0]
            # ind_gap_bnds = [ind_hull[ind_gap], ind_hull[ind_gap+1]]

            ind_gaps = np.where(np.diff(ind_hull)>1)[0]
            ind_gap_left = ind_gaps[0]
            ind_gap_right = ind_gaps[-1]
            ind_gap_bnds = [ind_hull[ind_gap_left], ind_hull[ind_gap_right+1]]


            X_gap_bnds = X_hires[ind_gap_bnds]

            X_hull = np.hstack((X_hires[:ind_gap_bnds[0]],
                                X_hires[ind_gap_bnds[1]:]))
            G_mix_hull = np.hstack((G_mix_hires[:ind_gap_bnds[0]],
                                    G_mix_hires[ind_gap_bnds[1]:]))

        if output_hull:
            output = {}
            output['X'] = X_hull
            output['G_mix'] = G_mix_hull

            return X_gap_bnds, output


        return X_gap_bnds

class RegularSoln(Soln):
    def __init__(self, Wh, Wv=0, Ws=0, mult=1):
        self.Wh = Wh
        self.Ws = Ws
        self.Wv = Wv
        self.mult = mult

    def _calc_RTloga_coef(self, X, T=300, P=1):
        W = self._calc_W(T=T, P=P)
        RTloga_coef = W*(1-X)**2
        return RTloga_coef


    def gibbs_excess(self, X, T=300, P=1):
        W = self._calc_W(T=T, P=P)
        G_XS = W*X*(1-X)
        return G_XS

    def _calc_W(self, T=300, P=1):
        Wh, Wv, Ws = self.Wh, self.Wv, self.Ws
        W = Wh + P*Wv - T*Ws
        return W

class SimpleRegularSoln(RegularSoln):
    def __init__(self, W, mult=1):
        self.W = W
        self.mult = mult

    def _calc_W(self, T=None, P=None):
        W = self.W
        return W

class AsymmSoln(RegularSoln):
    def __init__(self, Wh, Wv=0, Ws=0, alpha=1, mult=1):
        self.Wh = Wh
        self.Ws = Ws
        self.Wv = Wv
        self.alpha = alpha
        self.mult = mult

    def RTloga_coef(self, X, T=300, P=1, endmem=1):
        alpha = self.alpha
        if endmem==1:
            RTloga_coef = self._calc_RTloga_coef(X, alpha, T=T, P=P)
        elif endmem==2:
            RTloga_coef = self._calc_RTloga_coef(1-X, 1/alpha, T=T, P=P)
        else:
            assert False, 'endmem is not valid. endmem must be (1,2).'

        return RTloga_coef

    def gibbs_excess(self, X, T=300, P=1):

        loga_coef_1 = self.RTloga_coef(X, T=T, P=P, endmem=1)
        loga_coef_2 = self.RTloga_coef(X, T=T, P=P, endmem=2)

        G_XS = X*loga_coef_1 + (1-X)*loga_coef_2
        return G_XS

    def _calc_RTloga_coef(self, X, alpha, T=300, P=1):
        W = self._calc_W(T=T, P=P)

        W_asymm = W*2*1/(1+alpha)
        phi = X*1/(X*1+(1-X)*alpha)

        RTloga_coef = W_asymm*(1-phi)**2
        return RTloga_coef

#
def calc_activity(X, loga_coef):
    activity = X*np.exp(loga_coef)
    return activity

#
def _RTloga_coef_regular(X, W):
    RTloga_coef = W*(1-X)**2
    return RTloga_coef


def _RTloga_coef_vanlaar(X, W, alpha):
    # assumes alpha1=1

    # W_asymm = W*2*alpha1/(alpha1+alpha2)
    # phi1 = X*alpha1/(X*alpha1+(1-X)*alpha2)

    W_asymm = W*2*1/(1+alpha)
    phi = X*1/(X*1+(1-X)*alpha)
    RTloga_coef = W_asymm*(1-phi)**2
    return RTloga_coef

#
def calc_XlogX(X):
    XlogX = X*np.log(X)
    XlogX[X==0] = 0
    return XlogX

def calc_regular_soln(X, W):
    RTloga_coef = _RTloga_coef_regular(X, W)
    activity = calc_activity(X, RTloga_coef)
    RTloga_coef_2 = _RTloga_coef_regular(1-X, W)
    activity_2 = calc_activity(1-X, RTloga_coef_2)

    G_XS = W*X*(1-X)
    S_mix_ideal = -(calc_XlogX(X)+calc_XlogX(1-X))

    model = {}
    model['X'] = X
    model['W'] = W
    model['X1'] = X
    model['X2'] = 1-X
    model['loga_coef1'] = RTloga_coef
    model['loga_coef2'] = RTloga_coef_2
    model['a1'] = activity
    model['a2'] = activity_2

    model['G_XS'] = G_XS
    model['S_mix_ideal'] = S_mix_ideal

    return model

def calc_asymm_soln(X, W, alpha):
    RTloga_coef = _RTloga_coef_vanlaar(X, W, alpha)
    activity = calc_activity(X, RTloga_coef)

    RTloga_coef_2 = _RTloga_coef_vanlaar(1-X, W, 1/alpha)
    activity_2 = calc_activity(1-X, RTloga_coef_2)

    G_XS = X*RTloga_coef + (1-X)*RTloga_coef_2
    S_mix_ideal = -(calc_XlogX(X)+calc_XlogX(1-X))

    model = {}
    model['X'] = X
    model['W'] = W
    model['X1'] = X
    model['X2'] = 1-X
    model['loga_coef1'] = RTloga_coef
    model['loga_coef2'] = RTloga_coef_2
    model['a1'] = activity
    model['a2'] = activity_2

    model['G_XS'] = G_XS
    model['S_mix_ideal'] = S_mix_ideal

    return model


def _get_hull_ind(X, G_mix, fac=None):
    if fac is None:
        fac = 100

    max_G = np.max(np.abs(G_mix))

    upper_hull = [0.5, max_G*fac]
    XG = np.vstack((X, G_mix)).T
    XG = np.vstack((XG, upper_hull))
    hull = spatial.ConvexHull(XG)

    ind_hull = list(hull.vertices)
    ind_hull.remove(len(X))
    ind_hull = np.sort(ind_hull)

    return ind_hull

def get_solvus_bounds(X, G_mix, fac=None, N=1001, output_hull=False):
    f_hires = interpolate.interp1d(X, G_mix, kind='cubic')
    X_hires = np.linspace(X[0], X[-1], N)
    G_mix_hires = f_hires(X_hires)

    ind_hull = _get_hull_ind(X_hires, G_mix_hires, fac)


    if np.all(~(np.diff(ind_hull)>1)):
        X_gap_bnds = [np.nan, np.nan]
        X_hull = X_hires
        G_mix_hull = G_mix_hires

    else:
        ind_gap = np.where(np.diff(ind_hull)>1)[0][0]

        ind_gap_bnds = [ind_hull[ind_gap], ind_hull[ind_gap+1]]
        X_gap_bnds = X_hires[ind_gap_bnds]

        X_hull = np.hstack((X_hires[:ind_gap_bnds[0]],
                            X_hires[ind_gap_bnds[1]:]))
        G_mix_hull = np.hstack((G_mix_hires[:ind_gap_bnds[0]],
                                G_mix_hires[ind_gap_bnds[1]:]))

    if output_hull:
        output = {}
        output['X'] = X_hull
        output['G_mix'] = G_mix_hull

        return X_gap_bnds, output


    return X_gap_bnds

def calc_solvus(RT, mod, N=1001):
    X = mod['X']
    G_XS = mod['G_XS']
    S_mix_ideal = mod['S_mix_ideal']

    X_solvus = []
    for iRT in RT:
        G_mix_ideal = -iRT*S_mix_ideal
        G_mix = G_mix_ideal + G_XS

        X_gap_bnds = get_solvus_bounds(X, G_mix, N=N)
        X_solvus.append(X_gap_bnds)

    X_solvus = np.vstack(X_solvus)
    return X_solvus

def test_solvus_finder(mod, RT_W=.3):
    RT = mod['W']*RT_W
    X = mod['X']
    G_XS = mod['G_XS']
    S_mix_ideal = mod['S_mix_ideal']

    G_mix_ideal = -RT*S_mix_ideal
    G_mix = G_mix_ideal + G_XS
    X_gap_bnds, hull = get_solvus_bounds(X, G_mix, N=1001, output_hull=True)

    plt.figure()
    plt.plot(X, G_mix, 'k-')
    plt.plot(hull['X'], hull['G_mix'], 'r-')
