def humid(P, Pt, MMa, MMb):
    H = (P / (Pt - P)) * (MMa / MMb)
    return H

λ_water = lambda T: 2500 - 2.546 * T #in degC

def Pa(H, Pt, MMa, MMb, guessP):
    from scipy.optimize import fsolve
    def res(var):
        P = var
        res = H -  (P / (Pt - P)) * (MMa / MMb)
        return res
    guess = [guessP]
    P = fsolve(res, guessP)
    return P[0]

def Vh(MMa, MMb, H, Pt, T):
    R = 8.3145 #kJ/kmolK
    V = (1 / MMb + H / MMa) * ((R * T) / Pt)
    return V

def enthalpy(Cpa, Cpb, H, T):
    λ_w = λ(T)
    S = Cpb + Cpa * H
    Tref = 273.15 #K
    enth = S * (T - Tref) + λ_w * H
    return enth

def percH(pa, psat, T, Pt, MMa, MMb, dec):
    import numpy as np
    H = humid(pa, Pt, MMa, MMb)
    Ho = humid(psat, Pt, MMa, MMb)
    val = H / Ho
    perc = np.round(val, dec) * 100
    return perc