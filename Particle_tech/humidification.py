'''
This module can be used to perform typical Humidification calculations.
Functions available in this module:
- humid: Calculate the humidity using the partial pressure of the compound that can condense, the total pressure and the molar masses.
                                                H = (pa / (Pt - pa)) * (MMa / MMb)                              (1)
- Pa: Calculate the partial pressure of the compound using equation (1) and a inital guess pressure to solve for pressure.
- Vh: Calculate the specific humid volume V'.
                                                V = (1 / MMb + H / MMa) * ((R * T) / Pt)                        (2)
- λ_water: Calulate the latent heat of vaporization of water. Temperature is in degrees C
                                                λ_water = 2500 - 2.546 * T                                      (3)
- λ_tol: Calulate the latent heat of vaporization of toluene. Temperature is in degrees C
                                                λ_tol = 406.65 - 0.752 * T                                      (4)
- enthalpy: Calculate the enthalpy of the system. Current version only solves for water and toluene and Cp values needs to be inserted manually.
                                                enth = S * (T - Tref) + λ_w * H                                 (5)
- percH: Calculate the percentage humidity. The ratio of the humidty to the saturated humidity.
                                                %H = H / Ho                                                     (6)
- relH: Calculate the relative humidity. The ratio of the partial pressure of the vapour to the saturated humidity.
                                                Hr = pa/pa_sat                                                  (7)

Nomenclature:
Symbol      Name                                                    Units
- pa        Partial pressure of component that will condense        (kPa)
- Pt        Total pressure                                          (kPa)
- MMa       Molar mass of component that will condense              (kg/kmol)
- MMb       Molar mass of component that will condense              (kg/kmol)
- H         Humidity                                                (kg/kg)
- T         Temperature                                             (K)
- pa_sat    Saturated pressure of component that will condense      (kPa)
- dec       Decimal places                                          (-)

This module works well in combination with the Antoine module.
Created by Dewald Mare.

'''

def humid(pa, Pt, MMa, MMb):
    H = (pa / (Pt - pa)) * (MMa / MMb)
    return H

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

λ_water = lambda T: 2500 - 2.546 * T #in degC
λ_tol = lambda T: 406.65 - 0.752 * T #in degC

def enthalpy(name, Cpa, Cpb, H, T):
    if name == "Water":
        λ = λ_water(T - 273.15)
    elif name == "Toluene":
        λ = λ_tol(T - 273.15)
    S = Cpb + Cpa * H
    Tref = 273.15 #K
    enth = S * (T - Tref) + λ * H
    return enth

def percH(pa, psat, Pt, MMa, MMb, dec):
    import numpy as np
    H = humid(pa, Pt, MMa, MMb)
    Ho = humid(psat, Pt, MMa, MMb)
    val = H / Ho
    perc = np.round(val, dec) * 100
    return perc

def relH(pa, psat, dec):
    Hr = pa / psat
    perc = np.round(Hr, dec) * 100
    return perc