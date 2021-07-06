'''
This module can be used to perform typical calculations assossiated with Particle technology.
Classes available in the module:
    - Humidification
    - ParticleCharacterisation
    - ParticleMotion
    - Sedimentation
    - PackedBed
    - FluidisedBed
    - GasCyclones

This module works well in combination with the Antoine module.
Created by Dewald Mare.

'''
class Humidification:
    '''
    This class can be used to perform typical Humidification calculations.

    Functions available in this class:

    - humid: Calculate the humidity using the partial pressure of the compound that can condense, the total pressure and the molar masses.
                                                    H = (pa / (Pt - pa)) * (MMa / MMb)                              (1)
    - Pa: Calculate the partial pressure of the compound using equation (1) and a inital guess pressure to solve for pressure.
    - Vh: Calculate the specific humid volume V'.
                                                    V' = (1 / MMb + H / MMa) * ((R * T) / Pt)                       (2)
    - heatvap: Calculate the heat of vaporization of a compound. Compounds available:
        * Water
        * Toluene
    - enthalpy: Calculate the enthalpy of the system. Current version only solves for water and toluene and Cp values needs to be inserted manually.
                                                    enth = S * (T - Tref) + λ * H                                   (4)
      where
                                                    S = Cpa + Cpb * H                                               (5)
    - percH: Calculate the percentage humidity. The ratio of the humidty to the saturated humidity.
                                                    %H = H / Ho                                                     (6)
    - relH: Calculate the relative humidity. The ratio of the partial pressure of the vapour to the saturated humidity.
                                                    Hr = pa/pa_sat                                                  (7)

    Parameters
    ----------
    Symbol      Name                                                    Units
    - H         Humidity                                                (kg/kg)
    - pa        Partial pressure of component that will condense        (kPa)
    - Pt        Total pressure                                          (kPa)
    - MMa       Molar mass of component that will condense              (kg/kmol)
    - MMb       Molar mass of component that will not condense          (kg/kmol)
    - Vh        Humid volume                                            (m^3/kg)
    - R         Gas constant                                            (kJ/kmol K)
    - T         Temperature                                             (K)
    - enth      Enthalpy of system                                      (kJ/kg)
    - λ         Heat of vaporization                                    (kJ/kg)
    - Cpa       Heat capacity of compound that will condense            (kJ/kg K)
    - Cpb       Heat capacity of compound that will not condense        (kJ/kg K)
    - %H        Percentage humidity                                     (%)
    - Ho        Saturated humidity                                      (kg/kg)
    - Hr        Relative humidity                                       (-)
    - pa_sat    Saturated pressure of component that will condense      (kPa)
    - dec       Decimal places for answer                               (-)
    '''
    def humid(pa, Pt, MMa, MMb):
        '''
        Calculate the humidity using the partial pressure of the compound that can condense, the total pressure and the molar masses.
                                                H = (pa / (Pt - pa)) * (MMa / MMb)
        
        Parameters
        ----------
        - pa        Partial pressure of component that will condense        (kPa)
        - Pt        Total pressure                                          (kPa)
        - MMa       Molar mass of component that will condense              (kg/kmol)
        - MMb       Molar mass of component that will not condense          (kg/kmol)

        Returns
        -------
        - H         Humidity                                                (kg/kg)
        '''
        return (pa / (Pt - pa)) * (MMa / MMb)

    def Pa(H, Pt, MMa, MMb, guessP):
        '''
        Calculate the partial pressure of the compound using equation (1) and a inital guess pressure to solve for pressure.

        Parameters
        ----------
        - H         Humidity                                                (kg/kg)
        - Pt        Total pressure                                          (kPa)
        - MMa       Molar mass of component that will condense              (kg/kmol)
        - MMb       Molar mass of component that will not condense          (kg/kmol)
        - guessP    Initial guess pressure                                  (kPa)

        Returns
        -------
        - Pa        Partial pressure of component that will condense        (kPa)        
        '''
        from scipy.optimize import fsolve
        def res(P):
            res = H -  (P / (Pt - P)) * (MMa / MMb)
            return res
        return fsolve(res, [guessP])[0]

    def Vh(MMa, MMb, H, Pt, T):
        '''
        Calculate the specific humid volume V'.
                                                V = (1 / MMb + H / MMa) * ((R * T) / Pt)
        
        Parameters
        ----------
        - MMa       Molar mass of component that will condense              (kg/kmol)
        - MMb       Molar mass of component that will not condense          (kg/kmol)
        - H         Humidity                                                (kg/kg)
        - Pt        Total pressure                                          (kPa)
        - T         Temperature                                             (K)

        Returns
        -------
        - V'        Humid volume                                            (m^3/kg)        
        '''
        R = 8.3145 #kJ/kmolK
        return (1 / MMb + H / MMa) * ((R * T) / Pt)

    def heatvap(name, T):
        '''
        Calculate the heat of vaporization of a compound. Compounds available:
        - Water
        - Toluene

        Parameters
        ----------
        - name      Name of compound that will condense                     (-)
        - T         Temperature                                             (K)

        Returns
        -------
        - λ         Heat of vaporization of compound                        (kJ/kg)
        '''
        if name == "Water":
            return 2500 - 2.546 * (T - 273.15)
        elif name == "Toluene":
            return 406.65 - 0.752 * (T - 273.15)

    def enthalpy(name, Cpa, Cpb, H, T):
        '''
        Calculate the enthalpy of the system. Current version only solves for water and toluene and Cp values needs to be inserted manually.
                                                enth = S * (T - Tref) + λ * H

        Parameters
        ----------
        - name      Name of compound that will condense                     (-)
        - Cpa       Heat capacity of compound that will condense            (kJ/kg K)
        - Cpb       Heat capacity of compound that will not condense        (kJ/kg K)
        - H         Humidity                                                (kg/kg)
        - λ         Heat of vaporization of compound                        (kJ/kg)
        - T         Temperature                                             (K)

        Returns
        -------
        - enth      Enthalpy of system                                      (kJ/kg)
        '''
        if name == "Water":
            λ = 2500 - 2.546 * (T - 273.15)
        elif name == "Toluene":
            λ = 406.65 - 0.752 * (T - 273.15)
        S = Cpb + Cpa * H
        print('S = ', S, 'kJ/kgK')
        Tref = 273.15 #K
       
        return  S * (T - Tref) + λ * H

    def percH(pa, psat, Pt, MMa, MMb, dec):
        '''
        Calculate the percentage humidity. The ratio of the humidty to the saturated humidity.
                                                %H = H / Ho

        Parameters
        ----------
        - pa        Partial pressure of component that will condense        (kPa)
        - pa_sat    Saturated pressure of component that will condense      (kPa)
        - Pt        Total pressure                                          (kPa)
        - MMa       Molar mass of component that will condense              (kg/kmol)
        - MMb       Molar mass of component that will not condense          (kg/kmol)
        - dec       Decimal places for answer                               (-)

        Returns
        -------
        - %H        Percentage humidity                                     (%)
        '''
        import numpy as np
        H = (pa / (Pt - pa)) * (MMa / MMb)
        Ho = (psat / (Pt - pa)) * (MMa / MMb)
        val = H / Ho
        perc = np.round(val, dec) * 100
        return perc

    def relH(pa, psat, dec):
        '''
        Calculate the relative humidity. The ratio of the partial pressure of the vapour to the saturated humidity.
                                                Hr = pa/pa_sat

        Parameters
        ----------
        - pa        Partial pressure of component that will condense        (kPa)
        - pa_sat    Saturated pressure of component that will condense      (kPa)

        Returns
        -------
        - Hr        Relative humidity                                       (-)
        '''
        import numpy as np
        return np.round(pa / psat, dec) * 100

class ParticleCharacterisation:
    ''' 
    This class can be used to perform typical particle characterisation calculations.

    Functions available in this class:

    - dv: Calculate the equivalent volume diameter of a particle:
                                            dv = (6 * Vp / np.pi) ** (1 / 3)                            (1)
    - ds: Calculate the equivalent surface area diameter of a particle:
                                            ds = (Ap / np.pi) ** (1 / 2)                                (2)
    - dsv: Calculate the Sauter mean diameter of a particle:
                                            dsv = (6 * Vp / Ap)                                         (3)
    - sphericity: Calculate the sphericity of a particle:
                                            sph = (dv / ds) ** 2                                        (4)

    Parameters
    ----------
    Symbol          Name                                    Units
    - dv              Equivalent volume diameter              (m)
    - ds              Equivalent surface area diameter        (m)
    - dsv             Sauter mean diameter                    (m)
    - sph             Sphericity                              (-)
    '''
    def dv(Vp):
        import numpy as np
        '''
        Calculate the equivalent volume diameter of a particle:
                                            dv = (6 * Vp / np.pi) ** (1 / 3)
        
        Parameters
        ----------
        - Vp              Volume of particle                      (m^3)

        Returns
        -------
        - dv              Equivalent volume diameter              (m)
        '''
        return (6 * Vp / np.pi) ** (1 / 3)

    def ds(Ap):
        import numpy as np
        '''
        Calculate the equivalent surface area diameter of a particle:
                                            ds = (Ap / np.pi) ** (1 / 2)
        
        Parameters
        ----------
        - Ap              Area of particle                        (m^2)

        Returns
        -------
        - ds              Equivalent surface area diameter        (m)
        '''
        return (Ap / np.pi) ** (1 / 2)

    def dsv(Vp, Ap):
        '''
        Calculate the Sauter mean diameter of a particle:
                                            dsv = (6 * Vp / Ap)
        
        Parameters
        ----------
        - Vp              Volume of particle                      (m^3)
        - Ap              Area of particle                        (m^2)

        Returns
        -------
        - dsv             Sauter mean diameter                    (m)`
        '''
        return 6 * Vp / Ap

    def sphericity(dv, ds):
        '''
        Calculate the sphericity of a particle:
                                        sph = (dv / ds) ** 2
        Parameters
        ----------
        - dv              Equivalent volume diameter              (m)
        - ds              Equivalent surface area diameter        (m)

        Returns
        -------
        - sph             Sphericity                              (-)`
        '''
        return (dv / ds) ** 2

class ParticleMotion:
    '''
    This class can be used to perform typical particle motion calculations.

    Functions available in this class:
    - Ar: Calculate the Archimedes number:
                                        Ar = (g * ρ * (ρs - ρ) * d ** 3) / (μ ** 2)                 (1)
    - Re: Calculate the Reynolds number:
                                        Re = (d * u * ρ) / μ                                        (2)
    - Ar_Re: Use the Schiller and Neuman equation to find the Achimedes number, given Re:
                                        Ar = 18 * Re + 2.7 * Re ** 1.68                             (3)
    - Re_Ar: Use the Schiller and Neuman equation to implicitly find the Reynolds number, given Ar:
                                        Ar = 18 * Re + 2.7 * Re ** 1.68                             (4)
    - u_laminar: Calculate the linear velocity of a particle in the laminar region:
                                        u = (g * (ρs - ρ) * d ** 2) / (18 * μ)                      (5)

    Parameters
    ---------
    Symbol                    Name                        Units
    - Ar                      Archimedes number           (-)
    - Re                      Reynolds number             (-)
    - g                       Gravity constant            (m/s^2)
    - ρ                       Fluid density               (kg/m^3)
    - ρs                      Solid density               (kg/m^3)
    - d                       Diameter                    (m)
    - μ                       Viscosity                   (kg/ms)
    - u                       Linear velocity             (m/s)
    '''
    def Ar(ρ, ρs, d, μ):
        '''
        Calculate the Archimedes number:
                                        Ar = (g * ρ * (ρs - ρ) * d ** 3) / (μ ** 2)
        
        Parameters
        ---------
        - g                       Gravity constant            (m/s^2)
        - ρ                       Fluid density               (kg/m^3)
        - ρs                      Solid density               (kg/m^3)
        - d                       Diameter                    (m)
        - μ                       Viscosity                   (kg/ms)

        Returns
        -------
        - Ar                      Archimedes number           (-)                
        '''
        g = 9.81 #m/s2
        return (g * ρ * (ρs - ρ) * d ** 3) / (μ ** 2)

    def Re(d, u, ρ, μ):
        '''
        Calculate the Reynolds number:
                                        Re = (d * u * ρ) / μ
        
        Parameters
        ---------
        - d                       Diameter                    (m)
        - u                       Linear velocity             (m/s)
        - ρ                       Fluid density               (kg/m^3)
        - μ                       Viscosity                   (kg/ms)

        Returns
        -------
        - Re                      Reynolds number             (-)                
        '''
        return (d * u * ρ) / μ

    def Ar_Re(Re):
        '''
        Use the Schiller and Neuman equation to find the Achimedes number, given Re:
                                        Ar = 18 * Re + 2.7 * Re ** 1.68
        
        Parameters
        ---------
        - Re                       Reynolds number            (-)

        Returns
        -------
        - Ar                       Archimedes number          (-)

        '''
        return 18 * Re + 2.7 * Re ** 1.68

    def Re_Ar(Ar, guess):
        '''
        Use the Schiller and Neuman equation to implicitly find the Reynolds number, given Ar:
                                        Ar = 18 * Re + 2.7 * Re ** 1.68
        
        Parameters
        ---------
        - Ar                      Archimedes number           (-)

        Returns
        -------
        - Re                      Reynolds number             (-)

        '''
        from scipy.optimize import fsolve
        def sol(Re):
            res = Ar - 18 * Re - 2.7 * Re ** 1.68
            return res
        return fsolve(sol, [guess])[0]

    def u_laminar(g, ρ, ρs, d, μ):
        '''
        Calculate the linear velocity of a particle in the laminar region:
                                        u = (g * (ρs - ρ) * d ** 2) / (18 * μ)

        Parameters
        ---------
        - g                       Gravity constant            (m/s^2)
        - ρ                       Fluid density               (kg/m^3)
        - ρs                      Solid density               (kg/m^3)
        - d                       Diameter                    (m)
        - μ                       Viscosity                   (kg/ms)

        Returns
        -------
        - u                      Linear velocity              (m/s)
        '''
        g = 9.81 #m/s2
        return (g * (ρs - ρ) * d ** 2) / (18 * μ)

class Sedimentation:
    '''
    This class can be used to perform typical Sedimentation calculations.

    Functions available in this class:
    - uH: Calculate hindered settling velocity:
                                    uH = uS * ϵ ** n                                                            (1)
    - n: Implicitly solve for the constant relating hindered settling velocity to single-particle velocity:
                                    (4.8 - n) / (n - 2.4) = 0.043 * Ar ** 0.57 * (1 - 2.4 * (d / D) ** 0.27)    (2)
    - Amin: Calculate the minimum cross-sectional area of the thickner:
                                    Amin = Qf * Cf / Gfmax                                                      (3)
    - uo: Calculate overflow velocity, stating if it is larger or smaller than feed velocity, assuming Co = 0 kg/m^3:
                                    uo = Qo / A                                                                 (4)
    Parameters
    ----------
    Symbol          Name                            Units
    - uH            Hindered setteling velocity     (m/s)
    - uS            Single-particle velocity        (m/s)
    - ϵ             Voidage                         (-)
    - n             Relational constant             (-)
    - Ar            Archeimedes number              (-)
    - d             Particle diameter               (m)
    - D             Vessel diameter                 (m)
    - Amin          Minimum cross-sectional area    (m^2)
    - Qf            Feed volumetric flowrate        (m^3/s)
    - Cf            Feed concentration              (kg/m^3)
    - Gfmax         Feed flux                       (kg/m^2s)
    - uo            Linear overflow velocity        (m/s)
    - Qo            Overflow volumetric flowrate    (m^3/s)
    - Co            Overflow concentration          (kg/m^3)
    '''
    def uH(uS, ϵ, n):
        '''
        Calculate hindered settling velocity:
                                    uH = uS * ϵ ** n
        Parameters
        ----------
        - uS            Single-particle velocity        (m/s)
        - ϵ             Voidage                         (-)
        - n             Relational constant             (-)

        Returns
        ------
        - uH            Hindered setteling velocity     (m/s)
        '''
        return uS * ϵ ** n

    def n(Ar, d, Neglect_D):
        '''
        Implicitly solve for the constant relating hindered settling velocity to single-particle velocity:
                                    (4.8 - n) / (n - 2.4) = 0.043 * Ar ** 0.57 * (1 - 2.4 * (d / D) ** 0.27)
        
        Parameters
        ----------
        - Ar            Archeimedes number              (-)
        - d             Particle diameter               (m)
        - D             Vessel diameter                 (m)

        Returns
        ------
        - n             Relational constant             (-)
        '''
        if Neglect_D == True:
            return (4.8 + 2.4 * (0.043 * Ar ** 0.57)) / (0.043 * Ar ** 0.57 + 1)
        elif Neglect_D == False:
            D = float(input('What is the vessel diameter in m: '))
            return (4.8 + 2.4 * (0.043 * Ar ** 0.57 * (1 - 2.4 * (d / D) ** 0.27))) / ((0.043 * Ar ** 0.57 * (1 - 2.4 * (d / D) ** 0.27)) + 1)

    def Amin(Qf, Cf, Gfmax):
        '''
        Calculate the minimum cross-sectional area of the thickner:
                                    Amin = Qf * Cf / Gfmax

        Parameters
        ----------
        - Qf            Feed volumetric flowrate        (m^3/s)
        - Cf            Feed concentration              (kg/m^3)
        - Gfmax         Feed flux                       (kg/m^2s)

        Returns
        -------
        - Amin          Minimum cross-sectional area    (m^2)
        '''
        return Qf * Cf / Gfmax

    def uo(Qf, Cf, Cu, Gf, SF, Amin):
        '''  
        Calculate overflow velocity, stating if it is larger or smaller than feed velocity, assuming Co = 0 kg/m^3:
                                    uo = Qo / A
        
        Parameters
        ----------
        - SF            Safety factor                   (-)
        - Amin          Minimum cross-sectional area    (m^2)
        - Qf            Feed volumetric flowrate        (m^3/s)
        - Cf            Feed concentration              (kg/m^3)
        - Gfmax         Feed flux                       (kg/m^2s)
        - uo            linear velocity of overflow     (m/s)
        - Qo            Overflow volumetric flowrate    (m^3/s)
        - Co            Overflow concentration          (kg/m^3)

        Returns
        -------
        Either,
        - uo            Linear overflow velocity        (m/s)
        - usf           Linear feed settling velocity   (m/s)
        If usf < uo:
        - Act           Actual area of thickner         (m^2)
        '''
        #Assume Co = 0:
        Qu = Qf * Cf / Cu
        Qo = Qf - Qu
        uo = Qo / (SF * Amin)
        usf = Gf / Cf
        print(f'usf = {usf} m/s')
        print(f'uo = {uo} m/s')
        if usf < uo:
            print('Feed settling velocity is less than the velocity of the overflow, thus thickner will overflow, use usf to calculate Amin')
            Amin = Qo / usf
            return Amin * SF
        elif usf > uo:
            print('Feed settling velocity is greater than the velocity of the overflow, thus thickner will not overflow and Amin can be used')

class PackedBed:
    ''' 
    This class can be used to perform typical Packed bed calculations.

    Functions available in this class:
    - Ergun: Calculate the pressure drop across a packed bed reactor:
                    ΔP / L = 150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d)    (1)
    - Ergun_d: Calculate the diameter of the particles in a packed bed reactor given the pressure drop using the Ergun equation:
                    ΔP / L = 150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d)    (2)
    - Ergun_u: Calculate the velocity of the particles in a packed bed reactor given the pressure drop using the Ergun equation:
                    ΔP / L = 150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d)    (3)
    
    Parameters
    ----------
    Symbol                  Name                        Units
    - ΔP                    Pressure drop               (Pa)
    - L                     Length/Height of reactor    (m)
    - ϵ                     Porosity                    (-)
    - μ                     Viscosity                   (kg/ms)
    - uo                    Linear velocity             (m/s)
    - d                     Diameter of particles       (m)
    - ρ                     Density of fluid            (kg/m^3)
    '''
    def Ergun(L, ϵ, μ, uo, d, ρ):
        '''
        Calculate the pressure drop across a packed bed reactor:
                    ΔP / L = 150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d)
    
        Parameters
        ----------
        - L                     Length/Height of reactor    (m)
        - ϵ                     Porosity                    (-)
        - μ                     Viscosity                   (kg/ms)
        - uo                    Linear velocity             (m/s)
        - d                     Diameter of particles       (m)
        - ρ                     Density of fluid            (kg/m^3)

        Returns
        -------
        - ΔP                    Pressure drop               (Pa)
        '''
        return L * (150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d))

    def Ergun_d(L, ϵ, μ, uo, ΔP, ρ, guess):
        '''
        Calculate the diameter of the particles in a packed bed reactor given the pressure drop using the Ergun equation:
                    ΔP / L = 150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d)
    
        Parameters
        ----------
        - L                     Length/Height of reactor    (m)
        - ϵ                     Porosity                    (-)
        - μ                     Viscosity                   (kg/ms)
        - uo                    Linear velocity             (m/s)
        - ΔP                    Pressure drop               (Pa)
        - ρ                     Density of fluid            (kg/m^3)

        Returns
        -------
        - d                     Diameter of particles       (m)
        '''
        from scipy.optimize import fsolve
        def sol(d):
            res = ΔP / L - (150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d))
            return res
        return fsolve(sol, [guess])[0]
    
    def Ergun_u(L, ϵ, μ, d, ΔP, ρ, guess):
        '''
        Calculate the velocity of the particles in a packed bed reactor given the pressure drop using the Ergun equation:
                    ΔP / L = 150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d)
    
        Parameters
        ----------
        - L                     Length/Height of reactor    (m)
        - ϵ                     Porosity                    (-)
        - μ                     Viscosity                   (kg/ms)
        - d                     Diameter of particles       (m)
        - ΔP                    Pressure drop               (Pa)
        - ρ                     Density of fluid            (kg/m^3)

        Returns
        -------
        - uo                    Linear velocity             (m/s)
        '''
        from scipy.optimize import fsolve
        def sol(uo):
            res = ΔP / L - (150 * ((1 - ϵ) ** 2 / ϵ ** 3) * (μ * uo / d ** 2) + 1.75 * ((1 - ϵ) / ϵ ** 3) * (ρ * uo ** 2 / d))
            return res
        return fsolve(sol, [guess])[0]
    
class FluidisedBed:
    '''
    This class can be used to perform typical fluidized bed calculations.

    Functions available in this class:
    - delP: Calculate the pressure drop across a fluidized bed using a force balance:
                                    ΔP / L = (1 - ϵ) * (ρs - ρ) * g                                                     (1)
    - Find_ϵ: Find ϵ by solving a combined Ergun and force balance equation:
                                    Ar = 150 * ((1 - ϵ) / ϵ ** 3) * Re + 1.75 * (1 / ϵ ** 3) * Re ** 2                  (2)
    - Find_Re: Find Re by solving a combined Ergun and force balance equation:
                                    Ar = 150 * ((1 - ϵ) / ϵ ** 3) * Re + 1.75 * (1 / ϵ ** 3) * Re ** 2                  (3) 
    - umf_Baeyens: Calculate minimum fluidisation velocity using Baeyens equation:
                                    umf = ((ρs - ρ) ** 0.934 g ** 0.934 * dp ** 1.8) / (1110 * μ ** 0.87 * ρ ** 0.066)  (4)
    - umf_WenYu: Calculate minimum fluidisation velocity using Wen&Yu equation:
                                    umf = (d * u * ρ) / (Remf * μ)                                                      (5)
    - Fr: Calculate the type of fluidisation using Wilhelm and Kwauk correlation:
                                    Fr = umf ** 2 / (dp * g)                                                            (6)
    - umb: Calculate the superficial gas velocity where bubbles start to appear:
                                    umb = 2.07 * np.exp(0.716 * F) * ((dp * ρs ** 0.06) / μ ** 0.37)                    (7)

    Parameters
    ----------
    Symbol                  Name                            Units
    - ΔP                    Pressure drop                   (Pa)
    - L                     Length/Height of reactor        (m)
    - ϵ                     Porosity                        (-)
    - μ                     Viscosity                       (kg/ms)
    - uo                    Linear velocity                 (m/s)
    - dp                    Diameter of particles           (m)
    - ρ                     Density of fluid                (kg/m^3)
    - ρs                    Solid density of fluid          (kg/m^3)
    - Remf                  Reynolds number                 (-)
    - umf                   Minimum fluidisation velocity   (m/s)
    - g                     Gravitation constant            (m/s^2)
    - Fr                    Type of fluidisation            (-)
    - umb                   Bubbling velocity               (m/s)
    - F                     Fraction of powder < 45 μm      (-)
    '''
    def delP(L, ϵ, ρ, ρs):
        '''
        Calculate the pressure drop across a fluidized bed using a force balance:
                                    ΔP / L = (1 - ϵ) * (ρs - ρ) * g
        
        Parameters
        ----------
        - L                     Length/Height of reactor        (m)
        - ϵ                     Porosity                        (-)
        - ρ                     Density of fluid                (kg/m^3)
        - ρs                    Solid density of fluid          (kg/m^3)
        - g                     Gravitation constant            (m/s^2)

        Returns
        -------
        - ΔP                    Pressure drop                   (Pa)
        '''
        g = 9.81 #m/s2
        return L * (1 - ϵ) * (ρs - ρ) * g

    def Find_ϵ(Ar, dp, u, ρ, μ, guess):
        '''
        Find ϵ by solving a combined Ergun and force balance equation:
                                    Ar = 150 * ((1 - ϵ) / ϵ ** 3) * Re + 1.75 * (1 / ϵ ** 3) * Re ** 2

        Parameters
        ----------
        - Ar                    Archimedes number               (-)
        - dp                    Diameter of particles           (m)
        - u                     Linear velocity                 (m/s)
        - ρ                     Density of fluid                (kg/m^3)
        - μ                     Viscosity                       (kg/ms)

        Returns
        -------
        - ϵ                     Porosity                        (-) 
        '''
        from scipy.optimize import fsolve
        Re = (d * u * ρ) / μ
        def sol(ϵ):
            res = Ar - 150 * ((1 - ϵ) / ϵ ** 3) * Re - 1.75 * (1 / ϵ ** 3) * Re ** 2
            return res
        return fsolve(sol, [guess])[0]


    def Find_Re(Ar, ϵ, guess):
        '''
        Find Re by solving a combined Ergun and force balance equation:
                                    Ar = 150 * ((1 - ϵ) / ϵ ** 3) * Re + 1.75 * (1 / ϵ ** 3) * Re ** 2

        Parameters
        ----------
        - Ar                    Archimedes number               (-)
        - ϵ                     Porosity                        (-)

        Returns
        -------
        - Re                    Reynolds number                 (-)
        '''
        from scipy.optimize import fsolve
        def sol(Re):
            res = Ar - 150 * ((1 - ϵ) / ϵ ** 3) * Re - 1.75 * (1 / ϵ ** 3) * Re ** 2
            return res
        return fsolve(sol, [guess])[0]

    def umf_Baeyens(ρs, ρ, dp, μ):
        '''
        Minimum fluidisation velocity using Baeyens equation:
                                    umf = ((ρs - ρ) ** 0.934 * g ** 0.934 * dp ** 1.8) / (1110 * μ ** 0.87 * ρ ** 0.066)

        Parameters
        ----------
        - ρ                     Density of fluid                (kg/m^3)
        - ρs                    Solid density of fluid          (kg/m^3)
        - g                     Gravitation constant            (m/s^2)
        - μ                     Viscosity                       (kg/ms)
        - dp                    Diameter of particles           (m)

        Returns
        -------
        - umf                   Minimum fluidisation velocity   (m/s)
        '''
        g = 9.81 #m/s2
        return ((ρs - ρ) ** 0.934 * g ** 0.934 * dp ** 1.8) / (1110 * μ ** 0.87 * ρ ** 0.066)

    def umf_WenYu(ρs, ρ, d, μ):
        '''
        Calculate minimum fluidisation velocity using Wen&Yu equation:
                                    umf = (d * u * ρ) / (Remf * μ)
        
        Parameters
        ----------
        - ρ                     Density of fluid                (kg/m^3)
        - ρs                    Solid density of fluid          (kg/m^3)
        - g                     Gravitation constant            (m/s^2)
        - μ                     Viscosity                       (kg/ms)
        - dp                    Diameter of particles           (m)

        Returns
        -------
        - umf                   Minimum fluidisation velocity   (m/s)
        '''
        g = 9.81 #m/s2
        Ar = (g * ρ * (ρs - ρ) * d ** 3) / (μ ** 2)
        print(f'Ar = {Ar}')
        Remf = 33.7 * ((1 + 3.59e-5 * Ar) ** (0.5) - 1)
        print(f'Remf = {Remf}')
        umf = (Remf * μ) / (d * ρ)
        return umf
    
    def Fr(umf, dp):
        '''
        Calculate the type of fluidisation using Wilhelm and Kwauk correlation:
                                    Fr = umf ** 2 / (dp * g)
        Parameters
        ----------
        - umf                   Minimum fluidisation velocity   (m/s)
        - dp                    Diameter of particles           (m)

        Returns
        -------
        - Fr                    Type of fluidisation            (-)
        '''
        g = 9.81 #m/s2
        Fr = umf ** 2 / (dp * g)

        if Fr < 1:
            print('Particulate (or homogeneous) fluidisation. The bed behaves uniformly with the particles, and the fluid is distributed evenly through the expanded bed.')
            return Fr
        elif Fr > 1:
            print('The bed expands uniformly, but gas bubbles are observed. These bubbles have a lower density than the solids mixture and tend to rise rapidly to the bed surface.')
            return Fr

    def umb(F, dp, ρs, μ):
        '''
        Calculate the superficial gas velocity where bubbles start to appear:
                                    umb = 2.07 * np.exp(0.716 * F) * ((dp * ρs ** 0.06) / μ ** 0.37)

        Parameters
        ----------
        - F                     Fraction of powder < 45 μm      (-)
        - dp                    Diameter of particles           (m)
        - ρs                    Solid density of fluid          (kg/m^3)
        - μ                     Viscosity                       (kg/ms)

        Returns
        -------
        - umb                   Bubbling velocity               (m/s)
        '''
        import numpy as np
        return 2.07 * np.exp(0.716 * F) * ((dp * ρs ** 0.06) / μ ** 0.37)

class GasCyclones:
    '''
    This class can be used to perform typical gas cyclone calculations.

    Functions available in this class:
    - velocity: Calculate the characteristic velocity:
                                    Eu = delP / (1/2 * ρ * v ** 2)                                                         (1)
    - diameter: Calculate the diameter of a gas cyclone(s):
                                    v = (Q / n) / (np.pi/4 * D ** 2)                                                       (2)
    - x50: Calculate the cut diameter:
                                    Stk50 = (x50 ** 2 * ρs * v) / (18 * μ * D)                                             (3)
    - n: Calculate number of cyclones with above functions.
    Parameters
    ----------
    Symbol                      Name                            Units
    - Eu                        Euler number                    (-)
    - delP                      Pressure drop                   (Pa)
    - ρ                        Density                         (kg/m^3)
    - ρs                        Solid density                   (kg/m^3)
    - μ                         Viscosity                       (kg/ms)
    - v                         Velocity                        (m/s)
    - Q                         Volumetric flow rate            (m^3/s)
    - n                         Number of cyclones              (-)
    - D                         Cyclone diameter                (m)
    - Stk50                     Stokes number                   (-)
    - x50                       Cut diameter                    (m)
    '''
    def velocity(Eu, delP, ρ):
        '''
        Calculate the characteristic velocity:
                                    Eu = delP / (1/2 * ρ * v ** 2)
        Parameters
        ----------
        - Eu                        Euler number                    (-)
        - delP                      Pressure drop                   (Pa)
        - ρ                         Density                         (kg/m^3)

        Returns
        -------
        - v                         Velocity                        (m/s)
        '''
        import numpy as np
        return np.sqrt(delP / (1/2 * ρ * Eu))

    def diameter(v, Q, n):
        '''
        Calculate the diameter of a gas cyclone(s):
                                    v = (Q / n) / (np.pi/4 * D ** 2)
        Parameters
        ----------
        - v                         Velocity                        (m/s)
        - Q                         Volumetric flow rate            (m^3/s)
        - n                         Number of cyclones              (-)

        Returns
        -------
        - D                         Cyclone diameter                (m)
        '''
        import numpy as np
        return np.sqrt((Q / n) / (np.pi / 4 * v))

    def x50(Stk50, ρs, v, μ, D):
        '''
        Calculate the cut diameter:
                                    Stk50 = (x50 ** 2 * ρs * v) / (18 * μ * D)
        Parameters
        ----------
        - Stk50                     Stokes number                   (-)
        - ρs                        Solid density                   (kg/m^3)
        - v                         Velocity                        (m/s)
        - μ                         Viscosity                       (kg/ms)
        - D                         Cyclone diameter                (m)

        Returns
        -------
        - x50                       Cut diameter                    (m)
        '''
        import numpy as np
        return np.sqrt((Stk50 * 18 * μ * D) / (ρs * v))

    def n(Q, v, Stk50, x50, ρs, μ, D, guess):
        '''
        Calculate the number of cyclones in parallel.
        Parameters
        ----------
        - Eu                        Euler number                    (-)
        - delP                      Pressure drop                   (Pa)
        - ρ                         Density                         (kg/m^3)
        - ρs                        Solid density                   (kg/m^3)
        - μ                         Viscosity                       (kg/ms)
        - v                         Velocity                        (m/s)
        - Q                         Volumetric flow rate            (m^3/s)
        - n                         Number of cyclones              (-)
        - D                         Cyclone diameter                (m)
        - Stk50                     Stokes number                   (-)
        - x50                       Cut diameter                    (m)
        Returns
        -------
        - n                         Number of cyclones              (-)
        '''
        from scipy.optimize import fsolve
        import numpy as np
        def sol(n):
            D = ((Q / n) / (np.pi / 4 * v)) ** (1/2)
            res = Stk50 - (x50 ** 2 * ρs * v) / (18 * μ * D)
            return res
        n = fsolve(sol, [guess])[0]
        D = np.sqrt((Q / n) / (np.pi / 4 * v))
        print('D = ', D, 'm')
        return n
