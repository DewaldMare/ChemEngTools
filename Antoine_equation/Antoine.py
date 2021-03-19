def Antoine_eqn(name, T):
    '''
    This  solver computes the vapour pressure of a compound at a given temperature. It uses the Antoine equation in the following form:
    ln(P*) = A - B/(T + C) with P* in kPa and T in K.
    Available compounds:
    - Acetic acid
    - Acetone 
    - Ammonia 
    - Benzene
    - Carbon tetrachloride
    - Ethanol
    - Methanol
    - Toluene 
    - Water
    '''
    print("Available compounds: Acetic acid, Acetone, Ammonia, Benzene, Carbon tetrachloride, Ethanol, Methanol, Toluene, Water")
    import numpy as np
    Aa = [15.8667, 4097.86, -27.4937]    
    Ac = [14.7171, 2975.95, -34.5228]    
    Am = [15.494, 2363.24, -22.6207]    
    Be = [14.1603, 2948.78, -44.5633]    
    Ct = [14.6247, 3394.46, -10.2163]    
    Et = [16.1952, 3423.53, -55.7152]    
    Me = [16.4948, 3593.39, -35.2249]    
    To = [14.2515, 3242.38, -47.1806]    
    Wa = [16.5362, 3985.44, -38.9974]    

    if name == "Acetic acid":
        A = Aa[0]
        B = Aa[1]
        C = Aa[2]
    elif name == "Acetone":
        A = Ac[0]
        B = Ac[1]
        C = Ac[2]
    elif name == "Ammonia":
        A = Am[0]
        B = Am[1]
        C = Am[2]
    elif name == "Benzene":
        A = Be[0]
        B = Be[1]
        C = Be[2]
    elif name == "Carbon tetrachloride":
        A = Ct[0]
        B = Ct[1]
        C = Ct[2]
    elif name == "Ethanol":
        A = Et[0]
        B = Et[1]
        C = Et[2]
    elif name == "Methanol":
        A = Me[0]
        B = Me[1]
        C = Me[2]
    elif name == "Toluene":
        A = To[0]
        B = To[1]
        C = To[2]
    elif name == "Water":
        A = Wa[0]
        B = Wa[1]
        C = Wa[2]
    else:
        print("Please check if the compound is listed above and double check the spelling thereof.")
    
    Psat = np.exp(A - B/(T + C))
    return Psat