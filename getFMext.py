def getFMext(wave,R,source):
    
    """ 
        Input: wavelength in microns, nominal R = 3.1
            source='f99' for Fitzpatrick 1999 or ='fmunred' for fmunred.pro
            source='f99 tables' reproduces common misunderstanding from that paper
        Output: Al/EBV, so user must divide by R to get Al/AV !!!
    """
    
    from scipy.interpolate import interp1d
    
    import numpy as np

    x_anchors = 1.0E4 / np.array([np.inf, 26500., 12200., 6000., 5470., 4670., 4110.]) # microns
    
    if source=="f99 tables": # Don't use, for demonstration only
        a26500= 0.265
        a12200= 0.829
        a6000 = -0.426 +1.0044*R
        a5470 = -0.050 +1.0016*R
        a4670 =  0.701 +1.0016*R
        a4110 =  1.208 +1.0032*R -0.00033*R**2.0 # typo in the paper -KAL
    
    elif source=="f99": 
        a26500= 0.265*R/3.1
        a12200= 0.829*R/3.1
        a6000 = -0.426 +1.0044*R
        a5470 = -0.050 +1.0016*R
        a4670 =  0.701 +1.0016*R
        a4110 =  1.208 +1.0032*R -0.00033*R**2.0 # typo in the paper -KAL

    elif source=="fmunred":
        a26500= 0.26469*R/3.1
        a12200= 0.82925*R/3.1
        a6000 = -4.22809e-01 +1.00270*R+2.13572e-04*R**2.0
        a5470 = -5.13540e-02 +1.00216*R-7.35778e-05*R**2.0
        a4670 =  7.00127e-01 +1.00184*R-3.32598e-05*R**2.0
        a4110 =  1.19456 +1.01707*R -5.46959e-03*R**2.0+ 7.97809e-04*R**3.0 -4.45636e-05*R**4.0

    a_anchors = np.array([0.0, a26500, a12200, a6000, a5470, a4670, a4110])

    f=interp1d(x_anchors, a_anchors, kind='cubic')    
    
    
    return f(1.0/wave)