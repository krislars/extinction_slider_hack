def getMags_ext(wavelength,flux):
    
    import numpy as np
    import scipy.interpolate as interp
    import scipy.integrate as integrate
    
    """ This function converts an extinction curve to a delta magnitude.
        Because of this, the zero points are not needed.
        Inputs: band character ID, wavelength in microns, and flux is A(lambda)
    """
    
    # Optical filter curves from http://www.aip.de/en/research/facilities/stella/instruments/data
    #     Converting from nm to microns
    # SDSS centers/widths/F0 from http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
    # NOTE: ugriz filters are on the AB magnitude system, KsHJUBVRI are on the Vega system

    mags = {}

    for band in ['Ks','H','J','R','V','B']:
        
        if band=='Ks':
            bandwav,bandpass=np.loadtxt('filters/Ks_2MASS.txt',unpack=True)
            center=2.159        # micron
            F0=4.283E-10        # W m^-2 micron^-1
            dlambda=0.262

        elif band=='H':
            bandwav,bandpass=np.loadtxt('filters/H_2MASS.txt',unpack=True)
            center=1.662
            F0=1.133E-9
            dlambda=0.251

        elif band=='J':
            bandwav,bandpass=np.loadtxt('filters/J_2MASS.txt',unpack=True)
            center=1.235
            F0=3.129E-9
            dlambda=0.162

        # elif band=='U':
        #    bandwav_A,bandpass_100=np.loadtxt('filters/Bessel_U.txt',unpack=True)
        #    bandwav=bandwav_A/1.0E3
        #    bandpass=bandpass_100/100.0
        #    center=0.365
        #    F0=4.18E-8
        #    dlambda=0.066

        elif band=='B':
            bandwav_A,bandpass_100=np.loadtxt('filters/Bessel_B.txt',unpack=True)
            bandwav=bandwav_A/1.0E3
            bandpass=bandpass_100/100.0
            center=0.445
            F0=6.32E-8
            dlambda=0.094

        elif band=='V':
            bandwav_A,bandpass_100=np.loadtxt('filters/Bessel_V.txt',unpack=True)
            bandwav=bandwav_A/1.0E3
            bandpass=bandpass_100/100.0
            center=0.551
            F0=3.63E-8
            dlambda=0.088

        elif band=='R':
            bandwav_A,bandpass_100=np.loadtxt('filters/Bessel_R.txt',unpack=True)
            bandwav=bandwav_A/1.0E3
            bandpass=bandpass_100/100.0
            center=0.658
            F0=2.18E-8
            dlambda=0.138

        # elif band=='I':
        #    bandwav_A,bandpass_100=np.loadtxt('filters/Bessel_I.txt',unpack=True)
        #    bandwav=bandwav_A/1.0E3
        #    bandpass=bandpass_100/100.0
        #    center=0.806
        #    F0=1.13E-8
        #    dlambda=0.149
        

        filterband=np.zeros(wavelength.size)
        
        bandinterp=interp.interp1d(bandwav,bandpass)
            #1D function between bandwav(x) and bandpass(y), y=f(x)
        
        if bandwav[0]<bandwav[-1]: # increasing wavelength
            inband=np.logical_and(wavelength>bandwav[0],wavelength<bandwav[-1] )
        else: # decreasing wavelength
            inband=np.logical_and(wavelength<bandwav[0],wavelength>bandwav[-1] )
            #Array of 'True's where wavelength is in range of bandwav, 'False's outside of that range. Basically just an index of bandwav in wavelength array
        
        filterband[inband]=bandinterp(wavelength[inband])
            #Ignore all values of bandpass outside of filter range. Now filterband is the same shape as wavelength and flux
    
        dwav=np.zeros(wavelength.size)
        
        dwav[0:-1]=np.abs(wavelength[1:]-wavelength[0:-1])
    
        mags[band]=-2.5*np.log10((integrate.trapz(10.0**(-flux*0.4)*filterband,x=wavelength))/integrate.trapz(filterband,x=wavelength))


    return mags
