import numpy as np



def opticsdata(wavs):

    nlambda=len(wavs)

    class crystal():

        name=np.asarray([ 'al2o3'],dtype=str)
        refrindex=np.full((nlambda,2),0.0,dtype=complex)
        thermopticalcoeff=np.full((8,2),0.0,dtype=float)
        transmission=np.full((nlambda),0.0,dtype=complex)
        thermexpansioncoeff=np.full((2),0.0,dtype=float)

        def __init__(self,name=[''],refrindex=np.full((2,nlambda),0.0,dtype=complex),\
                     thermopticalcoeff=np.full((2,8),0.0,dtype=float), \
                     transmission=np.full((nlambda),0.0,dtype=complex),\
                     thermexpansioncoeff=np.full((2),0.0,dtype=float)):
            self.name=name
            self.N=refrindex
            self.dNdT=thermopticalcoeff
            self.T=transmission
            self.tec=thermexpansioncoeff

    optics=[]
    for crys in [ 'al2o3', 'bab2o4', 'caco3', 'mgf2', 'sio2', 'tio2', 'yvo4', '', '', 'artif' ]:
        optics.append(crystal(name=crys))
    print(optics[0].N.shape)
    #------------------------------------------------------------------------------
    # Sellmeier equation: N[lambda] = SQRT(P1 + P2*lambda^2 / (lambda^2 - P3) $
    #                                         + P4*lambda^2 / (lambda^2 - P5) $
    #                                         + ...
    #                                         + P(N-1)*lambda^2 / (lambda^2 - PN) )
    # N[*,0] : ordinairy refractive index
    # N[*,1] : extraordinairy refractive index
    #------------------------------------------------------------------------------

    # Al2O3, Sapphire, 0.22 - 5.0 micron
    # Malitson, I. H. and Dodge, M. J., Refractive index and birefringence of synthetic sapphire,
    # J.Opt. Soc. Am. 62, 1405 (1972).
    optics[0].N[0:7,0] = [1., 1.43134936, 0.0726631**2., 0.65054713, 0.1193242**2., 5.3414021, 18.028251**2.]
    optics[0].N[0:7,1] = [1., 1.5039759, 0.0740288**2., 0.55069141, 0.1216529**2.,6.59273791, 20.072248**2.]
    optics[0].dNdT[0,:] = [293., 1.755, -45.2665E-06, 83.5457E-06, 8.27, 5.85,  7.21, -2.4]
    optics[0].dNdT[1,:] = [293., 1.748, -39.8961E-06, 81.9579E-06, 8.00, 5.42,  6.47, -2.2]
    optics[0].tec = 1E-6*[7.21, 6.47]
    #optics[0].dndtor = 13.6*1e-6 # @ 589 nm, RT
    #optics[0].dndtex = 14.7*1e-6 # @ 589 nm, RT

    ## BaB2O4, BBO, 0.22 - 1.06 micron
    #optics[1].N[0:5,0] = SQRT(2.7405 + 0.0184*lambda^2 / (lambda^2 - 0.0179) - 0.0155 * lambda^2)
    #optics[1].N[0:5,1] = SQRT(2.3730 + 0.0128*lambda^2 / (lambda^2 - 0.0156) - 0.0044 * lambda^2)

    # beta-BaB2O4, BBO, 0.22 - 1.06 micron - dn/dT values of beta-BBO
    # Eimerl, D., Davis, L., and Velsko, S., Optical, mechanical, and thermal properties of barium borate,
    # J. Appl. Phys. 62, 1968 (1987).
    # Ghosh, G, Temperature dispersion of refractive indices in βBaB2O4 and LiB3O5 crystals for nonlinear optical devices
    # J. Appl. Phys. 78, 6752 (1995)
    optics[1].N[0,0:5] = [1., 1.73651, 0.0120649, 0.0758505, 264.897]
    optics[1].N[1,0:5] = [1., 1.36847, 0.0100179,   1.29495, 187.560]
    optics[1].dNdT[0,:] = [293., 1.610, -19.3007E-6, -34.9683E-6, 19.00, 6.43,  4.00,  1.40]
    optics[1].dNdT[1,:] = [293., 1.520, -141.421E-6, 110.8630E-6, 17.00, 6.43, 36.00, -5.40]
    optics[1].tec = 1E-6*[4.00, 36.00]

    # CaCO3, Calcite, 0.2 - 2.2 micron
    # Gray, D. E., (Ed.), American Institute of Physics Handbook, 3rd ed. (McGraw-Hill, New York, 1972).
    optics[2].N[0,0:9] = [1., 0.8559, 0.0588**2., 0.83913, 0.141**2., 0.0009, 0.197**2., 0.6845, 7.005**2.]
    optics[2].N[1,0:7] = [1., 1.0856, 0.07897**2., 0.0988, 0.142**2., 0.317, 1.468**2.]
    optics[2].dNdT[0,:] = [293., 1.613, -121.689E-06, 122.494E-06, 10.80, 10.00, 25.00, -7.60]
    optics[2].dNdT[1,:] = [293., 1.472,  12.7011E-06, 20.4803E-06,  9.05,  6.83, -3.70, -1.20]
    optics[2].tec = 1E-6*[25.00, -3.70]

    # MgF2, 0.4 - 3.1 micron
    # Dodge, M. J., Refractive properties of magnesium fluoride,
    # Appl. Opt. 23, 1980 (1984).
    optics[3].N[0,0:7] = [1., 0.48755108, 0.04338408**2., 0.39875031, 0.09461442**2., 2.3120353, 23.793604**2.]
    optics[3].N[1,0:7] = [1., 0.41344023, 0.03684262**2., 0.50497499, 0.09076162**2., 2.4904862, 23.771995**2.]
    optics[3].dNdT[0,:] = [293., 1.290, -37.2043E-06, 39.3186E-06, 13.10, 8.00,  9.3, -4.70]
    optics[3].dNdT[1,:] = [293., 1.290, -56.7859E-06, 57.3986E-06, 15.50, 8.00, 14.2, -6.90]
    optics[3].tec = 1E-6*[9.3, 14.2]
    # lattice/ionic contribution: (o: L=15.9824E-6, l_ip = 40.0) (e: L=8.9996E-6, l_ip = 40.0)
    #optics[3].dndtor = 1.12*1e-6 # @ 633 nm, RT
    #optics[3].dndtex = 0.58*1e-6 # @ 633 nm, RT

    # SiO2, Quartz, 0.18 - 0.71 micron
    # Radhakrishnan, T., Further studies on the temperature variation of # the refractive index of crystals,
    # Proc. Indian Acad. Sci., A33, 22 (1951).
    #optics[4].N[0:10,0] = [1., 0.663044, 0.060^2, 0.517852, 0.106^2, 0.175912, 0.119^2, 0.565380, 8.844^2, 1.675299, 20.742^2]
    #optics[4].N[0:10,1] = [1., 0.665721, 0.060^2, 0.503511, 0.106^2, 0.214792, 0.119^2, 0.539173, 8.792^2, 1.807613, 197.70^2]

    # SiO2, Quartz, 0.18 - 0.71 micron
    # Gosh 1999 (Halle)
    optics[4].N[0,0:5] = [1.28604141, 1.07044083, 1.00585997*0.01, 1.10202242, 100.]
    optics[4].N[1,0:5] = [1.28851804, 1.09509924, 1.02101864*0.01, 1.15662475, 100.]
    optics[4].dNdT[0,:] = [293., 1.515, -61.184E-06, 43.9990E-06, 10.30, 8.90,  6.88, -3.02]
    optics[4].dNdT[1,:] = [293., 1.520, -70.1182E-06, 49.2875E-06, 10.30, 8.90, 12.38, -3.32]
    optics[4].tec = 1E-6*[6.88, 12.38]
    # fit expansion coef. 15.70 (meas. 6.88) and 17.80 (meas. 12.38)

    # TiO2, Rutile, 0.43 - 1.5 micron
    # DeVore,J. R., Refractive index of rutile and sphalerite,
    # J. Opt. Soc. Am. 41, 416 (1951).
    optics[5].N[0,0:5] = [1., 4.913, 0.0, 0.2441, 0.0803]
    optics[5].N[1,0:5] = [1., 6.097, 0.0, 0.3322, 0.0843]
    optics[5].dNdT[0,:] = [293., 2.432, -132.253E-06, 64.5269E-06, 4.10, 3.50,  8.98, -0.46]
    optics[5].dNdT[1,:] = [293., 2.683, -127.565E-06, 45.2141E-06, 4.10, 3.50,  6.87, -0.26]
    optics[5].tec = 1E-6*[8.98, 6.87]

    # YVO4, 0.5 - 1.06 micron, HoOM Ref 112/113 - dn/dT values of TiO2!
    # Maunder, E. A. and DeShazer, L. G., Use of yttrium orthovanadate for # polarizers,
    # J. Opt. Soc. Am. 61, 684A (1971).
    # Lomheim, T. S. and DeShazer, L. G., Optical absorption intensities
    # of trivalent neodymium in the uniaxial crystal yttrium orthovanadate,
    # J. Appl. Phys. 49, 5517 (1978).
    optics[6].N[0,0:3] = [1., 2.7665, 0.026884]
    optics[6].N[1,0:3] = [1., 3.5930, 0.032103]
    optics[6].dNdT[0,:] = [293., 2.432, -132.253E-06, 64.5269E-06, 4.10, 3.50,  8.98, -0.46]
    optics[6].dNdT[0,:] = [293., 2.683, -127.565E-06, 45.2141E-06, 4.10, 3.50,  6.87, -0.26]
    optics[6].tec = 1E-6*[8.98, 6.87]

    # Artificial material with spectrally constant refractive index of 1.5 and birefringence of 0.01
    optics[9].N[0,0] = [2.25]
    optics[9].N[1,0] = [2.2801]
    optics[9].dNdT[0,:] = [293., 1.5, 0.0E-06, 0.0E-06, 4.10, 3.50,  8.98, -0.46]
    optics[9].dNdT[0,:] = [293., 1.51, 0.0E-06, 0.0E-06, 4.10, 3.50,  6.87, -0.26]
    optics[9].tec = 1E-6*[0.0, 0.0]
    
    print(vars(optics[0]).items())

    return optics

def demod(MORinput,lambdainp, spectra, fitout, conf=[1,0,0], \
          lrange=[390,820], aolp=0, winstep=0, nodeltafit=0):
    #print('lrange',lrange)
    #print('lambdainp',lambdainp[0],lambdainp[-1])
    #Select data in wavelength range
    nlambda=len(lambdainp)
    lrange = [lrange[0]-10,lrange[1]+20] #widen the user supplied wavelength range
    if lrange[0] < lambdainp[0]:
        print('   Warning: adjusting lambda_min to avoid edge effect to : %i nm, (was %i nm.)' % (lambdainp[0]+10,lambdainp[0]))
        lrange[0]=lambdainp[0]+10.
    if lrange[1] > lambdainp[-1]:
        print('   Warning: adjusting lambda_max to avoid edge effect to: %i nm, (was %i nm.)' % (lambdainp[-1]-20.,lambdainp[-1]))
        lrange[1]=lambdainp[-1]-20.

    lmin=np.min(np.argwhere(lambdainp >= (lrange[0]-10.) ))
    lmax=np.min(np.argwhere(lambdainp >= (lrange[1]+20.) ))
    lminout = np.min(np.argwhere(lambdainp >= (lrange[0])))
    lmaxout = np.min(np.argwhere(lambdainp >= (lrange[1])))
    Lambda = lambdainp[lmin:1+lmax]
    nlambda = len(Lambda)
    """
    print('lmin',lmin)
    print('lmax',lmax)
    print('lminout',lminout)
    print('lmaxout',lmaxout)
    print(Lambda)
    print(nlambda)
    """

    szspectra=spectra.shape
    if szspectra[0] > 1:
        nspect=szspectra[0]
    else:
        nspect=1

    optics =opticsdata(Lambda)
    
    return None
