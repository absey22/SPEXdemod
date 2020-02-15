import numpy as np
from classobjects import crystal,spectralwindow


def opticsdata(lambdainp):

    nlambda=np.size(lambdainp)
  
    #initialize instances of each crystal, to fill in with parameters below
    optics=[]
    for crys in [ 'al2o3', 'bab2o4', 'caco3', 'mgf2', 'sio2', 'tio2', 'yvo4', '', '', 'artif' ]:
        optics.append(crystal(name=crys,Nwavs=nlambda))
    
    #------------------------------------------------------------------------------
    # Sellmeier equation: N[lambda] = SQRT(P1 + P2*lambda^2 / (lambda^2 - P3) $
    #                                         + P4*lambda^2 / (lambda^2 - P5) $
    #                                         + ...
    #                                         + P(N-1)*lambda^2 / (lambda^2 - PN) )
    # N[0,:,0] : ordinairy refractive index
    # N[1,:,0] : extraordinairy refractive index
    #------------------------------------------------------------------------------

    # Al2O3, Sapphire, 0.22 - 5.0 micron
    # Malitson, I. H. and Dodge, M. J., Refractive index and birefringence of synthetic sapphire,
    # J.Opt. Soc. Am. 62, 1405 (1972).
    
    optics[0].N[0,0:7,0] = [1., 1.43134936, 0.0726631**2., 0.65054713, 0.1193242**2., 5.3414021, 18.028251**2.]
    optics[0].N[1,0:7,0] = [1., 1.5039759, 0.0740288**2., 0.55069141, 0.1216529**2.,6.59273791, 20.072248**2.]
    optics[0].dNdT[0,:] = [293., 1.755, -45.2665E-06, 83.5457E-06, 8.27, 5.85,  7.21, -2.4]
    optics[0].dNdT[1,:] = [293., 1.748, -39.8961E-06, 81.9579E-06, 8.00, 5.42,  6.47, -2.2]
    optics[0].tec = 1E-6*np.asarray([7.21, 6.47])
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
    optics[1].N[0,0:5,0] = [1., 1.73651, 0.0120649, 0.0758505, 264.897]
    optics[1].N[1,0:5,0] = [1., 1.36847, 0.0100179,   1.29495, 187.560]
    optics[1].dNdT[0,:] = [293., 1.610, -19.3007E-6, -34.9683E-6, 19.00, 6.43,  4.00,  1.40]
    optics[1].dNdT[1,:] = [293., 1.520, -141.421E-6, 110.8630E-6, 17.00, 6.43, 36.00, -5.40]
    optics[1].tec = 1E-6*np.asarray([4.00, 36.00])
    
    # CaCO3, Calcite, 0.2 - 2.2 micron
    # Gray, D. E., (Ed.), American Institute of Physics Handbook, 3rd ed. (McGraw-Hill, New York, 1972).
    optics[2].N[0,0:9,0] = [1., 0.8559, 0.0588**2., 0.83913, 0.141**2., 0.0009, 0.197**2., 0.6845, 7.005**2.]
    optics[2].N[1,0:7,0] = [1., 1.0856, 0.07897**2., 0.0988, 0.142**2., 0.317, 1.468**2.]
    optics[2].dNdT[0,:] = [293., 1.613, -121.689E-06, 122.494E-06, 10.80, 10.00, 25.00, -7.60]
    optics[2].dNdT[1,:] = [293., 1.472,  12.7011E-06, 20.4803E-06,  9.05,  6.83, -3.70, -1.20]
    optics[2].tec = 1E-6*np.asarray([25.00, -3.70])
    
    # MgF2, 0.4 - 3.1 micron
    # Dodge, M. J., Refractive properties of magnesium fluoride,
    # Appl. Opt. 23, 1980 (1984).
    optics[3].N[0,0:7,0] = [1., 0.48755108, 0.04338408**2., 0.39875031, 0.09461442**2., 2.3120353, 23.793604**2.]
    optics[3].N[1,0:7,0] = [1., 0.41344023, 0.03684262**2., 0.50497499, 0.09076162**2., 2.4904862, 23.771995**2.]
    optics[3].dNdT[0,:] = [293., 1.290, -37.2043E-06, 39.3186E-06, 13.10, 8.00,  9.3, -4.70]
    optics[3].dNdT[1,:] = [293., 1.290, -56.7859E-06, 57.3986E-06, 15.50, 8.00, 14.2, -6.90]
    optics[3].tec = 1E-6*np.asarray([9.3, 14.2])
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
    optics[4].N[0,0:5,0] = [1.28604141, 1.07044083, 1.00585997*0.01, 1.10202242, 100.]
    optics[4].N[1,0:5,0] = [1.28851804, 1.09509924, 1.02101864*0.01, 1.15662475, 100.]
    optics[4].dNdT[0,:] = [293., 1.515, -61.184E-06, 43.9990E-06, 10.30, 8.90,  6.88, -3.02]
    optics[4].dNdT[1,:] = [293., 1.520, -70.1182E-06, 49.2875E-06, 10.30, 8.90, 12.38, -3.32]
    optics[4].tec = 1E-6*np.asarray([6.88, 12.38])
    # fit expansion coef. 15.70 (meas. 6.88) and 17.80 (meas. 12.38)
    
    # TiO2, Rutile, 0.43 - 1.5 micron
    # DeVore,J. R., Refractive index of rutile and sphalerite,
    # J. Opt. Soc. Am. 41, 416 (1951).
    optics[5].N[0,0:5,0] = [1., 4.913, 0.0, 0.2441, 0.0803]
    optics[5].N[1,0:5,0] = [1., 6.097, 0.0, 0.3322, 0.0843]
    optics[5].dNdT[0,:] = [293., 2.432, -132.253E-06, 64.5269E-06, 4.10, 3.50,  8.98, -0.46]
    optics[5].dNdT[1,:] = [293., 2.683, -127.565E-06, 45.2141E-06, 4.10, 3.50,  6.87, -0.26]
    optics[5].tec = 1E-6*np.asarray([8.98, 6.87])
    
    # YVO4, 0.5 - 1.06 micron, HoOM Ref 112/113 - dn/dT values of TiO2!
    # Maunder, E. A. and DeShazer, L. G., Use of yttrium orthovanadate for # polarizers,
    # J. Opt. Soc. Am. 61, 684A (1971).
    # Lomheim, T. S. and DeShazer, L. G., Optical absorption intensities
    # of trivalent neodymium in the uniaxial crystal yttrium orthovanadate,
    # J. Appl. Phys. 49, 5517 (1978).
    optics[6].N[0,0:3,0] = [1., 2.7665, 0.026884]
    optics[6].N[1,0:3,0] = [1., 3.5930, 0.032103]
    optics[6].dNdT[0,:] = [293., 2.432, -132.253E-06, 64.5269E-06, 4.10, 3.50,  8.98, -0.46]
    optics[6].dNdT[0,:] = [293., 2.683, -127.565E-06, 45.2141E-06, 4.10, 3.50,  6.87, -0.26]
    optics[6].tec = 1E-6*np.asarray([8.98, 6.87])

    # Artificial material with spectrally constant refractive index of 1.5 and birefringence of 0.01
    optics[9].N[0,0,0] = 2.25
    optics[9].N[1,0,0] = 2.2801
    optics[9].dNdT[0,:] = [293., 1.5, 0.0E-06, 0.0E-06, 4.10, 3.50,  8.98, -0.46]
    optics[9].dNdT[0,:] = [293., 1.51, 0.0E-06, 0.0E-06, 4.10, 3.50,  6.87, -0.26]
    optics[9].tec = 1E-6*np.asarray([0.0, 0.0])
    
    return optics


def demod_refrindex(lambdainp,params1,params2,temperature):
    nlambda=np.size(lambdainp)
    lamb=lambdainp/1000.

    Nabs=np.zeros((2,nlambda,2))
    dNabsdT=np.zeros((2,nlambda))
    
    # N[0,:,0] : ordinairy refractive index
    # N[1,:,0] : extraordinairy refractive index
    #Calculate refractive index as a function of wavelength via Sellmeier equation
    Nabs[0,:,0] = params1[0,0,0]
    Nabs[1,:,0] = params1[1,0,0]
    
    #print(int((np.size(np.argwhere(params1[0,:,0] != 0))-1)/2))
    
    for jj in range(int((np.size(np.argwhere(params1[0,:,0] != 0))-1)/2)):
        
        Nabs[0,:,0] += params1[0,2*jj+1,0]*lamb**2. / (lamb**2. - params1[0,2*jj+2,0])
        Nabs[1,:,0] += params1[1,2*jj+1,0]*lamb**2. / (lamb**2. - params1[1,2*jj+2,0])

    
    Nabs=np.sqrt(Nabs)
    

    #Calculate temperature correction w.r.t. 293 K
    dT = temperature - params2[0,0]
    
    if params2[0,7] == 0:
        #Using Schott's Sellmeier-type equation for glasses
        dNabsdT[0,:] = (Nabs[0,:,0]**2. - 1.) / (2.*Nabs[0,:,0]) * (params2[0,1]*dT + params2[0,2]*dT**2. + params2[0,3]*dT**3. \
             + (params2[0,4]*dT + params2[0,5]*dT**2.) / (lamb**2. - (params2[0,6])**2.) )

        dNabsdT[1,:] = (Nabs[1,:,0]**2. - 1.) / (2.*Nabs[1,:,0]) * (params2[1,1]*dT + params2[1,2]*dT**2. + params2[1,3]*dT**3. \
             + (params2[1,4]*dT + params2[1,5]*dT**2.) / (lamb**2. - (params2[1,6])**2.) )
    else:
        #Using dNdT equation of Gosh
        Econv = 1.E9 * 6.626E-34 * 2.9979E8 / 1.6E-19
        lambdaig = Econv / params2[:,4]
        Ro = lambdainp**2. / (lambdainp**2. - lambdaig[0]**2.)
        Re = lambdainp**2. / (lambdainp**2. - lambdaig[1]**2.)
        
        dNabsdT[0,:] = (params2[0,2] * Ro + params2[0,3] * Ro**2.) / (2. * Nabs[0,:,0])
        dNabsdT[1,:] = (params2[1,2] * Re + params2[1,3] * Re**2.) / (2. * Nabs[1,:,0])

    N=Nabs+np.transpose(dNabsdT)[np.newaxis,:]*dT
    
    return N




def demod_get_delta(MORinput,lambdainp,optics,N=np.zeros((1))):
    """
    MORattrs=np.asarray(list(MORinput.__dict__.keys()))
    print(MORattrs)
    comp=np.where(np.any(MORattrs=='comp'),MORinput.comp,'al2o3')
    
    #doesn't work for values of ncomp other than 2...
    #comp=np.asarray(['aaaa','bbbb','cccc'])
    ncomp=np.size(comp)
    thick=np.where(np.any(MORattrs=='thick'),MORinput.thick,np.ones((ncomp)))
    #par=np.where(attrs=='par',np.ones(( 2,np.where(ncomp < 2, 2, ncomp))),MORinput.par)
    par=np.where(np.any(MORattrs=='par')&MORinput.par!=None,np.zeros((2,ncomp)),MORinput.par) #DEFAULTS
    rot=np.where(np.any(MORattrs=='rot'),np.zeros(ncomp),MORinput.rot)    #SWITCHED
    temp=np.where(np.any(MORattrs=='temp'),np.zeros((ncomp))+293.,MORinput.temp)
    """
    #collect config of MOR from input class structure
    comp=MORinput.comp
    ncomp=np.size(comp)
    thick=MORinput.thick
    #par = (WHERE(tags eq 'PAR') eq -1) ? FLTARR(ncomp > 2,2)+0 : input.par
    par=MORinput.par
    if par.shape[0]<ncomp:
        par=np.zeros((2,ncomp))
    rot=MORinput.rot
    temp=MORinput.temp

    
    if np.size(N)<2: 
        Nset=False
    else:
        Nset=True
    
    nlambda=np.size(lambdainp)
    delta=np.zeros(nlambda)

    #loop over optical elements
    crystalnames=np.asarray([optics[i].name for i in range(len(optics))])
    for jj in range(ncomp):
        
        #calculate refractive index
        i_crystal=np.argwhere(crystalnames==comp[jj])[0][0]
        print(comp[jj],'   ',optics[i_crystal].name)
        params1= np.where(Nset,0,np.copy(optics[i_crystal].N))
        params2= np.where(Nset,0,np.copy(optics[i_crystal].dNdT))
        N1=np.where(Nset,N, demod_refrindex(lambdainp,params1,params2,temp[jj]))

        #print('0',N1[0,0])
        #print('1',N1[1,0])

        #Calculate birefringence
        biref1 = N1[1,:,0] - N1[0,:,0]
        
        #Calculate thermal expansion
        a = optics[i_crystal].tec[0]
        
        d1 = 1E6*thick[jj]*(1.+a*(temp[jj]-293.))
        #print(d1)
        #Calculate retardance including field-of-view behavior and thermal expansion
        phi = par[0,jj]

        
        print(par[1,jj] , rot[jj])
        theta = np.where(biref1[0] > 0, par[1,jj] + rot[jj], par[1,jj] + np.pi/2. + rot[jj])
        
        c1 = np.sin(phi)*np.sin(theta)
        c2 = np.sin(phi)*np.cos(theta)
        p0 = N1[0,:,0] * np.sqrt(1. - (c1/N1[0,:,0])**2. - (c2/N1[0,:,0])**2.)
        p1 = N1[1,:,0] * np.sqrt(1. - (c1/N1[0,:,0])**2. - (c2/N1[1,:,0])**2.) #typo in N1?
        delta1 = np.cos(np.pi*rot[jj]/90.) * (biref1 * d1 / lambdainp) * ( p1 - p0 ) / biref1
  
        #Add retardance of each component        
        delta += delta1
        
  
    return (delta*lambdainp)*(1.-np.zeros(nlambda)/(nlambda-1.))

                
def demod_spc_win(lambdainp,delta,GvHmethod=True):
    nlambda=np.size(lambdainp)
                
    window=spectralwindow(nlambda)

    if GvHmethod:
        #GvH wing it method
        k=delta/lambdainp
        for i in range(nlambda):
            window.minimum[i]= np.argwhere( (np.abs(k-(k[i]+0.5)) == np.min(np.abs(k-(k[i]+0.5))) ) )
            window.maximum[i]=np.argwhere( (np.abs(k-(k[i]-0.5)) == np.min(np.abs(k-(k[i]-0.5))) ) )
        window.length=window.maximum-window.minimum+1.
    else:
        
        #Calculate spectral resolution of wavelength array
        res=(100./(lambdainp[np.clip((np.arange(nlambda)+1),a_min=None,a_max=(nlambda-1) )]-lambdainp))/100.
        res[nlambda-1]=res[nlambda-2]

        #Calculate spectral window length
        length=res*lambdainp**2./(np.abs(retardance)*(1.+(lambdainp/(2.*np.abs(retardance)))**2.))
        #Calculate position of left and right boundary of each spectral window
        window.minimum = np.clip(np.arange(nlambda)-length/2.,a_min=0.0,a_max=None)
        window.maximum = np.clip(np.arange(nlambda)+length/2.,a_min=None,a_max=(nlambda-1))
        window.length = window.maximum - window.minimum + 1
        
    return window


def demod(MORinput,lambdainp, spectra, fitout, conf=[1,0,0], \
          lrange=[390,820], aolp=0, winstep=0, nodeltafit=False, nlamb=False):
    #print('lrange',lrange)
    #print('lambdainp',lambdainp[0],lambdainp[-1])
    #Select data in wavelength range
    nlambda=np.size(lambdainp)
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
    nlambda = np.size(Lambda)
    spectra = spectra[:,:,lmin:1+lmax]
    

    szspectra=spectra.shape
    if szspectra[0] > 1:
        nspect=szspectra[0]
    else:
        nspect=1

    optics =opticsdata(Lambda)
    
    #for i in range(10):
    #    print(optics[i].name,optics[i].dNdT)

    #Calculate spectral retardance based on input composition of multiple-order retarder
    deltalit = demod_get_delta(MORinput, Lambda, optics)

    print("Mean retardance = %1.2f" % np.mean(deltalit))
    
    #Calculate spectral window lengths and positions
    spcwin = demod_spc_win(Lambda, deltalit)

    if nodeltafit:
        delta=deltalit
    else:
        #Calculate normalized spectrum of calibration spectrum
        spcnorm=demod_spc_norm(spcwin,Lambda,spectra)
        #Calculate number of point to be used in the interpolation in the full spectrum fit
        if not nlamb:
            nlamb=np.clip( np.abs(np.ceil( 0.6*np.mean(deltalit)*(1./Lambda[0]-1./Lambda[np.size(Lambda)-1]))), a_min=4,a_max=40)
        print(nlamb)
        #Perform full spectrum fit in order to obtain best fit of the spectral retardance
        #delta=demod_delta()
        #Recalculate spectral window lengths and positions
        #spcwin=demod_spc_win()
        
    return None




