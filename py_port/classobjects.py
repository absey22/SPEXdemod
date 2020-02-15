import numpy as np


#setting up parameters of the multiple order retarder, with various crystals (see demod_abs.py)
class MOR:

    def __init__(self,composition=None,thickness=None,par=None,rot=None,temp=None ):
        #if no crystal names are included when instantiating the class
        compdefault=['al2o3','mgf2']
        ncomp=len(compdefault)
        
        if composition is None:
            self.comp=compdefault
        else:
            self.comp=composition
        if thickness is None:
            self.thick=[1.2,2.88]
        else:
            self.thick=thickness
        #not specified with default in GvH code originally.
        if par is None:
            self.par=np.zeros((2,ncomp))
        else:
            self.par=par
        if rot is None:
            self.rot=np.zeros(ncomp)
        else:
            self.par=par
        if temp is None:
            self.temp=np.zeros((ncomp))+293.
        else:
            self.temp=temp


#object for constructing a crystal as part of the MOR in the SPEX optics (see demod_abs.py)
class crystal:
        #name=np.asarray([ 'al2o3'],dtype=str)
        #refrindex=np.full((nlambda,2),0.0,dtype=complex)
        #thermopticalcoeff=np.full((8,2),0.0,dtype=float)
        #transmission=np.full((nlambda),0.0,dtype=complex)
        #thermexpansioncoeff=np.full((2),0.0,dtype=float)

        #np.full((2,nlambda),0.0,dtype=complex)
        def __init__(self,name, Nwavs,refrindex=None,thermopticalcoeff=None, \
                     transmission=None,thermexpansioncoeff=None ):
            #None's used to counteract python 'mutable' types behavior
            self.name=name  #crystal name
            if refrindex is None:
                self.N=np.zeros((2,Nwavs,2))      #Refractive index
            else:
                self.N=refrindex
            if thermopticalcoeff is None:
                self.dNdT=np.zeros((2,8))   #Thermal optical coefficient
            else:
                seld.dNdT=thermopticalcoeff
            if transmission is None:
                self.T=np.zeros((2,Nwavs))      #Transmission
            else:
                self.T=transmission
            if thermexpansioncoeff is None:
                self.tec=np.zeros((2))    #Thermal expansion coefficient
            else:
                self.tec=thermexpansioncoeff


            

#object to hold spectral window sizes for fitting (see demod_abs.py)
class spectralwindow:
        def __init__(self,Nwavs=None):
            if Nwavs is None:
                print("Length of wavelength array not specified.")
                self.minimum=np.zeros(1)
                self.maximum=np.zeros(1)
                self.length=np.zeros(1)
            else:
                self.minimum=np.zeros(Nwavs)
                self.maximum=np.zeros(Nwavs)
                self.length=np.zeros(Nwavs)
