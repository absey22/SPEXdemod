import matplotlib.pyplot as plt
import numpy as np

import sys
import glob

from scipy.io import readsav


if len(sys.argv) == 1:
    parentdir='/home/seymour/Documents/SPEX/demodulation_pipeline/SPEXredux/ScriptRun_09072013 105035'
else:
    parentdir=sys.argv[1]

    
if str(parentdir)[-1]!='/':
    parentdir+='/'



#take a list of filenames and sorts them according to the info stored in the file names
#from AvaSpec
def sortfilenames(filenames):
    exps=np.zeros((2*len(texp)))
    channels=exps.copy()
    for i,fname in enumerate(filenames): #extract spectral channel and exposure 'index' info from file list
        channels[i]=fname[fname.find('U2_')-1] #key for sorting by channel
        exps[i]=fname.split('_')[-8] # key for sorting by exposure index number

    expkey=np.argsort(exps) #return indices ordered for sorting by the exposures

    filenames=filenames[expkey] #sort the file list in ascending exposure    
    channels=channels[expkey] #and reorder channelkey accordindly

    channelkey=np.zeros_like(expkey)
    for i in range(0,len(channels),2): #sort channels in each successive exposure
        channelkey[i:i+2]=i+np.argsort(channels[i:i+2])

    filenames=filenames[channelkey] #sort the (sorted) file list in ascending channels
    #print(np.asarray(filenames).reshape(len(texp),2))
    return filenames

#takes a list of filenames and creates a numpy array of rows with the data in each file
#RESHAPE to exposure rows, and spectral channel columns
def importfilenames(filenames):
    spectraldata=[]
    for fname in filenames:
        spectraldata.append(np.genfromtxt(fname,delimiter=',')[:-1])#cut of nan at the end
    #RESHAPE step
    spectraldata=np.asarray(spectraldata).reshape(len(texp),2,spectraldata[0].shape[-1])
    return spectraldata

def poly2d(x,y,coeffs):
    z = np.zeros((np.size(x),np.size(y)))
    x,y=[x],[y]
    for xx in range(np.size(x)):
       for yy in range(np.size(y)):
          for i in range(coeffs.shape[1]):
             for j in range(coeffs.shape[0]):
                z[yy,xx] += coeffs[i,j] * ((x[xx])**i) * ((y[yy])**j)
    return z


    
#================ IMPORT META DATA

metafile=glob.glob(parentdir+"Meta*.*")

with open(metafile[0], encoding="utf8", errors='ignore') as meta:
    meta=meta.read().splitlines()
  
pantilt=np.nan*np.empty((len(meta),2))
temperature=np.nan*np.empty((len(meta),2))
texp=np.nan*np.empty((len(meta)))

for i,line in enumerate(meta):
    line=line.split('!')
    if line[1].find('Calculated integration time')!=-1:
        texp[i]=float(line[0])
        
    if line[1].find('At horizontal/azimuthal pos (Pan)')!=-1:
        pantilt[i,0]=float(line[0])
    if line[1].find('At vertical/elevational pos (Tilt)')!=-1:
        pantilt[i,1]=float(line[0])

    if line[1].find('Detector 0 temp ')!=-1:
       temperature[i,0]=float(line[0])
    if line[1].find('Detector 1 temp ')!=-1:
        temperature[i,1]=float(line[0])
       
 #reformat 
texp=texp[~np.isnan(texp)]
pantilt=pantilt[~np.isnan(pantilt)].reshape(int(len(pantilt[~np.isnan(pantilt)])/2),2)
pantilt=pantilt[:len(texp),:] #crop down to the amount as per number of texp's
temperature=temperature[~np.isnan(temperature)].reshape(int(len(temperature[~np.isnan(temperature)])/2),2)
temperature=temperature[:len(texp),:] #crop down to the amount as per number of texp's

"""
print('pantilt')
print(pantilt)
print('temp')
print(temperature)
print('texp')
print(texp)
"""



#================ IMPORT SPECTRA

specfiles=np.asarray(glob.glob(parentdir+'Spectrometer_110516*U2_*pix.txt'))

#sort in ascending exposure, and then ascending channel
#(so that the order matches with the exp times and temps in the Meta data)
specfiles=sortfilenames(specfiles)

#import the spectral data into an multi-dim array
spec=importfilenames(specfiles)


#================ IMPORT DATE + TIMES

times=[]
for specfile in specfiles:
    t=specfile.split(' ')[2].split('_')[0] #extract the 6 digit timestamp
    t=float(t[0:2])+float(t[2:4])/60.+float(t[4:])/3600. #convert to a decimal type hour
    times.append(round(t,4))               
times=times[:len(texp)]

date=np.asarray([specfiles[0].split(' ')[0].split('_')[-1][:2], \
                 specfiles[0].split(' ')[0].split('_')[-1][2:4], \
                 specfiles[0].split(' ')[0].split('_')[-1][4:] ], dtype=float)#extract 8 digit day stamp


#================ IMPORT BLACK PIXELS

blackfiles=np.asarray(glob.glob(parentdir+'Spectrometer_110516*U2_*pix.dark13.txt'))

#sort in ascending exposure, and then ascending channel
blackfiles=sortfilenames(blackfiles)

#import the spectral data into an multi-dim array
black=importfilenames(blackfiles)



#================ DARK SUBTRACTION


darkmap=readsav('/home/seymour/Documents/SPEX/demodulation_pipeline/darkmap.sav')

darkmodblack=darkmap.darkmodblack
darkmodspec=darkmap.darkmodspec

#darkmodblack has python shape (#exps,13,5,5)
#darkmodblack has   IDL  shape (5,5,13,#exps)

darkblack=np.zeros((len(texp),2,darkmodblack.shape[1]))
darkspec=np.zeros((len(texp),2,darkmodspec.shape[1]))
#for each spectral channel (ch) of each exposure (e)
for e in range(len(texp)):
    for ch in range(2):
        for p in range(darkmodblack.shape[1]):
            darkblack[e,ch,p]=poly2d(texp[e],temperature[e,ch],darkmodblack[ch,p,:,:])
        for p in range(darkmodspec.shape[1]):
            darkspec[e,ch,p]=poly2d(texp[e],temperature[e,ch],darkmodspec[ch,p,:,:])
            spec[e,ch,p]=(spec[e,ch,p]-darkspec[e,ch,p])-np.mean(black[e,ch,:]-darkblack[e,ch,:])





#================ WAVELENGTH CALIBRATION
coeffs = np.asarray([ [355.688, 0.167436, -2.93242e-06, -2.22549e-10], \
          [360.071, 0.165454, -3.35036e-06, -1.88750e-10] ] )

#wavs = make_array(2, 3648)
#for c=0, 1 do wavs[c,*] = poly(findgen(3648),coefs[*,c])

from scipy.interpolate import interp1d
plt.plot(spec[0,0],'r')
plt.plot(spec[0,1],'b')


wavs=np.zeros((2,spec.shape[2]))
for ch in range(2): # in each spectral channel
    wavs[ch,:]=np.polyval(coeffs[ch,:][::-1],np.arange(spec.shape[2])) #reverse order of coeffs


for e in range(len(texp)):
    interpfunc=interp1d(wavs[:,0],spec[e,0,:],kind='linear')
    spec[e,0,:]=interpfunc(wavs[:,1])
    
plt.plot(spec[0,0],'r--')
plt.plot(spec[0,1],'b--')
plt.savefig('./temp.png')
