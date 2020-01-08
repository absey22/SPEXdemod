import matplotlib.pyplot as plt
import numpy as np

import sys
import glob


if len(sys.argv) == 1:
    dir='/home/seymour/Documents/SPEX/demodulation_pipeline/SPEXredux/ScriptRun_09072013_105035'
else:
    dir=sys.argv[1]

    
if str(dir)[-1]!='/':
    dir+='/'

    
#================ IMPORT META DATA

metafile=glob.glob(dir+"Meta*.*")

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
        pantilt[i-1,1]=float(line[0])

    if line[1].find('Detector 0 temp ')!=-1:
       temperature[i,0]=float(line[0])
    if line[1].find('Detector 1 temp ')!=-1:
        temperature[i-1,1]=float(line[0])
       
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

specfiles=glob.glob(dir+'Spectrometer_110516*U2_*pix.txt')

spec=[]
for i in range(len(specfiles)):
    spec.append(np.genfromtxt(specfiles[i],delimiter=',')[:-1])

spec=np.asarray(spec)

for s in spec:
    plt.plot(s)
plt.show()


#================ IMPORT DATE + TIMES

times=[]
for specfile in specfiles:
    t=specfile.split(' ')[1].split('_')[0]
    t=float(t[0:2])+float(t[2:4])/60.+float(t[4:])/3600.
    times.append(round(t,4))               
times=times[:len(texp)]

date=np.asarray([specfiles[0].split(' ')[0].split('_')[-1][:2], \
                 specfiles[0].split(' ')[0].split('_')[-1][2:4], \
                 specfiles[0].split(' ')[0].split('_')[-1][4:] ], dtype=float)


#================ IMPORT BLACK PIXELS

blackfiles=glob.glob(dir+'Spectrometer_110516*U2_*pix.dark13.txt')

black=[]
for i in range(len(blackfiles)):
    black.append(np.genfromtxt(blackfiles[i],delimiter=',')[:-1])

black=np.asarray(spec)

print(black)
