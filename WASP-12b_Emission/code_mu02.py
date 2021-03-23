#!/usr/bin/env python3
import sys
import numpy as np
from numpy import *
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.lines
from scipy import integrate
from scipy.integrate import simps
from matplotlib import pyplot
import emcee
import corner
#np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(threshold=False)
plt.rcParams['figure.figsize'] = [12,7]
#Data = pd.read_table('try.dat',sep = '\s+', dtype='unicode', header = None, index_col = None)
#print(Data)


################################################################################################################################################################################
####################################################################################### Source #################################################################################
################################################################################################################################################################################



######################################################### Wasp-12 (Source) @ 6300K ############################################################

myfile = open('Modified_lte063-4.5-0.0a+0.0.BT-Settl.spec.7_Intermediate', 'w')

with open("lte063-4.5-0.0a+0.0.BT-Settl.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 173:  
            lines = lines[0:173]+" "+lines[173:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()

myfile = open('Modified_lte063-4.5-0.0a+0.0.BT-Settl.spec.7', 'w')
with open("Modified_lte063-4.5-0.0a+0.0.BT-Settl.spec.7_Intermediate") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 128:  
            lines = lines[0:128]+" "+lines[128:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()

# Read in the file
with open('Modified_lte063-4.5-0.0a+0.0.BT-Settl.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_lte063-4.5-0.0a+0.0.BT-Settl.spec.7', 'w') as file:
  file.write(filedata)

Data = pd.read_table('Modified_lte063-4.5-0.0a+0.0.BT-Settl.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)


Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_6300 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_6300 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_6300 = 10**(flux_6300-8.)
#error_6300 = Spectra_Data[:,[16]]
#error_6300 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_6300).astype(float)
#Error_6300 = error_6300[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10) 
pyplot.yscale('log', base=10)
plt.plot(Wavelength_6300*(10**(-4)),Flux_6300,'y-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()





################################################################################################################################################################################
######################################################################################### mu=2 #################################################################################
################################################################################################################################################################################




######################################################### Wasp-12b (mu=2; Redistribution Factor=r2; C/O=0.27542) ############################################################



myfile = open('Modified_C_850_O_906_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("C_850_O_906_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 14:  
            lines = lines[0:14]+" "+lines[14:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()



# Read in the file
with open('Modified_C_850_O_906_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_C_850_O_906_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_C_850_O_906_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02_CO_0_27542_r2 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_mu02_CO_0_27542_r2 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02_CO_0_27542_r2 = 10**(flux_mu02_CO_0_27542_r2-8.)
#error_mu02_CO_0_27542_r2 = Spectra_Data[:,[16]]
#error_mu02_CO_0_27542_r2 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02_CO_0_27542_r2).astype(float)
#Error_mu02_CO_0_27542_r2 = error_mu02_CO_0_27542_r2[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10) 
pyplot.yscale('log', base=10)
plt.plot(Wavelength_mu02_CO_0_27542_r2*(10**(-4)), Flux_mu02_CO_0_27542_r2,'r-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()





######################################################### Wasp-12b (mu=2; Redistribution Factor=r2; C/O=0.34674) ############################################################



myfile = open('Modified_C_850_O_896_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("C_850_O_896_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 14:  
            lines = lines[0:14]+" "+lines[14:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()



# Read in the file
with open('Modified_C_850_O_896_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_C_850_O_896_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_C_850_O_896_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02_CO_0_34674_r2 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_mu02_CO_0_34674_r2 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02_CO_0_34674_r2 = 10**(flux_mu02_CO_0_34674_r2-8.)
#error_mu02_CO_0_34674_r2 = Spectra_Data[:,[16]]
#error_mu02_CO_0_34674_r2 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02_CO_0_34674_r2).astype(float)
#Error_mu02_CO_0_34674_r2 = error_mu02_CO_0_34674_r2[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10) 
pyplot.yscale('log', base=10)
plt.plot(Wavelength_mu02_CO_0_34674_r2*(10**(-4)), Flux_mu02_CO_0_34674_r2,'r-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()




######################################################### Wasp-12b (mu=2; Redistribution Factor=r2; C/O=0.43652) ############################################################




myfile = open('Modified_C_850_O_886_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("C_850_O_886_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 14:  
            lines = lines[0:14]+" "+lines[14:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()



# Read in the file
with open('Modified_C_850_O_886_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_C_850_O_886_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_C_850_O_886_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02_CO_0_43652_r2 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_mu02_CO_0_43652_r2 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02_CO_0_43652_r2 = 10**(flux_mu02_CO_0_43652_r2-8.)
#error_mu02_CO_0_43652_r2 = Spectra_Data[:,[16]]
#error_mu02_CO_0_43652_r2 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02_CO_0_43652_r2).astype(float)
#Error_mu02_CO_0_43652_r2 = error_mu02_CO_0_43652_r2[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10) 
pyplot.yscale('log', base=10)
plt.plot(Wavelength_mu02_CO_0_43652_r2*(10**(-4)), Flux_mu02_CO_0_43652_r2,'r-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()





######################################################### Wasp-12b (mu=2; Redistribution Factor=r2; C/O=0.5495) ############################################################



myfile = open('Modified_C_850_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("C_850_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 14:  
            lines = lines[0:14]+" "+lines[14:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()



# Read in the file
with open('Modified_C_850_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_C_850_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_C_850_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02_CO_0_5495_r2 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_mu02_CO_0_5495_r2 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02_CO_0_5495_r2 = 10**(flux_mu02_CO_0_5495_r2-8.)
#error_mu02_CO_0_5495_r2 = Spectra_Data[:,[16]]
#error_mu02_CO_0_5495_r2 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02_CO_0_5495_r2).astype(float)
#Error_mu02_CO_0_5495_r2 = error_mu02_CO_0_5495_r2[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10)
pyplot.yscale('log', base=10)
plt.plot(Wavelength_mu02_CO_0_5495_r2*(10**(-4)), Flux_mu02_CO_0_5495_r2,'r-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()






######################################################### Wasp-12b (mu=2; Redistribution Factor=r2; C/O=0.69183) ############################################################



myfile = open('Modified_C_860_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("C_860_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 14:  
            lines = lines[0:14]+" "+lines[14:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()



# Read in the file
with open('Modified_C_860_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_C_860_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_C_860_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02_CO_0_69183_r2 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_mu02_CO_0_69183_r2 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02_CO_0_69183_r2 = 10**(flux_mu02_CO_0_69183_r2-8.)
#error_mu02_CO_0_69183_r2 = Spectra_Data[:,[16]]
#error_mu02_CO_0_69183_r2 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02_CO_0_69183_r2).astype(float)
#Error_mu02_CO_0_69183_r2 = error_mu02_CO_0_69183_r2[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10)
pyplot.yscale('log', base=10)
plt.plot(Wavelength_mu02_CO_0_69183_r2*(10**(-4)), Flux_mu02_CO_0_69183_r2,'r-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()









######################################################### Wasp-12b (mu=2; Redistribution Factor=r2; C/O=0.8709) ############################################################



myfile = open('Modified_C_870_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("C_870_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 14:  
            lines = lines[0:14]+" "+lines[14:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()



# Read in the file
with open('Modified_C_870_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_C_870_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_C_870_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02_CO_0_8709_r2 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_mu02_CO_0_8709_r2 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02_CO_0_8709_r2 = 10**(flux_mu02_CO_0_8709_r2-8.)
#error_mu02_CO_0_8709_r2 = Spectra_Data[:,[16]]
#error_mu02_CO_0_8709_r2 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02_CO_0_8709_r2).astype(float)
#Error_mu02_CO_0_8709_r2 = error_mu02_CO_0_8709_r2[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10)
pyplot.yscale('log', base=10)
plt.plot(Wavelength_mu02_CO_0_8709_r2*(10**(-4)), Flux_mu02_CO_0_8709_r2,'r-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()





######################################################### Wasp-12b (mu=2; Redistribution Factor=r2; C/O=1.0964) ############################################################



myfile = open('Modified_C_880_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("C_880_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 14:  
            lines = lines[0:14]+" "+lines[14:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()



# Read in the file
with open('Modified_C_880_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_C_880_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_C_880_O_876_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02_CO_1_0964_r2 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_mu02_CO_1_0964_r2 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02_CO_1_0964_r2 = 10**(flux_mu02_CO_1_0964_r2-8.)
#error_mu02_CO_1_0964_r2 = Spectra_Data[:,[16]]
#error_mu02_CO_1_0964_r2 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02_CO_1_0964_r2).astype(float)
#Error_mu02_CO_1_0964_r2 = error_mu02_CO_1_0964_r2[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10)
pyplot.yscale('log', base=10)
plt.plot(Wavelength_mu02_CO_1_0964_r2*(10**(-4)), Flux_mu02_CO_1_0964_r2,'r-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()





######################################################### Wasp-12b (mu=2; Redistribution Factor=r2; C/O=1.51356) ############################################################



myfile = open('Modified_C_860_O_842_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("C_860_O_842_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
    for lines in data.readlines():
        lines = lines.rstrip()
        if len(lines) > 14:  
            lines = lines[0:14]+" "+lines[14:]
            myfile.write(lines)
            myfile.write("\n")
        else:
            lines = lines
            myfile.write(lines)
            myfile.write("\n")
myfile.close()



# Read in the file
with open('Modified_C_860_O_842_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_C_860_O_842_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_C_860_O_842_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02_CO_1_51356_r2 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
flux_mu02_CO_1_51356_r2 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02_CO_1_51356_r2 = 10**(flux_mu02_CO_1_51356_r2-8.)
#error_mu02_CO_1_51356_r2 = Spectra_Data[:,[16]]
#error_mu02_CO_1_51356_r2 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02_CO_1_51356_r2).astype(float)
#Error_mu02_CO_1_51356_r2 = error_mu02_CO_1_51356_r2[wavelength1[0][0]:wavelength2[0][0],[0]]

pyplot.xscale('log', base=10)
pyplot.yscale('log', base=10)
plt.plot(Wavelength_mu02_CO_1_51356_r2*(10**(-4)), Flux_mu02_CO_1_51356_r2,'r-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()





########################################################################################################################
#################################################### Source Interpolation ##############################################
########################################################################################################################

Rs = 1.630*(696340*((10)**3))
Rp = 1.790*(69911*((10)**3))
Ts = 6300
h = 6.62607*((10)**-34)
c = 299792458
k = 1.3806*((10)**-23)


Wavelength_Source_Interp = Wavelength_mu02_CO_0_27542_r2.flatten()

Wavelength_6300 = Wavelength_6300.flatten()
Flux_6300 = Flux_6300.flatten()
Flux_Source_Interp = np.interp(Wavelength_Source_Interp, Wavelength_6300, Flux_6300)

#Error_6300 = Error_6300.flatten()
#Error_Source_Interp = np.interp(Wavelength_Source_Interp,Wavelength_6300, Error_6300)
#Error_mu02_CO_0_5495_r2_Flatten = Error_mu02_CO_0_5495_r2.flatten()

pyplot.xscale('log', base=10)
pyplot.yscale('log', base=10)
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Source_Interp,'c-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()


#################################################################################
############################## Flux Ratio #######################################
#################################################################################


Flux_mu02_CO_0_27542_r2_Flatten = Flux_mu02_CO_0_27542_r2.flatten()
Flux_Ratio_mu02_CO_0_27542_r2 = ((Flux_mu02_CO_0_27542_r2_Flatten/Flux_Source_Interp)*((Rp/Rs)**2))
pyplot.xscale('log', base=10) 
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Ratio_mu02_CO_0_27542_r2,color = '#5c8281')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()


Flux_mu02_CO_0_34674_r2_Flatten = Flux_mu02_CO_0_34674_r2.flatten()
Flux_Ratio_mu02_CO_0_34674_r2 = ((Flux_mu02_CO_0_34674_r2_Flatten/Flux_Source_Interp)*((Rp/Rs)**2))
pyplot.xscale('log', base=10) 
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Ratio_mu02_CO_0_34674_r2,color = '#5c8281')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()


Flux_mu02_CO_0_43652_r2_Flatten = Flux_mu02_CO_0_43652_r2.flatten()
Flux_Ratio_mu02_CO_0_43652_r2 = ((Flux_mu02_CO_0_43652_r2_Flatten/Flux_Source_Interp)*((Rp/Rs)**2))
pyplot.xscale('log', base=10) 
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Ratio_mu02_CO_0_43652_r2,color = '#5c8281')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()


Flux_mu02_CO_0_5495_r2_Flatten = Flux_mu02_CO_0_5495_r2.flatten()
Flux_Ratio_mu02_CO_0_5495_r2 = ((Flux_mu02_CO_0_5495_r2_Flatten/Flux_Source_Interp)*((Rp/Rs)**2))
pyplot.xscale('log', base=10)
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Ratio_mu02_CO_0_5495_r2,color = '#5c8281')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()


Flux_mu02_CO_0_69183_r2_Flatten = Flux_mu02_CO_0_69183_r2.flatten()
Flux_Ratio_mu02_CO_0_69183_r2 = ((Flux_mu02_CO_0_69183_r2_Flatten/Flux_Source_Interp)*((Rp/Rs)**2))
pyplot.xscale('log', base=10)
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Ratio_mu02_CO_0_69183_r2,color = '#5c8281')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()


Flux_mu02_CO_0_8709_r2_Flatten = Flux_mu02_CO_0_8709_r2.flatten()
Flux_Ratio_mu02_CO_0_8709_r2 = ((Flux_mu02_CO_0_8709_r2_Flatten/Flux_Source_Interp)*((Rp/Rs)**2))
pyplot.xscale('log', base=10)
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Ratio_mu02_CO_0_8709_r2,color = '#5c8281')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()


Flux_mu02_CO_1_0964_r2_Flatten = Flux_mu02_CO_1_0964_r2.flatten()
Flux_Ratio_mu02_CO_1_0964_r2 = ((Flux_mu02_CO_1_0964_r2_Flatten/Flux_Source_Interp)*((Rp/Rs)**2))
pyplot.xscale('log', base=10)
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Ratio_mu02_CO_1_0964_r2,color = '#5c8281')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()


Flux_mu02_CO_1_51356_r2_Flatten = Flux_mu02_CO_1_51356_r2.flatten()
Flux_Ratio_mu02_CO_1_51356_r2 = ((Flux_mu02_CO_1_51356_r2_Flatten/Flux_Source_Interp)*((Rp/Rs)**2))
pyplot.xscale('log', base=10)
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Source_Interp*(10**(-4)), Flux_Ratio_mu02_CO_1_51356_r2,color = '#5c8281')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()






##################################################################################################################
#################################################### Interpolation ###############################################
##################################################################################################################



x = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
results = []
for i in range(0, len(Wavelength_Source_Interp)):
    xp = [10**(8.50-9.06), 10**(8.50-8.96), 10**(8.50-8.86), 10**(8.50-8.76), 10**(8.60-8.76), 10**(8.70-8.76), 10**(8.80-8.76), 10**(8.60-8.42)]
    fp = [Flux_Ratio_mu02_CO_0_27542_r2[i], Flux_Ratio_mu02_CO_0_34674_r2[i], Flux_Ratio_mu02_CO_0_43652_r2[i], Flux_Ratio_mu02_CO_0_5495_r2[i], Flux_Ratio_mu02_CO_0_69183_r2[i], Flux_Ratio_mu02_CO_0_8709_r2[i], Flux_Ratio_mu02_CO_1_0964_r2[i], Flux_Ratio_mu02_CO_1_51356_r2[i]]
    f = np.interp(x, xp, fp)
    results.append(f)

Flux_Ratio_Table_mu02 = np.stack(results, axis=0) # For the C/O ratio given in x
Flux_Ratio_Table_mu02 = np.array([Flux_Ratio_mu02_CO_0_27542_r2, Flux_Ratio_Table_mu02[:,0], Flux_Ratio_mu02_CO_0_34674_r2, Flux_Ratio_Table_mu02[:,1], Flux_Ratio_mu02_CO_0_43652_r2, Flux_Ratio_Table_mu02[:,2], Flux_Ratio_mu02_CO_0_5495_r2, Flux_Ratio_Table_mu02[:,3], Flux_Ratio_mu02_CO_0_69183_r2, Flux_Ratio_Table_mu02[:,4], Flux_Ratio_Table_mu02[:,5], Flux_Ratio_mu02_CO_0_8709_r2, Flux_Ratio_Table_mu02[:,6], Flux_Ratio_Table_mu02[:,7], Flux_Ratio_mu02_CO_1_0964_r2, Flux_Ratio_Table_mu02[:,8], Flux_Ratio_Table_mu02[:,9], Flux_Ratio_Table_mu02[:,10], Flux_Ratio_Table_mu02[:,11], Flux_Ratio_Table_mu02[:,12], Flux_Ratio_mu02_CO_1_51356_r2]).T


results2 = []
for i in range(0, len(Wavelength_Source_Interp)):
    xp = [10**(8.50-9.06), 10**(8.50-8.96), 10**(8.50-8.86), 10**(8.50-8.76), 10**(8.60-8.76), 10**(8.70-8.76), 10**(8.80-8.76), 10**(8.60-8.42)]
    fp = [Flux_mu02_CO_0_27542_r2[i].item(), Flux_mu02_CO_0_34674_r2[i].item(), Flux_mu02_CO_0_43652_r2[i].item(), Flux_mu02_CO_0_5495_r2[i].item(), Flux_mu02_CO_0_69183_r2[i].item(), Flux_mu02_CO_0_8709_r2[i].item(), Flux_mu02_CO_1_0964_r2[i].item(), Flux_mu02_CO_1_51356_r2[i].item()]
    f = np.interp(x, xp, fp)
    results2.append(f)

Planet_Flux_Table_mu02 = np.stack(results2, axis=0) # For the C/O ratio given in x
Planet_Flux_Table_mu02 = np.array([Flux_mu02_CO_0_27542_r2.flatten(), Planet_Flux_Table_mu02[:,0], Flux_mu02_CO_0_34674_r2.flatten(), Planet_Flux_Table_mu02[:,1], Flux_mu02_CO_0_43652_r2.flatten(), Planet_Flux_Table_mu02[:,2], Flux_mu02_CO_0_5495_r2.flatten(), Planet_Flux_Table_mu02[:,3], Flux_mu02_CO_0_69183_r2.flatten(), Planet_Flux_Table_mu02[:,4], Planet_Flux_Table_mu02[:,5], Flux_mu02_CO_0_8709_r2.flatten(), Planet_Flux_Table_mu02[:,6], Planet_Flux_Table_mu02[:,7], Flux_mu02_CO_1_0964_r2.flatten(), Planet_Flux_Table_mu02[:,8], Planet_Flux_Table_mu02[:,9], Planet_Flux_Table_mu02[:,10], Planet_Flux_Table_mu02[:,11], Planet_Flux_Table_mu02[:,12], Flux_mu02_CO_1_51356_r2.flatten()]).T

X = np.sort(np.append(x,xp))

##################################################################################################################
###################################################### Smoothing #################################################
##################################################################################################################



def smooth(Flux_Ratio, box_pts):
    box = np.ones(box_pts)/box_pts
    Flux_Ratio_Smooth = np.convolve(Flux_Ratio, box, mode='full')
    return Flux_Ratio_Smooth


m=100
n=250
o=500

#################################################### Box Width = 100 ##############################################

Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,0],m)
Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_100_Full[m-1:-(m-1)]
wavelength_source_interp_100_trim = Wavelength_Source_Interp[m-1:]


Flux_Ratio_mu02_CO_0_3_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,1],m)
Flux_Ratio_mu02_CO_0_3_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_3_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,2],m)
Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_4_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,3],m)
Flux_Ratio_mu02_CO_0_4_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_4_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,4],m)
Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_100_Full[m-1:-(m-1)]

Flux_Ratio_mu02_CO_0_5_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,5],m)
Flux_Ratio_mu02_CO_0_5_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_5_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,6],m)
Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_6_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,7],m)
Flux_Ratio_mu02_CO_0_6_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_6_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,8],m)
Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_7_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,9],m)
Flux_Ratio_mu02_CO_0_7_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_7_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_8_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,10],m)
Flux_Ratio_mu02_CO_0_8_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_8_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,11],m)
Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_0_9_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,12],m)
Flux_Ratio_mu02_CO_0_9_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_0_9_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_1_0_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,13],m)
Flux_Ratio_mu02_CO_1_0_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_1_0_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,14],m)
Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_1_1_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,15],m)
Flux_Ratio_mu02_CO_1_1_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_1_1_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_1_2_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,16],m)
Flux_Ratio_mu02_CO_1_2_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_1_2_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_1_3_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,17],m)
Flux_Ratio_mu02_CO_1_3_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_1_3_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_1_4_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,18],m)
Flux_Ratio_mu02_CO_1_4_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_1_4_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_1_5_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,19],m)
Flux_Ratio_mu02_CO_1_5_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_1_5_r2_Boxcar_100_Full[m-1:-(m-1)]


Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_100_Full = smooth(Flux_Ratio_Table_mu02[:,20],m)
Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_100_Trim = Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_100_Full[m-1:-(m-1)]




#################################################### Box Width = 250 ##############################################






Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,0],n)
Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_250_Full[n-1:-(n-1)]
wavelength_source_interp_250_trim = Wavelength_Source_Interp[n-1:]


Flux_Ratio_mu02_CO_0_3_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,1],n)
Flux_Ratio_mu02_CO_0_3_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_3_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,2],n)
Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_4_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,3],n)
Flux_Ratio_mu02_CO_0_4_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_4_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,4],n)
Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_250_Full[n-1:-(n-1)]

Flux_Ratio_mu02_CO_0_5_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,5],n)
Flux_Ratio_mu02_CO_0_5_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_5_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,6],n)
Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_6_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,7],n)
Flux_Ratio_mu02_CO_0_6_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_6_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,8],n)
Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_7_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,9],n)
Flux_Ratio_mu02_CO_0_7_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_7_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_8_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,10],n)
Flux_Ratio_mu02_CO_0_8_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_8_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,11],n)
Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_0_9_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,12],n)
Flux_Ratio_mu02_CO_0_9_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_0_9_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_1_0_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,13],n)
Flux_Ratio_mu02_CO_1_0_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_1_0_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,14],n)
Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_1_1_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,15],n)
Flux_Ratio_mu02_CO_1_1_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_1_1_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_1_2_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,16],n)
Flux_Ratio_mu02_CO_1_2_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_1_2_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_1_3_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,17],n)
Flux_Ratio_mu02_CO_1_3_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_1_3_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_1_4_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,18],n)
Flux_Ratio_mu02_CO_1_4_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_1_4_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_1_5_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,19],n)
Flux_Ratio_mu02_CO_1_5_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_1_5_r2_Boxcar_250_Full[n-1:-(n-1)]


Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_250_Full = smooth(Flux_Ratio_Table_mu02[:,20],n)
Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_250_Trim = Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_250_Full[n-1:-(n-1)]







#################################################### Box Width = 500 ##############################################



Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,0],o)
Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_500_Full[o-1:-(o-1)]
wavelength_source_interp_500_trim = Wavelength_Source_Interp[o-1:]


Flux_Ratio_mu02_CO_0_3_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,1],o)
Flux_Ratio_mu02_CO_0_3_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_3_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,2],o)
Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_4_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,3],o)
Flux_Ratio_mu02_CO_0_4_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_4_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,4],o)
Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_500_Full[o-1:-(o-1)]

Flux_Ratio_mu02_CO_0_5_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,5],o)
Flux_Ratio_mu02_CO_0_5_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_5_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,6],o)
Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_6_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,7],o)
Flux_Ratio_mu02_CO_0_6_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_6_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,8],o)
Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_7_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,9],o)
Flux_Ratio_mu02_CO_0_7_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_7_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_8_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,10],o)
Flux_Ratio_mu02_CO_0_8_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_8_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,11],o)
Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_0_9_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,12],o)
Flux_Ratio_mu02_CO_0_9_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_0_9_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_1_0_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,13],o)
Flux_Ratio_mu02_CO_1_0_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_1_0_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,14],o)
Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_1_1_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,15],o)
Flux_Ratio_mu02_CO_1_1_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_1_1_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_1_2_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,16],o)
Flux_Ratio_mu02_CO_1_2_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_1_2_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_1_3_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,17],o)
Flux_Ratio_mu02_CO_1_3_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_1_3_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_1_4_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,18],o)
Flux_Ratio_mu02_CO_1_4_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_1_4_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_1_5_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,19],o)
Flux_Ratio_mu02_CO_1_5_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_1_5_r2_Boxcar_500_Full[o-1:-(o-1)]


Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_500_Full = smooth(Flux_Ratio_Table_mu02[:,20],o)
Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_500_Trim = Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_500_Full[o-1:-(o-1)]






#################################################################################
############################### Smoothing Output ################################
#################################################################################


Flux_Ratio_mu02_Boxcar_100_Trim = np.array([Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_3_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_4_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_5_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_6_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_7_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_8_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_0_9_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_1_0_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_1_1_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_1_2_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_1_3_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_1_4_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_1_5_r2_Boxcar_100_Trim, Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_100_Trim]).T


Flux_Ratio_mu02_Boxcar_250_Trim = np.array([Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_3_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_4_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_5_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_6_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_7_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_8_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_0_9_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_1_0_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_1_1_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_1_2_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_1_3_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_1_4_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_1_5_r2_Boxcar_250_Trim, Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_250_Trim]).T


Flux_Ratio_mu02_Boxcar_500_Trim = np.array([Flux_Ratio_mu02_CO_0_27542_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_3_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_34674_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_4_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_43652_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_5_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_5495_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_6_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_69183_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_7_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_8_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_8709_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_0_9_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_1_0_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_1_0964_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_1_1_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_1_2_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_1_3_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_1_4_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_1_5_r2_Boxcar_500_Trim, Flux_Ratio_mu02_CO_1_51356_r2_Boxcar_500_Trim]).T


Wavelength_Source_Interp_100_Trim = np.tile(wavelength_source_interp_100_trim,(len(X),1)).T

Wavelength_Source_Interp_250_Trim = np.tile(wavelength_source_interp_250_trim,(len(X),1)).T

Wavelength_Source_Interp_500_Trim = np.tile(wavelength_source_interp_500_trim,(len(X),1)).T



##################################################################################################################################################################
############################################# Photometry in Spitzer Bands (Using Definition from Charbonneau et al., 2005) #######################################
##################################################################################################################################################################



################################################################################
######################### IRAC/SPITZER Transmission ############################
################################################################################

IRAC_CH_1 = pd.read_table('201125ch1trans_full.txt',sep = '\s+', dtype='unicode', header = None, index_col = None)
IRAC_CH_2 = pd.read_table('201125ch2trans_full.txt',sep = '\s+', dtype='unicode', header = None, index_col = None)
IRAC_CH_3 = pd.read_table('201125ch3trans_full.txt',sep = '\s+', dtype='unicode', header = None, index_col = None)
IRAC_CH_4 = pd.read_table('201125ch4trans_full.txt',sep = '\s+', dtype='unicode', header = None, index_col = None)


# Band-1 Wavelengths (Channel Centred at 3.6 um)

CH_1_L = 30810.60 #In Angstrom
CH_1_U = 40103.80 #In Angstrom


# Band-2 Wavelengths (Channel Centred at 4.5 um)

CH_2_L = 37225.0#37224.90
CH_2_U = 52219.80


# Band-3 Wavelengths (Channel Centred at 5.8 um)

CH_3_L = 47442.20#47442.10
CH_3_U = 66225.00#66225.10


# Band-4 Wavelengths (Channel Centred at 8.5 um)

CH_4_L = 61511.50
CH_4_U = 97287.50



###########################################################################
######################### Channel-1 Photometry ############################
###########################################################################

Index1_Channel_1 = np.where(Wavelength_Source_Interp == CH_1_L)
Index2_Channel_1 = np.where(Wavelength_Source_Interp == CH_1_U)

IRAC_CH_1_Wavelength_Source_Interp = Wavelength_Source_Interp[Index1_Channel_1[0][0]:Index2_Channel_1[0][0]]
IRAC_CH_1_wavelength_center = 3.6000 #In microns

IRAC_CH_1_Data = np.array(IRAC_CH_1)

IRAC_CH_1_wavelength = IRAC_CH_1_Data[:,[0]]
IRAC_CH_1_transmission = IRAC_CH_1_Data[:,[1]]
IRAC_CH_1_Transmission = np.frompyfunc(lambda x: x.replace(',',''),1,1)(IRAC_CH_1_transmission).astype(float)
IRAC_CH_1_Wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(IRAC_CH_1_wavelength).astype(float)*(10**(4))

IRAC_CH_1_Wavelength = IRAC_CH_1_Wavelength.flatten()
IRAC_CH_1_Transmission = IRAC_CH_1_Transmission.flatten()
IRAC_CH_1_Transmission_Interp = np.interp(IRAC_CH_1_Wavelength_Source_Interp, IRAC_CH_1_Wavelength, IRAC_CH_1_Transmission)


plt.plot(IRAC_CH_1_Wavelength_Source_Interp*(10**(-4)), IRAC_CH_1_Transmission_Interp,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

IRAC_CH_1_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Channel_1[0][0]:Index2_Channel_1[0][0]]
IRAC_CH_1_Flux_Source_Interp = Flux_Source_Interp[Index1_Channel_1[0][0]:Index2_Channel_1[0][0]]


f2_CH_1 = (IRAC_CH_1_Flux_Source_Interp)*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)


f1_CH_1_mu02_CO_0_27542_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,0])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_27542_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_3_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,1])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_3_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_34674_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,2])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_34674_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_4_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,3])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_4_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_43652_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,4])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_43652_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_5_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,5])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_5_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_5495_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,6])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_5495_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_6_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,7])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_6_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_69183_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,8])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_69183_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_7_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,9])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_7_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_8_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,10])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_8_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_8709_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,11])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_8709_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_0_9_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,12])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_0_9_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_1_0_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,13])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_1_0_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_1_0964_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,14])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_1_0964_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_1_1_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,15])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_1_1_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_1_2_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,16])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_1_2_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_1_3_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,17])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_1_3_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_1_4_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,18])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_1_4_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_1_5_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,19])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_1_5_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))

f1_CH_1_mu02_CO_1_51356_r2 = (IRAC_CH_1_Planet_Flux_Table_mu02[:,20])*(IRAC_CH_1_Transmission_Interp)*(IRAC_CH_1_Wavelength_Source_Interp)/(h*c)
IRAC_CH_1_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_1_mu02_CO_1_51356_r2, IRAC_CH_1_Wavelength_Source_Interp))/(simps(f2_CH_1, IRAC_CH_1_Wavelength_Source_Interp))






###########################################################################
######################### Channel-2 Photometry ############################
###########################################################################

Index1_Channel_2 = np.where(Wavelength_Source_Interp == CH_2_L)
Index2_Channel_2 = np.where(Wavelength_Source_Interp == CH_2_U)

IRAC_CH_2_Wavelength_Source_Interp = Wavelength_Source_Interp[Index1_Channel_2[0][0]:Index2_Channel_2[0][0]]
IRAC_CH_2_wavelength_center = 4.5000 #In microns

IRAC_CH_2_Data = np.array(IRAC_CH_2)

IRAC_CH_2_wavelength = IRAC_CH_2_Data[:,[0]]
IRAC_CH_2_transmission = IRAC_CH_2_Data[:,[1]]
IRAC_CH_2_Transmission = np.frompyfunc(lambda x: x.replace(',',''),1,1)(IRAC_CH_2_transmission).astype(float)
IRAC_CH_2_Wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(IRAC_CH_2_wavelength).astype(float)*(10**(4))

IRAC_CH_2_Wavelength = IRAC_CH_2_Wavelength.flatten()
IRAC_CH_2_Transmission = IRAC_CH_2_Transmission.flatten()
IRAC_CH_2_Transmission_Interp = np.interp(IRAC_CH_2_Wavelength_Source_Interp, IRAC_CH_2_Wavelength, IRAC_CH_2_Transmission)


plt.plot(IRAC_CH_2_Wavelength_Source_Interp*(10**(-4)), IRAC_CH_2_Transmission_Interp,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

IRAC_CH_2_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Channel_2[0][0]:Index2_Channel_2[0][0]]
IRAC_CH_2_Flux_Source_Interp = Flux_Source_Interp[Index1_Channel_2[0][0]:Index2_Channel_2[0][0]]


f2_CH_2 = (IRAC_CH_2_Flux_Source_Interp)*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)



f1_CH_2_mu02_CO_0_27542_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,0])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_27542_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_3_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,1])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_3_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_34674_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,2])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_34674_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_4_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,3])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_4_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_43652_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,4])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_43652_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_5_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,5])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_5_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_5495_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,6])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_5495_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_6_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,7])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_6_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_69183_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,8])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_69183_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_7_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,9])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_7_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_8_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,10])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_8_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_8709_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,11])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_8709_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_0_9_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,12])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_0_9_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_1_0_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,13])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_1_0_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_1_0964_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,14])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_1_0964_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_1_1_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,15])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_1_1_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_1_2_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,16])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_1_2_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_1_3_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,17])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_1_3_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_1_4_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,18])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_1_4_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_1_5_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,19])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_1_5_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))

f1_CH_2_mu02_CO_1_51356_r2 = (IRAC_CH_2_Planet_Flux_Table_mu02[:,20])*(IRAC_CH_2_Transmission_Interp)*(IRAC_CH_2_Wavelength_Source_Interp)/(h*c)
IRAC_CH_2_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_2_mu02_CO_1_51356_r2, IRAC_CH_2_Wavelength_Source_Interp))/(simps(f2_CH_2, IRAC_CH_2_Wavelength_Source_Interp))





###########################################################################
######################### Channel-3 Photometry ############################
###########################################################################

Index1_Channel_3 = np.where(Wavelength_Source_Interp == CH_3_L)
Index2_Channel_3 = np.where(Wavelength_Source_Interp == CH_3_U)

IRAC_CH_3_Wavelength_Source_Interp = Wavelength_Source_Interp[Index1_Channel_3[0][0]:Index2_Channel_3[0][0]]
IRAC_CH_3_wavelength_center = 5.8000 #In microns

IRAC_CH_3_Data = np.array(IRAC_CH_3)

IRAC_CH_3_wavelength = IRAC_CH_3_Data[:,[0]]
IRAC_CH_3_transmission = IRAC_CH_3_Data[:,[1]]
IRAC_CH_3_Transmission = np.frompyfunc(lambda x: x.replace(',',''),1,1)(IRAC_CH_3_transmission).astype(float)
IRAC_CH_3_Wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(IRAC_CH_3_wavelength).astype(float)*(10**(4))

IRAC_CH_3_Wavelength = IRAC_CH_3_Wavelength.flatten()
IRAC_CH_3_Transmission = IRAC_CH_3_Transmission.flatten()
IRAC_CH_3_Transmission_Interp = np.interp(IRAC_CH_3_Wavelength_Source_Interp, IRAC_CH_3_Wavelength, IRAC_CH_3_Transmission)


plt.plot(IRAC_CH_3_Wavelength_Source_Interp*(10**(-4)), IRAC_CH_3_Transmission_Interp,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

IRAC_CH_3_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Channel_3[0][0]:Index2_Channel_3[0][0]]
IRAC_CH_3_Flux_Source_Interp = Flux_Source_Interp[Index1_Channel_3[0][0]:Index2_Channel_3[0][0]]


f2_CH_3 = (IRAC_CH_3_Flux_Source_Interp)*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)




f1_CH_3_mu02_CO_0_27542_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,0])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_27542_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_3_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,1])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_3_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_34674_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,2])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_34674_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_4_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,3])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_4_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_43652_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,4])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_43652_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_5_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,5])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_5_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_5495_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,6])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_5495_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_6_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,7])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_6_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_69183_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,8])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_69183_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_7_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,9])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_7_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_8_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,10])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_8_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_8709_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,11])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_8709_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_0_9_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,12])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_0_9_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_1_0_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,13])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_1_0_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_1_0964_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,14])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_1_0964_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_1_1_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,15])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_1_1_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_1_2_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,16])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_1_2_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_1_3_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,17])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_1_3_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_1_4_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,18])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_1_4_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_1_5_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,19])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_1_5_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))

f1_CH_3_mu02_CO_1_51356_r2 = (IRAC_CH_3_Planet_Flux_Table_mu02[:,20])*(IRAC_CH_3_Transmission_Interp)*(IRAC_CH_3_Wavelength_Source_Interp)/(h*c)
IRAC_CH_3_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_3_mu02_CO_1_51356_r2, IRAC_CH_3_Wavelength_Source_Interp))/(simps(f2_CH_3, IRAC_CH_3_Wavelength_Source_Interp))





###########################################################################
######################### Channel-4 Photometry ############################
###########################################################################

Index1_Channel_4 = np.where(Wavelength_Source_Interp == CH_4_L)
Index2_Channel_4 = np.where(Wavelength_Source_Interp == CH_4_U)

IRAC_CH_4_Wavelength_Source_Interp = Wavelength_Source_Interp[Index1_Channel_4[0][0]:Index2_Channel_4[0][0]]
IRAC_CH_4_wavelength_center = 8.0000 #In microns

IRAC_CH_4_Data = np.array(IRAC_CH_4)

IRAC_CH_4_wavelength = IRAC_CH_4_Data[:,[0]]
IRAC_CH_4_transmission = IRAC_CH_4_Data[:,[1]]
IRAC_CH_4_Transmission = np.frompyfunc(lambda x: x.replace(',',''),1,1)(IRAC_CH_4_transmission).astype(float)
IRAC_CH_4_Wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(IRAC_CH_4_wavelength).astype(float)*(10**(4))

IRAC_CH_4_Wavelength = IRAC_CH_4_Wavelength.flatten()
IRAC_CH_4_Transmission = IRAC_CH_4_Transmission.flatten()
IRAC_CH_4_Transmission_Interp = np.interp(IRAC_CH_4_Wavelength_Source_Interp, IRAC_CH_4_Wavelength, IRAC_CH_4_Transmission)


plt.plot(IRAC_CH_4_Wavelength_Source_Interp*(10**(-4)), IRAC_CH_4_Transmission_Interp,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

IRAC_CH_4_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Channel_4[0][0]:Index2_Channel_4[0][0]]
IRAC_CH_4_Flux_Source_Interp = Flux_Source_Interp[Index1_Channel_4[0][0]:Index2_Channel_4[0][0]]


f2_CH_4 = (IRAC_CH_4_Flux_Source_Interp)*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)





f1_CH_4_mu02_CO_0_27542_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,0])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_27542_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_3_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,1])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_3_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_34674_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,2])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_34674_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_4_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,3])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_4_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_43652_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,4])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_43652_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_5_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,5])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_5_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_5495_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,6])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_5495_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_6_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,7])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_6_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_69183_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,8])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_69183_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_7_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,9])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_7_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_8_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,10])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_8_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_8709_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,11])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_8709_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_0_9_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,12])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_0_9_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_1_0_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,13])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_1_0_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_1_0964_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,14])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_1_0964_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_1_1_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,15])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_1_1_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_1_2_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,16])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_1_2_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_1_3_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,17])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_1_3_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_1_4_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,18])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_1_4_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_1_5_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,19])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_1_5_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))

f1_CH_4_mu02_CO_1_51356_r2 = (IRAC_CH_4_Planet_Flux_Table_mu02[:,20])*(IRAC_CH_4_Transmission_Interp)*(IRAC_CH_4_Wavelength_Source_Interp)/(h*c)
IRAC_CH_4_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_CH_4_mu02_CO_1_51356_r2, IRAC_CH_4_Wavelength_Source_Interp))/(simps(f2_CH_4, IRAC_CH_4_Wavelength_Source_Interp))





#################################################################################
######################### Spitzer Photometric Output ############################
#################################################################################

IRAC_CH_1_mu02_Photometry = np.array([IRAC_CH_1_mu02_CO_0_27542_r2_Photometry, IRAC_CH_1_mu02_CO_0_3_r2_Photometry, IRAC_CH_1_mu02_CO_0_34674_r2_Photometry, IRAC_CH_1_mu02_CO_0_4_r2_Photometry, IRAC_CH_1_mu02_CO_0_43652_r2_Photometry, IRAC_CH_1_mu02_CO_0_5_r2_Photometry, IRAC_CH_1_mu02_CO_0_5495_r2_Photometry, IRAC_CH_1_mu02_CO_0_6_r2_Photometry, IRAC_CH_1_mu02_CO_0_69183_r2_Photometry, IRAC_CH_1_mu02_CO_0_7_r2_Photometry, IRAC_CH_1_mu02_CO_0_8_r2_Photometry, IRAC_CH_1_mu02_CO_0_8709_r2_Photometry, IRAC_CH_1_mu02_CO_0_9_r2_Photometry, IRAC_CH_1_mu02_CO_1_0_r2_Photometry, IRAC_CH_1_mu02_CO_1_0964_r2_Photometry, IRAC_CH_1_mu02_CO_1_1_r2_Photometry, IRAC_CH_1_mu02_CO_1_2_r2_Photometry, IRAC_CH_1_mu02_CO_1_3_r2_Photometry, IRAC_CH_1_mu02_CO_1_4_r2_Photometry, IRAC_CH_1_mu02_CO_1_5_r2_Photometry, IRAC_CH_1_mu02_CO_1_51356_r2_Photometry]);
IRAC_CH_2_mu02_Photometry = np.array([IRAC_CH_2_mu02_CO_0_27542_r2_Photometry, IRAC_CH_2_mu02_CO_0_3_r2_Photometry, IRAC_CH_2_mu02_CO_0_34674_r2_Photometry, IRAC_CH_2_mu02_CO_0_4_r2_Photometry, IRAC_CH_2_mu02_CO_0_43652_r2_Photometry, IRAC_CH_2_mu02_CO_0_5_r2_Photometry, IRAC_CH_2_mu02_CO_0_5495_r2_Photometry, IRAC_CH_2_mu02_CO_0_6_r2_Photometry, IRAC_CH_2_mu02_CO_0_69183_r2_Photometry, IRAC_CH_2_mu02_CO_0_7_r2_Photometry, IRAC_CH_2_mu02_CO_0_8_r2_Photometry, IRAC_CH_2_mu02_CO_0_8709_r2_Photometry, IRAC_CH_2_mu02_CO_0_9_r2_Photometry, IRAC_CH_2_mu02_CO_1_0_r2_Photometry, IRAC_CH_2_mu02_CO_1_0964_r2_Photometry, IRAC_CH_2_mu02_CO_1_1_r2_Photometry, IRAC_CH_2_mu02_CO_1_2_r2_Photometry, IRAC_CH_2_mu02_CO_1_3_r2_Photometry, IRAC_CH_2_mu02_CO_1_4_r2_Photometry, IRAC_CH_2_mu02_CO_1_5_r2_Photometry, IRAC_CH_2_mu02_CO_1_51356_r2_Photometry]);
IRAC_CH_3_mu02_Photometry = np.array([IRAC_CH_3_mu02_CO_0_27542_r2_Photometry, IRAC_CH_3_mu02_CO_0_3_r2_Photometry, IRAC_CH_3_mu02_CO_0_34674_r2_Photometry, IRAC_CH_3_mu02_CO_0_4_r2_Photometry, IRAC_CH_3_mu02_CO_0_43652_r2_Photometry, IRAC_CH_3_mu02_CO_0_5_r2_Photometry, IRAC_CH_3_mu02_CO_0_5495_r2_Photometry, IRAC_CH_3_mu02_CO_0_6_r2_Photometry, IRAC_CH_3_mu02_CO_0_69183_r2_Photometry, IRAC_CH_3_mu02_CO_0_7_r2_Photometry, IRAC_CH_3_mu02_CO_0_8_r2_Photometry, IRAC_CH_3_mu02_CO_0_8709_r2_Photometry, IRAC_CH_3_mu02_CO_0_9_r2_Photometry, IRAC_CH_3_mu02_CO_1_0_r2_Photometry, IRAC_CH_3_mu02_CO_1_0964_r2_Photometry, IRAC_CH_3_mu02_CO_1_1_r2_Photometry, IRAC_CH_3_mu02_CO_1_2_r2_Photometry, IRAC_CH_3_mu02_CO_1_3_r2_Photometry, IRAC_CH_3_mu02_CO_1_4_r2_Photometry, IRAC_CH_3_mu02_CO_1_5_r2_Photometry, IRAC_CH_3_mu02_CO_1_51356_r2_Photometry]);
IRAC_CH_4_mu02_Photometry = np.array([IRAC_CH_4_mu02_CO_0_27542_r2_Photometry, IRAC_CH_4_mu02_CO_0_3_r2_Photometry, IRAC_CH_4_mu02_CO_0_34674_r2_Photometry, IRAC_CH_4_mu02_CO_0_4_r2_Photometry, IRAC_CH_4_mu02_CO_0_43652_r2_Photometry, IRAC_CH_4_mu02_CO_0_5_r2_Photometry, IRAC_CH_4_mu02_CO_0_5495_r2_Photometry, IRAC_CH_4_mu02_CO_0_6_r2_Photometry, IRAC_CH_4_mu02_CO_0_69183_r2_Photometry, IRAC_CH_4_mu02_CO_0_7_r2_Photometry, IRAC_CH_4_mu02_CO_0_8_r2_Photometry, IRAC_CH_4_mu02_CO_0_8709_r2_Photometry, IRAC_CH_4_mu02_CO_0_9_r2_Photometry, IRAC_CH_4_mu02_CO_1_0_r2_Photometry, IRAC_CH_4_mu02_CO_1_0964_r2_Photometry, IRAC_CH_4_mu02_CO_1_1_r2_Photometry, IRAC_CH_4_mu02_CO_1_2_r2_Photometry, IRAC_CH_4_mu02_CO_1_3_r2_Photometry, IRAC_CH_4_mu02_CO_1_4_r2_Photometry, IRAC_CH_4_mu02_CO_1_5_r2_Photometry, IRAC_CH_4_mu02_CO_1_51356_r2_Photometry])



IRAC_mu02_Photometry = np.array([IRAC_CH_1_mu02_Photometry,IRAC_CH_2_mu02_Photometry,IRAC_CH_3_mu02_Photometry,IRAC_CH_4_mu02_Photometry]).T


IRAC_CH_1_Wavelength_Center = np.repeat(IRAC_CH_1_wavelength_center,len(X)); IRAC_CH_2_Wavelength_Center = np.repeat(IRAC_CH_2_wavelength_center,len(X)); IRAC_CH_3_Wavelength_Center = np.repeat(IRAC_CH_3_wavelength_center,len(X)); IRAC_CH_4_Wavelength_Center = np.repeat(IRAC_CH_4_wavelength_center,len(X))




####################################################################################################################################################################
############################################# Photometry in NB2315 Band (Using Definition from Charbonneau et al., 2005) ###########################################
####################################################################################################################################################################





################################################################################
############################ NB2315 Transmission ###############################
################################################################################


NB2315 = pd.read_table('NB2315.txt',sep = '\s+', dtype='unicode', header = None, index_col = None)

# Band Wavelengths (Channel Centred at 2.315 um)

NB2315_L = 22649.2 #In Angstrom
NB2315_U = 23649.2 #In Angstrom





###########################################################################
########################### NB2315 Photometry #############################
###########################################################################

Index1_NB2315 = np.where(Wavelength_Source_Interp == NB2315_L)
Index2_NB2315 = np.where(Wavelength_Source_Interp == NB2315_U)


NB2315_Wavelength_Source_Interp = Wavelength_Source_Interp[Index1_NB2315[0][0]:Index2_NB2315[0][0]]
NB2315_wavelength_center = 2.315 #In microns

NB2315_Data = np.array(NB2315)
NB2315_Data = NB2315_Data[::-1,:]




NB2315_wavelength = NB2315_Data[:,[0]]
NB2315_transmission = NB2315_Data[:,[1]] #Select the transmission column out of the two.
NB2315_Transmission = np.frompyfunc(lambda x: x.replace(',',''),1,1)(NB2315_transmission).astype(float)
NB2315_Wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(NB2315_wavelength).astype(float)*(10**(1))

NB2315_Wavelength = NB2315_Wavelength.flatten()
NB2315_Transmission = NB2315_Transmission.flatten()
NB2315_Transmission_Interp = np.interp(NB2315_Wavelength_Source_Interp, NB2315_Wavelength, NB2315_Transmission)


plt.plot(NB2315_Wavelength_Source_Interp*(10**(-4)), NB2315_Transmission_Interp,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

NB2315_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_NB2315[0][0]:Index2_NB2315[0][0]]
NB2315_Flux_Source_Interp = Flux_Source_Interp[Index1_NB2315[0][0]:Index2_NB2315[0][0]]


f2_NB2315 = (NB2315_Flux_Source_Interp)*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)



f1_NB2315_mu02_CO_0_27542_r2 = (NB2315_Planet_Flux_Table_mu02[:,0])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_27542_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_3_r2 = (NB2315_Planet_Flux_Table_mu02[:,1])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_3_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_34674_r2 = (NB2315_Planet_Flux_Table_mu02[:,2])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_34674_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_43652_r2 = (NB2315_Planet_Flux_Table_mu02[:,3])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_43652_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_4_r2 = (NB2315_Planet_Flux_Table_mu02[:,4])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_4_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_5_r2 = (NB2315_Planet_Flux_Table_mu02[:,5])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_5_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_5495_r2 = (NB2315_Planet_Flux_Table_mu02[:,6])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_5495_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_6_r2 = (NB2315_Planet_Flux_Table_mu02[:,7])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_6_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_69183_r2 = (NB2315_Planet_Flux_Table_mu02[:,8])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_69183_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_7_r2 = (NB2315_Planet_Flux_Table_mu02[:,9])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_7_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_8_r2 = (NB2315_Planet_Flux_Table_mu02[:,10])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_8_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_8709_r2 = (NB2315_Planet_Flux_Table_mu02[:,11])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_8709_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_0_9_r2 = (NB2315_Planet_Flux_Table_mu02[:,12])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_0_9_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_1_0_r2 = (NB2315_Planet_Flux_Table_mu02[:,13])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_1_0_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_1_0964_r2 = (NB2315_Planet_Flux_Table_mu02[:,14])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_1_0964_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_1_1_r2 = (NB2315_Planet_Flux_Table_mu02[:,15])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_1_1_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_1_2_r2 = (NB2315_Planet_Flux_Table_mu02[:,16])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_1_2_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_1_3_r2 = (NB2315_Planet_Flux_Table_mu02[:,17])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_1_3_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_1_4_r2 = (NB2315_Planet_Flux_Table_mu02[:,18])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_1_4_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_1_5_r2 = (NB2315_Planet_Flux_Table_mu02[:,19])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_1_5_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))


f1_NB2315_mu02_CO_1_51356_r2 = (NB2315_Planet_Flux_Table_mu02[:,20])*(NB2315_Transmission_Interp)*(NB2315_Wavelength_Source_Interp)/(h*c)
NB2315_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_NB2315_mu02_CO_1_51356_r2, NB2315_Wavelength_Source_Interp))/(simps(f2_NB2315, NB2315_Wavelength_Source_Interp))






#################################################################################
######################### NB2315 Photometric Output #############################
#################################################################################

NB2315_mu02_Photometry = np.array([[NB2315_mu02_CO_0_27542_r2_Photometry, NB2315_mu02_CO_0_3_r2_Photometry, NB2315_mu02_CO_0_34674_r2_Photometry, NB2315_mu02_CO_0_4_r2_Photometry, NB2315_mu02_CO_0_43652_r2_Photometry, NB2315_mu02_CO_0_5_r2_Photometry, NB2315_mu02_CO_0_5495_r2_Photometry, NB2315_mu02_CO_0_6_r2_Photometry, NB2315_mu02_CO_0_69183_r2_Photometry, NB2315_mu02_CO_0_7_r2_Photometry, NB2315_mu02_CO_0_8_r2_Photometry, NB2315_mu02_CO_0_8709_r2_Photometry, NB2315_mu02_CO_0_9_r2_Photometry, NB2315_mu02_CO_1_0_r2_Photometry, NB2315_mu02_CO_1_0964_r2_Photometry, NB2315_mu02_CO_1_1_r2_Photometry, NB2315_mu02_CO_1_2_r2_Photometry, NB2315_mu02_CO_1_3_r2_Photometry, NB2315_mu02_CO_1_4_r2_Photometry, NB2315_mu02_CO_1_5_r2_Photometry, NB2315_mu02_CO_1_51356_r2_Photometry]]).T



NB2315_Wavelength_Center = np.repeat(NB2315_wavelength_center,len(X))







##################################################################################################################################################################
############################################# Photometry in HST Bands (Using Definition from Charbonneau et al., 2005) ###########################################
##################################################################################################################################################################



################################################################################
######################### WFC3/HST Transmission ################################
################################################################################



# Grating-1 Wavelengths (Channel Centred at 11250.0 Angstrom)

Grating_1_L = 11000.0  #In Angstrom
Grating_1_U = 11500.0 #In Angstrom


# Grating-2 Wavelengths (Channel Centred at 11750.0 Angstrom)

Grating_2_L = 11500.0
Grating_2_U = 12000.0


# Grating-3 Wavelengths (Channel Centred at 12250.0 Angstrom)

Grating_3_L = 12000.0
Grating_3_U = 12500.0


# Grating-4 Wavelengths (Channel Centred at 12750.0 Angstrom)

Grating_4_L = 12500.0
Grating_4_U = 13000.0


# Grating-5 Wavelengths (Channel Centred at 13250.0 Angstrom)

Grating_5_L = 13000.0
Grating_5_U = 13500.0


# Grating-6 Wavelengths (Channel Centred at 13750.0 Angstrom)

Grating_6_L = 13500.0
Grating_6_U = 14000.0


# Grating-7 Wavelengths (Channel Centred at 14250.0 Angstrom)

Grating_7_L = 14000.0
Grating_7_U = 14500.0


# Grating-8 Wavelengths (Channel Centred at 14750.0 Angstrom)

Grating_8_L = 14500.0
Grating_8_U = 15000.0


# Grating-9 Wavelengths (Channel Centred at 15250.0 Angstrom)

Grating_9_L = 15000.0
Grating_9_U = 15500.0


# Grating-10 Wavelengths (Channel Centred at 15750.0 Angstrom)

Grating_10_L = 15500.0
Grating_10_U = 16000.0


# Grating-11 Wavelengths (Channel Centred at 16250.0 Angstrom)

Grating_11_L = 16000.0
Grating_11_U = 16500.0


def uniform(wl,wl_l,wl_u):
    if(wl >= wl_l and wl <= wl_u):
        return wl/wl
    else:
        return 0.0/wl



###############################################################################
######################### Grating-1 Photometry ################################
###############################################################################

Index1_Grating_1 = np.where(Wavelength_Source_Interp == Grating_1_L)
Index2_Grating_1 = np.where(Wavelength_Source_Interp == Grating_1_U)

WFC3_Grating_1_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_1[0][0]:Index2_Grating_1[0][0]]
WFC3_Grating_1_wavelength_center = 1.1250 #In microns


WFC3_Grating_1_Transmission = []

for i in range(0,len(WFC3_Grating_1_Wavelength_Source)):
    WFC3_Grating_1_Transmission.append(uniform(wl=WFC3_Grating_1_Wavelength_Source[i],wl_l=Grating_1_L,wl_u=Grating_1_U))

WFC3_Grating_1_Transmission = np.array(WFC3_Grating_1_Transmission)


plt.plot(WFC3_Grating_1_Wavelength_Source*(10**(-4)), WFC3_Grating_1_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_1_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_1[0][0]:Index2_Grating_1[0][0]]
WFC3_Grating_1_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_1[0][0]:Index2_Grating_1[0][0]]


f2_Grating_1 = (WFC3_Grating_1_Flux_Source_Interp)*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)



f1_Grating_1_mu02_CO_0_27542_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_27542_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_3_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_3_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_34674_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_34674_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_4_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_4_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_43652_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_43652_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_5_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_5_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_5495_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_5495_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_6_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_6_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_69183_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_69183_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_7_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_7_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_8_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_8_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_8709_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_8709_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_0_9_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_0_9_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_1_0_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_1_0_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_1_0964_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_1_0964_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_1_1_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_1_1_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_1_2_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_1_2_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_1_3_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_1_3_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_1_4_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_1_4_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_1_5_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_1_5_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))


f1_Grating_1_mu02_CO_1_51356_r2 = (WFC3_Grating_1_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_1_Transmission)*(WFC3_Grating_1_Wavelength_Source)/(h*c)
WFC3_Grating_1_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_1_mu02_CO_1_51356_r2, WFC3_Grating_1_Wavelength_Source))/(simps(f2_Grating_1, WFC3_Grating_1_Wavelength_Source))
















###############################################################################
######################### Grating-2 Photometry ################################
###############################################################################

Index1_Grating_2 = np.where(Wavelength_Source_Interp == Grating_2_L)
Index2_Grating_2 = np.where(Wavelength_Source_Interp == Grating_2_U)

WFC3_Grating_2_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_2[0][0]:Index2_Grating_2[0][0]]
WFC3_Grating_2_wavelength_center = 1.1750 #In microns


WFC3_Grating_2_Transmission = []

for i in range(0,len(WFC3_Grating_2_Wavelength_Source)):
    WFC3_Grating_2_Transmission.append(uniform(wl=WFC3_Grating_2_Wavelength_Source[i],wl_l=Grating_2_L,wl_u=Grating_2_U))

WFC3_Grating_2_Transmission = np.array(WFC3_Grating_2_Transmission)


plt.plot(WFC3_Grating_2_Wavelength_Source*(10**(-4)), WFC3_Grating_2_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_2_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_2[0][0]:Index2_Grating_2[0][0]]
WFC3_Grating_2_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_2[0][0]:Index2_Grating_2[0][0]]


f2_Grating_2 = (WFC3_Grating_2_Flux_Source_Interp)*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)


f1_Grating_2_mu02_CO_0_27542_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_27542_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_3_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_3_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_34674_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_34674_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_4_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_4_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_43652_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_43652_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_5_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_5_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_5495_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_5495_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_6_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_6_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_69183_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_69183_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_7_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_7_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_8_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_8_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_8709_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_8709_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_0_9_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_0_9_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_1_0_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_1_0_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_1_0964_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_1_0964_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_1_1_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_1_1_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_1_2_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_1_2_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_1_3_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_1_3_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_1_4_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_1_4_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_1_5_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_1_5_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))


f1_Grating_2_mu02_CO_1_51356_r2 = (WFC3_Grating_2_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_2_Transmission)*(WFC3_Grating_2_Wavelength_Source)/(h*c)
WFC3_Grating_2_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_2_mu02_CO_1_51356_r2, WFC3_Grating_2_Wavelength_Source))/(simps(f2_Grating_2, WFC3_Grating_2_Wavelength_Source))




###############################################################################
######################### Grating-3 Photometry ################################
###############################################################################

Index1_Grating_3 = np.where(Wavelength_Source_Interp == Grating_3_L)
Index2_Grating_3 = np.where(Wavelength_Source_Interp == Grating_3_U)

WFC3_Grating_3_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_3[0][0]:Index2_Grating_3[0][0]]
WFC3_Grating_3_wavelength_center = 1.2250 #In microns


WFC3_Grating_3_Transmission = []

for i in range(0,len(WFC3_Grating_3_Wavelength_Source)):
    WFC3_Grating_3_Transmission.append(uniform(wl=WFC3_Grating_3_Wavelength_Source[i],wl_l=Grating_3_L,wl_u=Grating_3_U))

WFC3_Grating_3_Transmission = np.array(WFC3_Grating_3_Transmission)


plt.plot(WFC3_Grating_3_Wavelength_Source*(10**(-4)), WFC3_Grating_3_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_3_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_3[0][0]:Index2_Grating_3[0][0]]
WFC3_Grating_3_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_3[0][0]:Index2_Grating_3[0][0]]


f2_Grating_3 = (WFC3_Grating_3_Flux_Source_Interp)*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)



f1_Grating_3_mu02_CO_0_27542_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_27542_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_3_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_3_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_34674_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_34674_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_4_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_4_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_43652_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_43652_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_5_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_5_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_5495_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_5495_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_6_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_6_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_69183_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_69183_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_7_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_7_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_8_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_8_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_8709_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_8709_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_0_9_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_0_9_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_1_0_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_1_0_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_1_0964_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_1_0964_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_1_1_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_1_1_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_1_2_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_1_2_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_1_3_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_1_3_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_1_4_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_1_4_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_1_5_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_1_5_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))


f1_Grating_3_mu02_CO_1_51356_r2 = (WFC3_Grating_3_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_3_Transmission)*(WFC3_Grating_3_Wavelength_Source)/(h*c)
WFC3_Grating_3_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_3_mu02_CO_1_51356_r2, WFC3_Grating_3_Wavelength_Source))/(simps(f2_Grating_3, WFC3_Grating_3_Wavelength_Source))





###############################################################################
######################### Grating-4 Photometry ################################
###############################################################################

Index1_Grating_4 = np.where(Wavelength_Source_Interp == Grating_4_L)
Index2_Grating_4 = np.where(Wavelength_Source_Interp == Grating_4_U)

WFC3_Grating_4_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_4[0][0]:Index2_Grating_4[0][0]]
WFC3_Grating_4_wavelength_center = 1.2750 #In microns


WFC3_Grating_4_Transmission = []

for i in range(0,len(WFC3_Grating_4_Wavelength_Source)):
    WFC3_Grating_4_Transmission.append(uniform(wl=WFC3_Grating_4_Wavelength_Source[i],wl_l=Grating_4_L,wl_u=Grating_4_U))

WFC3_Grating_4_Transmission = np.array(WFC3_Grating_4_Transmission)


plt.plot(WFC3_Grating_4_Wavelength_Source*(10**(-4)), WFC3_Grating_4_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_4_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_4[0][0]:Index2_Grating_4[0][0]]
WFC3_Grating_4_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_4[0][0]:Index2_Grating_4[0][0]]


f2_Grating_4 = (WFC3_Grating_4_Flux_Source_Interp)*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)





f1_Grating_4_mu02_CO_0_27542_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_27542_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_3_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_3_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_34674_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_34674_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_4_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_4_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_43652_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_43652_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_5_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_5_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_5495_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_5495_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_6_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_6_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_69183_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_69183_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_7_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_7_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_8_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_8_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_8709_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_8709_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_0_9_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_0_9_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_1_0_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_1_0_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_1_0964_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_1_0964_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_1_1_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_1_1_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_1_2_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_1_2_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_1_3_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_1_3_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_1_4_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_1_4_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_1_5_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_1_5_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))


f1_Grating_4_mu02_CO_1_51356_r2 = (WFC3_Grating_4_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_4_Transmission)*(WFC3_Grating_4_Wavelength_Source)/(h*c)
WFC3_Grating_4_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_4_mu02_CO_1_51356_r2, WFC3_Grating_4_Wavelength_Source))/(simps(f2_Grating_4, WFC3_Grating_4_Wavelength_Source))







###############################################################################
######################### Grating-5 Photometry ################################
###############################################################################

Index1_Grating_5 = np.where(Wavelength_Source_Interp == Grating_5_L)
Index2_Grating_5 = np.where(Wavelength_Source_Interp == Grating_5_U)

WFC3_Grating_5_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_5[0][0]:Index2_Grating_5[0][0]]
WFC3_Grating_5_wavelength_center = 1.3250 #In microns


WFC3_Grating_5_Transmission = []

for i in range(0,len(WFC3_Grating_5_Wavelength_Source)):
    WFC3_Grating_5_Transmission.append(uniform(wl=WFC3_Grating_5_Wavelength_Source[i],wl_l=Grating_5_L,wl_u=Grating_5_U))

WFC3_Grating_5_Transmission = np.array(WFC3_Grating_5_Transmission)


plt.plot(WFC3_Grating_5_Wavelength_Source*(10**(-4)), WFC3_Grating_5_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_5_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_5[0][0]:Index2_Grating_5[0][0]]
WFC3_Grating_5_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_5[0][0]:Index2_Grating_5[0][0]]


f2_Grating_5 = (WFC3_Grating_5_Flux_Source_Interp)*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)



f1_Grating_5_mu02_CO_0_27542_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_27542_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_3_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_3_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_34674_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_34674_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_4_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_4_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_43652_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_43652_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_5_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_5_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_5495_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_5495_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_6_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_6_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_69183_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_69183_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_7_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_7_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_8_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_8_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_8709_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_8709_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_0_9_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_0_9_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_1_0_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_1_0_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_1_0964_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_1_0964_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_1_1_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_1_1_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_1_2_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_1_2_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_1_3_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_1_3_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_1_4_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_1_4_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_1_5_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_1_5_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))


f1_Grating_5_mu02_CO_1_51356_r2 = (WFC3_Grating_5_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_5_Transmission)*(WFC3_Grating_5_Wavelength_Source)/(h*c)
WFC3_Grating_5_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_5_mu02_CO_1_51356_r2, WFC3_Grating_5_Wavelength_Source))/(simps(f2_Grating_5, WFC3_Grating_5_Wavelength_Source))





###############################################################################
######################### Grating-6 Photometry ################################
###############################################################################

Index1_Grating_6 = np.where(Wavelength_Source_Interp == Grating_6_L)
Index2_Grating_6 = np.where(Wavelength_Source_Interp == Grating_6_U)

WFC3_Grating_6_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_6[0][0]:Index2_Grating_6[0][0]]
WFC3_Grating_6_wavelength_center = 1.3750 #In microns


WFC3_Grating_6_Transmission = []

for i in range(0,len(WFC3_Grating_6_Wavelength_Source)):
    WFC3_Grating_6_Transmission.append(uniform(wl=WFC3_Grating_6_Wavelength_Source[i],wl_l=Grating_6_L,wl_u=Grating_6_U))

WFC3_Grating_6_Transmission = np.array(WFC3_Grating_6_Transmission)


plt.plot(WFC3_Grating_6_Wavelength_Source*(10**(-4)), WFC3_Grating_6_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_6_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_6[0][0]:Index2_Grating_6[0][0]]
WFC3_Grating_6_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_6[0][0]:Index2_Grating_6[0][0]]


f2_Grating_6 = (WFC3_Grating_6_Flux_Source_Interp)*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)



f1_Grating_6_mu02_CO_0_27542_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_27542_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_3_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_3_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_34674_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_34674_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_4_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_4_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_43652_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_43652_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_5_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_5_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_5495_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_5495_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_6_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_6_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_69183_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_69183_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_7_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_7_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_8_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_8_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_8709_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_8709_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_0_9_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_0_9_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_1_0_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_1_0_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_1_0964_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_1_0964_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_1_1_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_1_1_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_1_2_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_1_2_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_1_3_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_1_3_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_1_4_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_1_4_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_1_5_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_1_5_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))


f1_Grating_6_mu02_CO_1_51356_r2 = (WFC3_Grating_6_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_6_Transmission)*(WFC3_Grating_6_Wavelength_Source)/(h*c)
WFC3_Grating_6_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_6_mu02_CO_1_51356_r2, WFC3_Grating_6_Wavelength_Source))/(simps(f2_Grating_6, WFC3_Grating_6_Wavelength_Source))






###############################################################################
######################### Grating-7 Photometry ################################
###############################################################################

Index1_Grating_7 = np.where(Wavelength_Source_Interp == Grating_7_L)
Index2_Grating_7 = np.where(Wavelength_Source_Interp == Grating_7_U)

WFC3_Grating_7_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_7[0][0]:Index2_Grating_7[0][0]]
WFC3_Grating_7_wavelength_center = 1.4250 #In microns


WFC3_Grating_7_Transmission = []

for i in range(0,len(WFC3_Grating_7_Wavelength_Source)):
    WFC3_Grating_7_Transmission.append(uniform(wl=WFC3_Grating_7_Wavelength_Source[i],wl_l=Grating_7_L,wl_u=Grating_7_U))

WFC3_Grating_7_Transmission = np.array(WFC3_Grating_7_Transmission)


plt.plot(WFC3_Grating_7_Wavelength_Source*(10**(-4)), WFC3_Grating_7_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_7_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_7[0][0]:Index2_Grating_7[0][0]]
WFC3_Grating_7_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_7[0][0]:Index2_Grating_7[0][0]]


f2_Grating_7 = (WFC3_Grating_7_Flux_Source_Interp)*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)



f1_Grating_7_mu02_CO_0_27542_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_27542_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_3_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_3_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_34674_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_34674_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_4_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_4_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_43652_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_43652_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_5_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_5_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_5495_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_5495_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_6_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_6_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_69183_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_69183_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_7_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_7_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_8_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_8_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_8709_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_8709_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_0_9_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_0_9_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_1_0_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_1_0_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_1_0964_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_1_0964_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_1_1_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_1_1_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_1_2_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_1_2_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_1_3_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_1_3_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_1_4_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_1_4_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_1_5_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_1_5_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))


f1_Grating_7_mu02_CO_1_51356_r2 = (WFC3_Grating_7_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_7_Transmission)*(WFC3_Grating_7_Wavelength_Source)/(h*c)
WFC3_Grating_7_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_7_mu02_CO_1_51356_r2, WFC3_Grating_7_Wavelength_Source))/(simps(f2_Grating_7, WFC3_Grating_7_Wavelength_Source))






###############################################################################
######################### Grating-8 Photometry ################################
###############################################################################

Index1_Grating_8 = np.where(Wavelength_Source_Interp == Grating_8_L)
Index2_Grating_8 = np.where(Wavelength_Source_Interp == Grating_8_U)

WFC3_Grating_8_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_8[0][0]:Index2_Grating_8[0][0]]
WFC3_Grating_8_wavelength_center = 1.4750 #In microns


WFC3_Grating_8_Transmission = []

for i in range(0,len(WFC3_Grating_8_Wavelength_Source)):
    WFC3_Grating_8_Transmission.append(uniform(wl=WFC3_Grating_8_Wavelength_Source[i],wl_l=Grating_8_L,wl_u=Grating_8_U))

WFC3_Grating_8_Transmission = np.array(WFC3_Grating_8_Transmission)


plt.plot(WFC3_Grating_8_Wavelength_Source*(10**(-4)), WFC3_Grating_8_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_8_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_8[0][0]:Index2_Grating_8[0][0]]
WFC3_Grating_8_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_8[0][0]:Index2_Grating_8[0][0]]


f2_Grating_8 = (WFC3_Grating_8_Flux_Source_Interp)*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)



f1_Grating_8_mu02_CO_0_27542_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_27542_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_3_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_3_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_34674_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_34674_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_4_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_4_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_43652_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_43652_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_5_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_5_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_5495_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_5495_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_6_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_6_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_69183_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_69183_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_7_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_7_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_8_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_8_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_8709_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_8709_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_0_9_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_0_9_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_1_0_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_1_0_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_1_0964_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_1_0964_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_1_1_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_1_1_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_1_2_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_1_2_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_1_3_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_1_3_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_1_4_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_1_4_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_1_5_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_1_5_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))


f1_Grating_8_mu02_CO_1_51356_r2 = (WFC3_Grating_8_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_8_Transmission)*(WFC3_Grating_8_Wavelength_Source)/(h*c)
WFC3_Grating_8_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_8_mu02_CO_1_51356_r2, WFC3_Grating_8_Wavelength_Source))/(simps(f2_Grating_8, WFC3_Grating_8_Wavelength_Source))





###############################################################################
######################### Grating-9 Photometry ################################
###############################################################################

Index1_Grating_9 = np.where(Wavelength_Source_Interp == Grating_9_L)
Index2_Grating_9 = np.where(Wavelength_Source_Interp == Grating_9_U)

WFC3_Grating_9_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_9[0][0]:Index2_Grating_9[0][0]]
WFC3_Grating_9_wavelength_center = 1.5250 #In microns


WFC3_Grating_9_Transmission = []

for i in range(0,len(WFC3_Grating_9_Wavelength_Source)):
    WFC3_Grating_9_Transmission.append(uniform(wl=WFC3_Grating_9_Wavelength_Source[i],wl_l=Grating_9_L,wl_u=Grating_9_U))

WFC3_Grating_9_Transmission = np.array(WFC3_Grating_9_Transmission)


plt.plot(WFC3_Grating_9_Wavelength_Source*(10**(-4)), WFC3_Grating_9_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_9_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_9[0][0]:Index2_Grating_9[0][0]]
WFC3_Grating_9_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_9[0][0]:Index2_Grating_9[0][0]]


f2_Grating_9 = (WFC3_Grating_9_Flux_Source_Interp)*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)



f1_Grating_9_mu02_CO_0_27542_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_27542_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_3_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_3_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_34674_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_34674_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_4_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_4_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_43652_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_43652_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_5_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_5_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_5495_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_5495_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_6_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_6_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_69183_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_69183_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_7_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_7_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_8_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_8_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_8709_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_8709_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_0_9_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_0_9_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_1_0_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_1_0_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_1_0964_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_1_0964_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_1_1_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_1_1_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_1_2_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_1_2_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_1_3_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_1_3_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_1_4_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_1_4_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_1_5_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_1_5_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))


f1_Grating_9_mu02_CO_1_51356_r2 = (WFC3_Grating_9_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_9_Transmission)*(WFC3_Grating_9_Wavelength_Source)/(h*c)
WFC3_Grating_9_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_9_mu02_CO_1_51356_r2, WFC3_Grating_9_Wavelength_Source))/(simps(f2_Grating_9, WFC3_Grating_9_Wavelength_Source))






###############################################################################
######################### Grating-10 Photometry ################################
###############################################################################

Index1_Grating_10 = np.where(Wavelength_Source_Interp == Grating_10_L)
Index2_Grating_10 = np.where(Wavelength_Source_Interp == Grating_10_U)

WFC3_Grating_10_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_10[0][0]:Index2_Grating_10[0][0]]
WFC3_Grating_10_wavelength_center = 1.5750 #In microns


WFC3_Grating_10_Transmission = []

for i in range(0,len(WFC3_Grating_10_Wavelength_Source)):
    WFC3_Grating_10_Transmission.append(uniform(wl=WFC3_Grating_10_Wavelength_Source[i],wl_l=Grating_10_L,wl_u=Grating_10_U))

WFC3_Grating_10_Transmission = np.array(WFC3_Grating_10_Transmission)


plt.plot(WFC3_Grating_10_Wavelength_Source*(10**(-4)), WFC3_Grating_10_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_10_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_10[0][0]:Index2_Grating_10[0][0]]
WFC3_Grating_10_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_10[0][0]:Index2_Grating_10[0][0]]


f2_Grating_10 = (WFC3_Grating_10_Flux_Source_Interp)*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)



f1_Grating_10_mu02_CO_0_27542_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_27542_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_3_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_3_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_34674_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_34674_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_4_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_4_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_43652_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_43652_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_5_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_5_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_5495_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_5495_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_6_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_6_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_69183_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_69183_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_7_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_7_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_8_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_8_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_8709_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_8709_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_0_9_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_0_9_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_1_0_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_1_0_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_1_0964_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_1_0964_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_1_1_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_1_1_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_1_2_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_1_2_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_1_3_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_1_3_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_1_4_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_1_4_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_1_5_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_1_5_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))


f1_Grating_10_mu02_CO_1_51356_r2 = (WFC3_Grating_10_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_10_Transmission)*(WFC3_Grating_10_Wavelength_Source)/(h*c)
WFC3_Grating_10_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_10_mu02_CO_1_51356_r2, WFC3_Grating_10_Wavelength_Source))/(simps(f2_Grating_10, WFC3_Grating_10_Wavelength_Source))






################################################################################
######################### Grating-11 Photometry ################################
################################################################################

Index1_Grating_11 = np.where(Wavelength_Source_Interp == Grating_11_L)
Index2_Grating_11 = np.where(Wavelength_Source_Interp == Grating_11_U)

WFC3_Grating_11_Wavelength_Source = Wavelength_Source_Interp[Index1_Grating_11[0][0]:Index2_Grating_11[0][0]]
WFC3_Grating_11_wavelength_center = 1.6250 #In microns


WFC3_Grating_11_Transmission = []

for i in range(0,len(WFC3_Grating_11_Wavelength_Source)):
    WFC3_Grating_11_Transmission.append(uniform(wl=WFC3_Grating_11_Wavelength_Source[i],wl_l=Grating_11_L,wl_u=Grating_11_U))

WFC3_Grating_11_Transmission = np.array(WFC3_Grating_11_Transmission)


plt.plot(WFC3_Grating_11_Wavelength_Source*(10**(-4)), WFC3_Grating_11_Transmission,'g-')
plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Transmission')
plt.title('')
plt.grid(True)
plt.show()

WFC3_Grating_11_Planet_Flux_Table_mu02 = Planet_Flux_Table_mu02[Index1_Grating_11[0][0]:Index2_Grating_11[0][0]]
WFC3_Grating_11_Flux_Source_Interp = Flux_Source_Interp[Index1_Grating_11[0][0]:Index2_Grating_11[0][0]]


f2_Grating_11 = (WFC3_Grating_11_Flux_Source_Interp)*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)





f1_Grating_11_mu02_CO_0_27542_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,0])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_27542_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_27542_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_3_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,1])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_3_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_34674_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,2])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_34674_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_34674_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_4_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,3])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_4_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_43652_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,4])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_43652_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_43652_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_5_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,5])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_5_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_5495_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,6])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_5495_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_5495_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_6_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,7])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_6_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_6_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_69183_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,8])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_69183_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_69183_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_7_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,9])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_7_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_7_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_8_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,10])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_8_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_8_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_8709_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,11])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_8709_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_8709_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_0_9_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,12])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_0_9_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_0_9_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_1_0_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,13])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_1_0_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_1_0_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_1_0964_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,14])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_1_0964_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_1_0964_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_1_1_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,15])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_1_1_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_1_1_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_1_2_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,16])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_1_2_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_1_2_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_1_3_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,17])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_1_3_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_1_3_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_1_4_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,18])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_1_4_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_1_4_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_1_5_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,19])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_1_5_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_1_5_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))


f1_Grating_11_mu02_CO_1_51356_r2 = (WFC3_Grating_11_Planet_Flux_Table_mu02[:,20])*(WFC3_Grating_11_Transmission)*(WFC3_Grating_11_Wavelength_Source)/(h*c)
WFC3_Grating_11_mu02_CO_1_51356_r2_Photometry = ((Rp/Rs)**2)*(simps(f1_Grating_11_mu02_CO_1_51356_r2, WFC3_Grating_11_Wavelength_Source))/(simps(f2_Grating_11, WFC3_Grating_11_Wavelength_Source))







#################################################################################
######################### HST Photometric Output ################################
#################################################################################



WFC3_Grating_1_mu02_Photometry = np.array([WFC3_Grating_1_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_1_mu02_CO_0_3_r2_Photometry, WFC3_Grating_1_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_1_mu02_CO_0_4_r2_Photometry, WFC3_Grating_1_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_1_mu02_CO_0_5_r2_Photometry, WFC3_Grating_1_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_1_mu02_CO_0_6_r2_Photometry, WFC3_Grating_1_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_1_mu02_CO_0_7_r2_Photometry, WFC3_Grating_1_mu02_CO_0_8_r2_Photometry, WFC3_Grating_1_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_1_mu02_CO_0_9_r2_Photometry, WFC3_Grating_1_mu02_CO_1_0_r2_Photometry, WFC3_Grating_1_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_1_mu02_CO_1_1_r2_Photometry, WFC3_Grating_1_mu02_CO_1_2_r2_Photometry, WFC3_Grating_1_mu02_CO_1_3_r2_Photometry, WFC3_Grating_1_mu02_CO_1_4_r2_Photometry, WFC3_Grating_1_mu02_CO_1_5_r2_Photometry, WFC3_Grating_1_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_2_mu02_Photometry = np.array([WFC3_Grating_2_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_2_mu02_CO_0_3_r2_Photometry, WFC3_Grating_2_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_2_mu02_CO_0_4_r2_Photometry, WFC3_Grating_2_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_2_mu02_CO_0_5_r2_Photometry, WFC3_Grating_2_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_2_mu02_CO_0_6_r2_Photometry, WFC3_Grating_2_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_2_mu02_CO_0_7_r2_Photometry, WFC3_Grating_2_mu02_CO_0_8_r2_Photometry, WFC3_Grating_2_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_2_mu02_CO_0_9_r2_Photometry, WFC3_Grating_2_mu02_CO_1_0_r2_Photometry, WFC3_Grating_2_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_2_mu02_CO_1_1_r2_Photometry, WFC3_Grating_2_mu02_CO_1_2_r2_Photometry, WFC3_Grating_2_mu02_CO_1_3_r2_Photometry, WFC3_Grating_2_mu02_CO_1_4_r2_Photometry, WFC3_Grating_2_mu02_CO_1_5_r2_Photometry, WFC3_Grating_2_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_3_mu02_Photometry = np.array([WFC3_Grating_3_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_3_mu02_CO_0_3_r2_Photometry, WFC3_Grating_3_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_3_mu02_CO_0_4_r2_Photometry, WFC3_Grating_3_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_3_mu02_CO_0_5_r2_Photometry, WFC3_Grating_3_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_3_mu02_CO_0_6_r2_Photometry, WFC3_Grating_3_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_3_mu02_CO_0_7_r2_Photometry, WFC3_Grating_3_mu02_CO_0_8_r2_Photometry, WFC3_Grating_3_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_3_mu02_CO_0_9_r2_Photometry, WFC3_Grating_3_mu02_CO_1_0_r2_Photometry, WFC3_Grating_3_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_3_mu02_CO_1_1_r2_Photometry, WFC3_Grating_3_mu02_CO_1_2_r2_Photometry, WFC3_Grating_3_mu02_CO_1_3_r2_Photometry, WFC3_Grating_3_mu02_CO_1_4_r2_Photometry, WFC3_Grating_3_mu02_CO_1_5_r2_Photometry, WFC3_Grating_3_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_4_mu02_Photometry = np.array([WFC3_Grating_4_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_4_mu02_CO_0_3_r2_Photometry, WFC3_Grating_4_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_4_mu02_CO_0_4_r2_Photometry, WFC3_Grating_4_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_4_mu02_CO_0_5_r2_Photometry, WFC3_Grating_4_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_4_mu02_CO_0_6_r2_Photometry, WFC3_Grating_4_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_4_mu02_CO_0_7_r2_Photometry, WFC3_Grating_4_mu02_CO_0_8_r2_Photometry, WFC3_Grating_4_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_4_mu02_CO_0_9_r2_Photometry, WFC3_Grating_4_mu02_CO_1_0_r2_Photometry, WFC3_Grating_4_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_4_mu02_CO_1_1_r2_Photometry, WFC3_Grating_4_mu02_CO_1_2_r2_Photometry, WFC3_Grating_4_mu02_CO_1_3_r2_Photometry, WFC3_Grating_4_mu02_CO_1_4_r2_Photometry, WFC3_Grating_4_mu02_CO_1_5_r2_Photometry, WFC3_Grating_4_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_5_mu02_Photometry = np.array([WFC3_Grating_5_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_5_mu02_CO_0_3_r2_Photometry, WFC3_Grating_5_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_5_mu02_CO_0_4_r2_Photometry, WFC3_Grating_5_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_5_mu02_CO_0_5_r2_Photometry, WFC3_Grating_5_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_5_mu02_CO_0_6_r2_Photometry, WFC3_Grating_5_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_5_mu02_CO_0_7_r2_Photometry, WFC3_Grating_5_mu02_CO_0_8_r2_Photometry, WFC3_Grating_5_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_5_mu02_CO_0_9_r2_Photometry, WFC3_Grating_5_mu02_CO_1_0_r2_Photometry, WFC3_Grating_5_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_5_mu02_CO_1_1_r2_Photometry, WFC3_Grating_5_mu02_CO_1_2_r2_Photometry, WFC3_Grating_5_mu02_CO_1_3_r2_Photometry, WFC3_Grating_5_mu02_CO_1_4_r2_Photometry, WFC3_Grating_5_mu02_CO_1_5_r2_Photometry, WFC3_Grating_5_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_6_mu02_Photometry = np.array([WFC3_Grating_6_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_6_mu02_CO_0_3_r2_Photometry, WFC3_Grating_6_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_6_mu02_CO_0_4_r2_Photometry, WFC3_Grating_6_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_6_mu02_CO_0_5_r2_Photometry, WFC3_Grating_6_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_6_mu02_CO_0_6_r2_Photometry, WFC3_Grating_6_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_6_mu02_CO_0_7_r2_Photometry, WFC3_Grating_6_mu02_CO_0_8_r2_Photometry, WFC3_Grating_6_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_6_mu02_CO_0_9_r2_Photometry, WFC3_Grating_6_mu02_CO_1_0_r2_Photometry, WFC3_Grating_6_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_6_mu02_CO_1_1_r2_Photometry, WFC3_Grating_6_mu02_CO_1_2_r2_Photometry, WFC3_Grating_6_mu02_CO_1_3_r2_Photometry, WFC3_Grating_6_mu02_CO_1_4_r2_Photometry, WFC3_Grating_6_mu02_CO_1_5_r2_Photometry, WFC3_Grating_6_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_7_mu02_Photometry = np.array([WFC3_Grating_7_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_7_mu02_CO_0_3_r2_Photometry, WFC3_Grating_7_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_7_mu02_CO_0_4_r2_Photometry, WFC3_Grating_7_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_7_mu02_CO_0_5_r2_Photometry, WFC3_Grating_7_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_7_mu02_CO_0_6_r2_Photometry, WFC3_Grating_7_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_7_mu02_CO_0_7_r2_Photometry, WFC3_Grating_7_mu02_CO_0_8_r2_Photometry, WFC3_Grating_7_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_7_mu02_CO_0_9_r2_Photometry, WFC3_Grating_7_mu02_CO_1_0_r2_Photometry, WFC3_Grating_7_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_7_mu02_CO_1_1_r2_Photometry, WFC3_Grating_7_mu02_CO_1_2_r2_Photometry, WFC3_Grating_7_mu02_CO_1_3_r2_Photometry, WFC3_Grating_7_mu02_CO_1_4_r2_Photometry, WFC3_Grating_7_mu02_CO_1_5_r2_Photometry, WFC3_Grating_7_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_8_mu02_Photometry = np.array([WFC3_Grating_8_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_8_mu02_CO_0_3_r2_Photometry, WFC3_Grating_8_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_8_mu02_CO_0_4_r2_Photometry, WFC3_Grating_8_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_8_mu02_CO_0_5_r2_Photometry, WFC3_Grating_8_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_8_mu02_CO_0_6_r2_Photometry, WFC3_Grating_8_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_8_mu02_CO_0_7_r2_Photometry, WFC3_Grating_8_mu02_CO_0_8_r2_Photometry, WFC3_Grating_8_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_8_mu02_CO_0_9_r2_Photometry, WFC3_Grating_8_mu02_CO_1_0_r2_Photometry, WFC3_Grating_8_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_8_mu02_CO_1_1_r2_Photometry, WFC3_Grating_8_mu02_CO_1_2_r2_Photometry, WFC3_Grating_8_mu02_CO_1_3_r2_Photometry, WFC3_Grating_8_mu02_CO_1_4_r2_Photometry, WFC3_Grating_8_mu02_CO_1_5_r2_Photometry, WFC3_Grating_8_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_9_mu02_Photometry = np.array([WFC3_Grating_9_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_9_mu02_CO_0_3_r2_Photometry, WFC3_Grating_9_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_9_mu02_CO_0_4_r2_Photometry, WFC3_Grating_9_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_9_mu02_CO_0_5_r2_Photometry, WFC3_Grating_9_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_9_mu02_CO_0_6_r2_Photometry, WFC3_Grating_9_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_9_mu02_CO_0_7_r2_Photometry, WFC3_Grating_9_mu02_CO_0_8_r2_Photometry, WFC3_Grating_9_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_9_mu02_CO_0_9_r2_Photometry, WFC3_Grating_9_mu02_CO_1_0_r2_Photometry, WFC3_Grating_9_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_9_mu02_CO_1_1_r2_Photometry, WFC3_Grating_9_mu02_CO_1_2_r2_Photometry, WFC3_Grating_9_mu02_CO_1_3_r2_Photometry, WFC3_Grating_9_mu02_CO_1_4_r2_Photometry, WFC3_Grating_9_mu02_CO_1_5_r2_Photometry, WFC3_Grating_9_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_10_mu02_Photometry = np.array([WFC3_Grating_10_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_10_mu02_CO_0_3_r2_Photometry, WFC3_Grating_10_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_10_mu02_CO_0_4_r2_Photometry, WFC3_Grating_10_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_10_mu02_CO_0_5_r2_Photometry, WFC3_Grating_10_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_10_mu02_CO_0_6_r2_Photometry, WFC3_Grating_10_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_10_mu02_CO_0_7_r2_Photometry, WFC3_Grating_10_mu02_CO_0_8_r2_Photometry, WFC3_Grating_10_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_10_mu02_CO_0_9_r2_Photometry, WFC3_Grating_10_mu02_CO_1_0_r2_Photometry, WFC3_Grating_10_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_10_mu02_CO_1_1_r2_Photometry, WFC3_Grating_10_mu02_CO_1_2_r2_Photometry, WFC3_Grating_10_mu02_CO_1_3_r2_Photometry, WFC3_Grating_10_mu02_CO_1_4_r2_Photometry, WFC3_Grating_10_mu02_CO_1_5_r2_Photometry, WFC3_Grating_10_mu02_CO_1_51356_r2_Photometry]);
WFC3_Grating_11_mu02_Photometry = np.array([WFC3_Grating_11_mu02_CO_0_27542_r2_Photometry, WFC3_Grating_11_mu02_CO_0_3_r2_Photometry, WFC3_Grating_11_mu02_CO_0_34674_r2_Photometry, WFC3_Grating_11_mu02_CO_0_4_r2_Photometry, WFC3_Grating_11_mu02_CO_0_43652_r2_Photometry, WFC3_Grating_11_mu02_CO_0_5_r2_Photometry, WFC3_Grating_11_mu02_CO_0_5495_r2_Photometry, WFC3_Grating_11_mu02_CO_0_6_r2_Photometry, WFC3_Grating_11_mu02_CO_0_69183_r2_Photometry, WFC3_Grating_11_mu02_CO_0_7_r2_Photometry, WFC3_Grating_11_mu02_CO_0_8_r2_Photometry, WFC3_Grating_11_mu02_CO_0_8709_r2_Photometry, WFC3_Grating_11_mu02_CO_0_9_r2_Photometry, WFC3_Grating_11_mu02_CO_1_0_r2_Photometry, WFC3_Grating_11_mu02_CO_1_0964_r2_Photometry, WFC3_Grating_11_mu02_CO_1_1_r2_Photometry, WFC3_Grating_11_mu02_CO_1_2_r2_Photometry, WFC3_Grating_11_mu02_CO_1_3_r2_Photometry, WFC3_Grating_11_mu02_CO_1_4_r2_Photometry, WFC3_Grating_11_mu02_CO_1_5_r2_Photometry, WFC3_Grating_11_mu02_CO_1_51356_r2_Photometry])



WFC3_mu02_Photometry = np.array([WFC3_Grating_1_mu02_Photometry,WFC3_Grating_2_mu02_Photometry, WFC3_Grating_3_mu02_Photometry, WFC3_Grating_4_mu02_Photometry, WFC3_Grating_5_mu02_Photometry,WFC3_Grating_6_mu02_Photometry, WFC3_Grating_7_mu02_Photometry, WFC3_Grating_8_mu02_Photometry, WFC3_Grating_9_mu02_Photometry, WFC3_Grating_10_mu02_Photometry, WFC3_Grating_11_mu02_Photometry]).T


WFC3_Grating_1_Wavelength_Center = np.repeat(WFC3_Grating_1_wavelength_center,len(X));  WFC3_Grating_2_Wavelength_Center = np.repeat(WFC3_Grating_2_wavelength_center,len(X)); WFC3_Grating_3_Wavelength_Center = np.repeat(WFC3_Grating_3_wavelength_center,len(X)); WFC3_Grating_4_Wavelength_Center = np.repeat(WFC3_Grating_4_wavelength_center,len(X)); WFC3_Grating_5_Wavelength_Center = np.repeat(WFC3_Grating_5_wavelength_center,len(X)); WFC3_Grating_6_Wavelength_Center = np.repeat(WFC3_Grating_6_wavelength_center,len(X)); WFC3_Grating_7_Wavelength_Center = np.repeat(WFC3_Grating_7_wavelength_center,len(X)); WFC3_Grating_8_Wavelength_Center = np.repeat(WFC3_Grating_8_wavelength_center,len(X)); WFC3_Grating_9_Wavelength_Center = np.repeat(WFC3_Grating_9_wavelength_center,len(X)); WFC3_Grating_10_Wavelength_Center = np.repeat(WFC3_Grating_10_wavelength_center,len(X)); WFC3_Grating_11_Wavelength_Center = np.repeat(WFC3_Grating_11_wavelength_center,len(X))






###########################################################################
############################# Observed Data ###############################
###########################################################################


Croll_2011_Data = np.array(pd.read_table('croll_2011_data.txt',sep = '\s+', dtype='unicode', header = None, index_col = None), dtype=np.float64)

Crossfield_2012_Data = np.array(pd.read_table('crossfiled_2012_data.txt',sep = '\s+', dtype='unicode', header = None, index_col = None), dtype=np.float64)

Fohring_2013_Data = np.array(pd.read_table('fohring_2013_data.txt',sep = '\s+', dtype='unicode', header = None, index_col = None), dtype=np.float64)

HST_Data = np.array(pd.read_table('hst_data.txt',sep = '\s+', dtype='unicode', header = None, index_col = None), dtype=np.float64)

Lopez_Morales_2010_Data = np.array(pd.read_table('lopez-morales_2010.txt',sep = '\s+', dtype='unicode', header = None, index_col = None), dtype=np.float64)

Spitzer_Data = np.array(pd.read_table('spitzer_data.txt',sep = '\s+', dtype='unicode', header = None, index_col = None), dtype=np.float64)

Stevenson_2014_Data = np.array(pd.read_table('stevenson_2014.txt',sep = '\s+', dtype='unicode', header = None, index_col = None), dtype=np.float64)




##################################################################################################################
####################################################### Plotting #################################################
##################################################################################################################



#################################################### Box Width = 100 #############################################


pyplot.xscale('log', base=10)

color = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'gray', '#ff5733', '#668ccf', '#aaff80', '#065535', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

for i in range(0, len(X)):
    plt.plot((10**(-4))*Wavelength_Source_Interp_100_Trim[:,i], Flux_Ratio_mu02_Boxcar_100_Trim[:,i], color=color[i], lw=2, label ='$\dfrac{C}{O}$ = %s' %X[i])



plt.errorbar((10**(-4))*Spitzer_Data[:,0],(10**-2)*Spitzer_Data[:,1], (10**-2)*Spitzer_Data[:,2], color='black', fmt='o', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Spitzer/IRAC Data', zorder=3)
plt.errorbar((10**(-4))*Crossfield_2012_Data[0,:][0],(10**-2)*Crossfield_2012_Data[0,:][1], (10**-2)*Crossfield_2012_Data[0,:][2], color='black', fmt='p', ms = 9, capsize=5, capthick=1, ecolor='black', label ='Crossfield et al., 2012', zorder=3)
plt.errorbar((10**(-4))*Croll_2011_Data[:,0],(10**-2)*Croll_2011_Data[:,1], (10**-2)*Croll_2011_Data[:,2], color='black', fmt='^', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Croll et al., 2011', zorder=3)
plt.errorbar((10**(-4))*HST_Data[:,0],(10**-2)*HST_Data[:,1], (10**-2)*HST_Data[:,2], color='black', fmt='D', ms = 6, capsize=5, capthick=1, ecolor='black', label ='HST Data', zorder=3)
plt.errorbar((10**(-4))*Fohring_2013_Data[:,0],(10**-2)*Fohring_2013_Data[:,1], (10**-2)*Fohring_2013_Data[:,2], color='black', fmt='v', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Fohring et al., 2013 / Stevenson et al., 2014', zorder=3)
plt.errorbar((10**(-4))*Lopez_Morales_2010_Data[:,0],(10**-2)*Lopez_Morales_2010_Data[:,1], (10**-2)*Lopez_Morales_2010_Data[:,2], color='black', fmt='s', ms = 6, capsize=5, capthick=1, ecolor='black', label ='Lopez Morales et al., 2010', zorder=3)



class SymHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        xx=0.6*height
        return super(SymHandler, self).create_artists(legend, orig_handle, xdescent, xx, width, height, fontsize, trans)



# SPITZER

for i in range(0, len(X)):
    plt.plot(IRAC_CH_1_Wavelength_Center[i], IRAC_CH_1_mu02_Photometry[i], marker='o', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(IRAC_CH_2_Wavelength_Center, IRAC_CH_2_mu02_Photometry, marker='o', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(IRAC_CH_3_Wavelength_Center[i], IRAC_CH_3_mu02_Photometry[i], marker='o', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(IRAC_CH_4_Wavelength_Center, IRAC_CH_4_mu02_Photometry, marker='o', edgecolors='black', s=50, color=color, zorder=3)


# HST

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_1_Wavelength_Center[i], WFC3_Grating_1_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_2_Wavelength_Center, WFC3_Grating_2_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_3_Wavelength_Center[i], WFC3_Grating_3_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_4_Wavelength_Center, WFC3_Grating_4_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_5_Wavelength_Center[i], WFC3_Grating_5_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_6_Wavelength_Center, WFC3_Grating_6_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_7_Wavelength_Center[i], WFC3_Grating_7_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_8_Wavelength_Center, WFC3_Grating_8_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_9_Wavelength_Center[i], WFC3_Grating_9_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_10_Wavelength_Center, WFC3_Grating_10_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_11_Wavelength_Center[i], WFC3_Grating_11_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])


# NB2315

for i in range(0, len(X)):
    plt.plot(NB2315_Wavelength_Center[i], NB2315_mu02_Photometry[i], marker='p', ms = 7, markeredgecolor='black', mfc=color[i])



plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio')
plt.title('Planet-to-Star Flux Ratio of Wasp-12 System with mu 02 and different $\dfrac{C}{O}$, Boxcar Width = 100')
plt.grid(True)
leg = plt.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, fontsize='medium', ncol=2, handleheight=1.3, labelspacing=0.7)
plt.show()





#################################################### Box Width = 250 ##############################################

pyplot.xscale('log', base=10)


color = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'gray', '#ff5733', '#668ccf', '#aaff80', '#065535', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

for i in range(0, len(X)):
    plt.plot((10**(-4))*Wavelength_Source_Interp_250_Trim[:,i], Flux_Ratio_mu02_Boxcar_250_Trim[:,i], color=color[i], lw=2, label ='$\dfrac{C}{O}$ = %s' %X[i])



plt.errorbar((10**(-4))*Spitzer_Data[:,0],(10**-2)*Spitzer_Data[:,1], (10**-2)*Spitzer_Data[:,2], color='black', fmt='o', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Spitzer/IRAC Data', zorder=3)
plt.errorbar((10**(-4))*Crossfield_2012_Data[0,:][0],(10**-2)*Crossfield_2012_Data[0,:][1], (10**-2)*Crossfield_2012_Data[0,:][2], color='black', fmt='p', ms = 9, capsize=5, capthick=1, ecolor='black', label ='Crossfield et al., 2012', zorder=3)
plt.errorbar((10**(-4))*Croll_2011_Data[:,0],(10**-2)*Croll_2011_Data[:,1], (10**-2)*Croll_2011_Data[:,2], color='black', fmt='^', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Croll et al., 2011', zorder=3)
plt.errorbar((10**(-4))*HST_Data[:,0],(10**-2)*HST_Data[:,1], (10**-2)*HST_Data[:,2], color='black', fmt='D', ms = 6, capsize=5, capthick=1, ecolor='black', label ='HST Data', zorder=3)
plt.errorbar((10**(-4))*Fohring_2013_Data[:,0],(10**-2)*Fohring_2013_Data[:,1], (10**-2)*Fohring_2013_Data[:,2], color='black', fmt='v', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Fohring et al., 2013 / Stevenson et al., 2014', zorder=3)
plt.errorbar((10**(-4))*Lopez_Morales_2010_Data[:,0],(10**-2)*Lopez_Morales_2010_Data[:,1], (10**-2)*Lopez_Morales_2010_Data[:,2], color='black', fmt='s', ms = 6, capsize=5, capthick=1, ecolor='black', label ='Lopez Morales et al., 2010', zorder=3)



class SymHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        xx=0.6*height
        return super(SymHandler, self).create_artists(legend, orig_handle, xdescent, xx, width, height, fontsize, trans)


# SPITZER

for i in range(0, len(X)):
    plt.plot(IRAC_CH_1_Wavelength_Center[i], IRAC_CH_1_mu02_Photometry[i], marker='o', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(IRAC_CH_2_Wavelength_Center, IRAC_CH_2_mu02_Photometry, marker='o', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(IRAC_CH_3_Wavelength_Center[i], IRAC_CH_3_mu02_Photometry[i], marker='o', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(IRAC_CH_4_Wavelength_Center, IRAC_CH_4_mu02_Photometry, marker='o', edgecolors='black', s=50, color=color, zorder=3)


# HST

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_1_Wavelength_Center[i], WFC3_Grating_1_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_2_Wavelength_Center, WFC3_Grating_2_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_3_Wavelength_Center[i], WFC3_Grating_3_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_4_Wavelength_Center, WFC3_Grating_4_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_5_Wavelength_Center[i], WFC3_Grating_5_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_6_Wavelength_Center, WFC3_Grating_6_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_7_Wavelength_Center[i], WFC3_Grating_7_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_8_Wavelength_Center, WFC3_Grating_8_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_9_Wavelength_Center[i], WFC3_Grating_9_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_10_Wavelength_Center, WFC3_Grating_10_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_11_Wavelength_Center[i], WFC3_Grating_11_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])


# NB2315

for i in range(0, len(X)):
    plt.plot(NB2315_Wavelength_Center[i], NB2315_mu02_Photometry[i], marker='p', ms = 7, markeredgecolor='black', mfc=color[i])




plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio')
plt.title('Planet-to-Star Flux Ratio of Wasp-12 System with mu 02 and different $\dfrac{C}{O}$, Boxcar Width = 250')
plt.grid(True)
leg = plt.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, fontsize='medium', ncol=2, handleheight=1.3, labelspacing=0.7)
plt.show()




#################################################### Box Width = 500 ##############################################

pyplot.xscale('log', base=10)


color = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'gray', '#ff5733', '#668ccf', '#aaff80', '#065535', '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

for i in range(0, len(X)):
    plt.plot((10**(-4))*Wavelength_Source_Interp_500_Trim[:,i], Flux_Ratio_mu02_Boxcar_500_Trim[:,i], color=color[i], lw=2, label ='$\dfrac{C}{O}$ = %s' %X[i])



plt.errorbar((10**(-4))*Spitzer_Data[:,0],(10**-2)*Spitzer_Data[:,1], (10**-2)*Spitzer_Data[:,2], color='black', fmt='o', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Spitzer/IRAC Data', zorder=3)
plt.errorbar((10**(-4))*Crossfield_2012_Data[0,:][0],(10**-2)*Crossfield_2012_Data[0,:][1], (10**-2)*Crossfield_2012_Data[0,:][2], color='black', fmt='p', ms = 9, capsize=5, capthick=1, ecolor='black', label ='Crossfield et al., 2012', zorder=3)
plt.errorbar((10**(-4))*Croll_2011_Data[:,0],(10**-2)*Croll_2011_Data[:,1], (10**-2)*Croll_2011_Data[:,2], color='black', fmt='^', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Croll et al., 2011', zorder=3)
plt.errorbar((10**(-4))*HST_Data[:,0],(10**-2)*HST_Data[:,1], (10**-2)*HST_Data[:,2], color='black', fmt='D', ms = 6, capsize=5, capthick=1, ecolor='black', label ='HST Data', zorder=3)
plt.errorbar((10**(-4))*Fohring_2013_Data[:,0],(10**-2)*Fohring_2013_Data[:,1], (10**-2)*Fohring_2013_Data[:,2], color='black', fmt='v', ms = 8, capsize=5, capthick=1, ecolor='black', label ='Fohring et al., 2013 / Stevenson et al., 2014', zorder=3)
plt.errorbar((10**(-4))*Lopez_Morales_2010_Data[:,0],(10**-2)*Lopez_Morales_2010_Data[:,1], (10**-2)*Lopez_Morales_2010_Data[:,2], color='black', fmt='s', ms = 6, capsize=5, capthick=1, ecolor='black', label ='Lopez Morales et al., 2010', zorder=3)



class SymHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        xx=0.6*height
        return super(SymHandler, self).create_artists(legend, orig_handle, xdescent, xx, width, height, fontsize, trans)

# SPITZER

for i in range(0, len(X)):
    plt.plot(IRAC_CH_1_Wavelength_Center[i], IRAC_CH_1_mu02_Photometry[i], marker='o', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(IRAC_CH_2_Wavelength_Center, IRAC_CH_2_mu02_Photometry, marker='o', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(IRAC_CH_3_Wavelength_Center[i], IRAC_CH_3_mu02_Photometry[i], marker='o', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(IRAC_CH_4_Wavelength_Center, IRAC_CH_4_mu02_Photometry, marker='o', edgecolors='black', s=50, color=color, zorder=3)


# HST

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_1_Wavelength_Center[i], WFC3_Grating_1_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_2_Wavelength_Center, WFC3_Grating_2_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_3_Wavelength_Center[i], WFC3_Grating_3_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_4_Wavelength_Center, WFC3_Grating_4_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_5_Wavelength_Center[i], WFC3_Grating_5_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_6_Wavelength_Center, WFC3_Grating_6_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_7_Wavelength_Center[i], WFC3_Grating_7_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_8_Wavelength_Center, WFC3_Grating_8_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_9_Wavelength_Center[i], WFC3_Grating_9_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])

matplotlib.pyplot.scatter(WFC3_Grating_10_Wavelength_Center, WFC3_Grating_10_mu02_Photometry, marker='D', edgecolors='black', s=50, color=color, zorder=3)

for i in range(0, len(X)):
    plt.plot(WFC3_Grating_11_Wavelength_Center[i], WFC3_Grating_11_mu02_Photometry[i], marker='D', ms = 7, markeredgecolor='black', mfc=color[i])


# NB2315

for i in range(0, len(X)):
    plt.plot(NB2315_Wavelength_Center[i], NB2315_mu02_Photometry[i], marker='p', ms = 7, markeredgecolor='black', mfc=color[i])




plt.xlabel('Wavelength ($\mu m$)')
plt.ylabel('Flux Ratio')
plt.title('Planet-to-Star Flux Ratio of Wasp-12 System with mu 02 and different $\dfrac{C}{O}$, Boxcar Width = 500')
plt.grid(True)
leg = plt.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, fontsize='medium', ncol=2, handleheight=1.3, labelspacing=0.7)
plt.show()






############################################################################################################################
####################################################### 𝛘² Goodness of Fit #################################################
############################################################################################################################


# HST Photometry

CS_HST = []
for i in range(0, len(X)):
    cs_hst = np.sum((((HST_Data[:,1]*(10**(-2)))-WFC3_mu02_Photometry[i,:])**2)/((HST_Data[:,2]*(10**(-2)))**2))
    CS_HST.append(cs_hst)

Chi_Square_HST = np.stack((X, CS_HST), axis=-1)

Reduced_Chi_Square_HST = np.stack((X, [r / (HST_Data.shape[0]-0) for r in CS_HST]), axis=-1)


# Spitzer Photometry

CS_Spitzer = []
for i in range(0, len(X)):
    cs_spitzer = np.sum((((Spitzer_Data[:,1]*(10**(-2)))-IRAC_mu02_Photometry[i,:])**2)/((Spitzer_Data[:,2]*(10**(-2)))**2))
    CS_Spitzer.append(cs_spitzer)

Chi_Square_Spitzer = np.stack((X, CS_Spitzer), axis=-1)

Reduced_Chi_Square_Spitzer = np.stack((X, [r / (Spitzer_Data.shape[0]-0) for r in CS_Spitzer]), axis=-1)



# NB2315 Photometry

CS_NB2315 = []
for i in range(0, len(X)):
    cs_nb2315 = np.sum((((Crossfield_2012_Data[:,1]*(10**(-2)))-NB2315_mu02_Photometry[i,:])**2)/((Crossfield_2012_Data[:,2]*(10**(-2)))**2))
    CS_NB2315.append(cs_nb2315)

Chi_Square_NB2315 = np.stack((X, CS_NB2315), axis=-1)

Reduced_Chi_NB2315 = np.stack((X, [r / (Crossfield_2012_Data.shape[0]-0) for r in CS_NB2315]), axis=-1)



# Combined Photometry

Chi_Square_Combined = np.stack((X, Chi_Square_HST[:,1]+Chi_Square_Spitzer[:,1]+Chi_Square_NB2315[:,1]), axis=-1)
Reduced_Chi_Square_Combined = np.stack((X, Chi_Square_Combined[:,1]/(HST_Data.shape[0]+Spitzer_Data.shape[0]+Crossfield_2012_Data.shape[0]-0)), axis=-1)





###################################################################################################################################
####################################################### Write 𝛘² values in a file #################################################
###################################################################################################################################


np.savetxt('Chi_Square_mu02.dat', Chi_Square_Combined, delimiter='	', header='', comments='')

np.savetxt('Reduced_Chi_Square_mu02.dat', Reduced_Chi_Square_Combined, delimiter='	', header='', comments='')




##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################




