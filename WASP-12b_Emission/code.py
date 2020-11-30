#!/usr/bin/env python3
import sys
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import emcee
import corner
np.set_printoptions(threshold=sys.maxsize)
plt.rcParams['figure.figsize'] = [12,7]
#Data = pd.read_table('try.dat',sep = '\s+', dtype='unicode', header = None, index_col = None)
#print(Data)

--------------------------------------------------------------------------

#Use the code below
# Initialize a list for the data to be stored
data = []

# Iterate through your file to read the data
with open("try.dat") as f:
    for line in f.readlines():
        # Use .rstrip() to get rid of the newline character at the end
        data.append(line.rstrip("\r\n"))
# Assumes that data is the result from the above code
data = [i[:14] + " " + i[14:] if len(i) > 14 else i for i in data]

---------------------------------------------------------------------

myfile = open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w')

with open("lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7") as data:
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
with open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7', 'w') as file:
  file.write(filedata)



Data = pd.read_table('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu02.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu02 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu02 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
error_mu02 = Spectra_Data[:,[16]]
error_mu02 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu02).astype(float)
Error_mu02 = error_mu02[wavelength1[0][0]:wavelength2[0][0],[0]]

plt.plot(Wavelength_mu02,Flux_mu02,'r-')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()




-------------------------------------------------------------------------------------------
myfile = open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu03.spec.7', 'w')

with open("lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu03.spec.7") as data:
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
with open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu03.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu03.spec.7', 'w') as file:
  file.write(filedata)


Data = pd.read_table('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu03.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu03 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu03 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
error_mu03 = Spectra_Data[:,[16]]
error_mu03 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu03).astype(float)
Error_mu03 = error_mu03[wavelength1[0][0]:wavelength2[0][0],[0]]

plt.plot(Wavelength_mu03,Flux_mu03,'g-')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()


----------------------------------------------------------------------------------------------------

myfile = open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu04.spec.7', 'w')

with open("lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu04.spec.7") as data:
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
with open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu04.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu04.spec.7', 'w') as file:
  file.write(filedata)
  
  
Data = pd.read_table('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu04.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu04 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu04 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
error_mu04 = Spectra_Data[:,[16]]
error_mu04 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu04).astype(float)
Error_mu04 = error_mu04[wavelength1[0][0]:wavelength2[0][0],[0]]

plt.plot(Wavelength_mu04,Flux_mu04,'b-')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()


------------------------------------------------------------------------------------------------------

myfile = open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu05.spec.7', 'w')

with open("lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu05.spec.7") as data:
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
with open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu05.spec.7', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('D', 'e')

# Write the file out again
with open('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu05.spec.7', 'w') as file:
  file.write(filedata)

Data = pd.read_table('Modified_lte001-2.99-0.0a+0.0.BT-Settl.CIFIST2011_2017.D-1.0e8-irr6300_0.022-r2-mu05.spec.7',sep = '\s+', dtype='unicode', header = None, index_col = None)
Spectra_Data = np.array(Data)


wavelength = Spectra_Data[:,[0]]
flux = Spectra_Data[:,[1]]
flux = np.frompyfunc(lambda x: x.replace(',',''),1,1)(flux).astype(float)
wavelength = np.frompyfunc(lambda x: x.replace(',',''),1,1)(wavelength).astype(float)
wavelength1 = np.where(wavelength == 1000.)
wavelength2 = np.where(wavelength == 100000.)
Wavelength_mu05 = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux_mu05 = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
error_mu05 = Spectra_Data[:,[16]]
error_mu05 = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error_mu05).astype(float)
Error_mu05 = error_mu05[wavelength1[0][0]:wavelength2[0][0],[0]]

plt.plot(Wavelength_mu05,Flux_mu05,'k-')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()



------------------------------------------------------------------------------------------------------

myfile = open('Modified_lte063-4.5-0.0a+0.0.BT-Settl.spec.7_1', 'w')

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
with open("Modified_lte063-4.5-0.0a+0.0.BT-Settl.spec.7_1") as data:
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
Wavelength = wavelength[wavelength1[0][0]:wavelength2[0][0],[0]]
Flux = flux[wavelength1[0][0]:wavelength2[0][0],[0]]
error = Spectra_Data[:,[16]]
error = np.frompyfunc(lambda x: x.replace(',',''),1,1)(error).astype(float)
Error = error[wavelength1[0][0]:wavelength2[0][0],[0]]


plt.plot(Wavelength,Flux,'y-')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()














--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Rs = 1.630*(696340*((10)**3))
Rp = 1.790*(69911*((10)**3))
Ts = 6300
h = 6.62607*((10)**-34)
c = 299792458
k = 1.3806*((10)**-23)
Flux_mu05_Flatten = Flux_mu05.flatten()
Wavelength_Interp = Wavelength_mu05.flatten()
Wavelength = Wavelength.flatten()
Flux = Flux.flatten()
Flux_Interp = np.interp(Wavelength_Interp,Wavelength, Flux)

Error = Error.flatten()
Error_Interp = np.interp(Wavelength_Interp,Wavelength, Error)
Error_mu05_Flatten = Error_mu05.flatten()


#############################################################################
##############################Plotting#######################################
#############################################################################


plt.plot(Wavelength_Interp,Flux_Interp,'c-')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Flux (FLAM)') 
plt.title('') 
plt.grid(True) 
plt.show()

Flux_Ratio = Flux_mu05_Flatten/Flux_Interp
#axes = plt.gca()
#axes.set_xlim([1000,100000])
plt.plot(Wavelength_Interp,Flux_Ratio,color = '#5c8281')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Flux Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()


Brightness_Ratio = Flux_Ratio*((Rs/Rp)**2)

plt.plot(Wavelength_Interp,Brightness_Ratio,color = '#a5b557')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Brightness Ratio') 
plt.title('') 
plt.grid(True) 
plt.show()

Planet_Planck_Function = ((((np.exp((h*c)/(Wavelength_Interp*(10**(-10))*k*Ts)))-1)**(-1))*Brightness_Ratio*((2*h*c**2)/(Wavelength_Interp**5)))

#Fitting_Function = ((2*h*c**2)/(Wavelength_Interp**5))*(((np.exp((h*c)/(Wavelength_Interp*(10**(-10))*k*3000)))-1)**(-1))*(10**(7))
plt.plot(Wavelength_Interp,Planet_Planck_Function,color = '#a5b557')
#plt.plot(Wavelength_Interp,Fitting_Function,color = '#8b63f7')
plt.xlabel('Wavelength ($\AA$)') 
plt.ylabel('Planet Brightness') 
plt.title('') 
plt.grid(True) 
plt.show()



# Error Propogation

#Combined_Error = Flux_Ratio*(((((Error_mu05_Flatten)/(Flux_mu05_Flatten))**2)+(((Error_Interp)/(Flux_Interp))**2))**0.5)



























-----------------------------------------------------------------------------------------------
filters = ['HST/WFC3_UVIS1.F300X', 'HST/WFC3_UVIS1.F606W', 'HST/WFC3_UVIS1.F845M', 'HST/WFC3_IR.F098M', 'HST/WFC3_IR.F160W']
























