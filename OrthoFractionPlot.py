#Title: OrthoFractionPlot.py
#Author: Greg Wallace
#Organization: Washington State University, HYPER Lab
#Date Created: July 29, 2021
#Date Last Edited: September 6, 2021
#Purpose: Plot the equilibrium orthohydrogen fraction as a function of temperature for my thesis


import math #for "floor" function
import matplotlib.pyplot as plt
import numpy as np


#This function returns equilibrium ortho fraction based off temperature and number of iterations to use.  If no iterations are given, defaults to 7, which computs accurately from 0 to 300K up to 16 decimals
#At 20K, it should be about 0.00 (0%)
#At 77K, it shoul dbe about 0.49 (49%)
#At 300K and above, it should be 0.75 (75%)
#Source of Formula: Brandt Pedrow's thesis: https://wsuwp-uploads.s3.amazonaws.com/uploads/sites/44/2016/05/b_pedrow_011423571.pdf (pdf page 21, or numbered page 8)
#Citation: Brandt Patrick Pedrow, May 2016, "Parahydrogen-Orthohydrogen Conversion on Catalyst Loaded Scrim for Vapor Cooled Shielding of Cryogenic Storage Vessels", Master Thesis, Washington State University, Pullman.
def Yo_equilib(T,N=7,T_rot=85.4):
    if N<1 or (type(N) != int):#returns if N is not an int of 1 or above
        print("Yo_equilib requires integer N value of 1 or above")
        return None
    if T<=0:#returns if T is negative or 0, because negative temps don't exist
        print("Yo_equilib requires temperature must be in positive Kelvin")
        return None
    J=0
    K_para=0.0
    K_ortho=0.0
    #T_rot=85.4#characteristic rotational temperature of hydrogen, 85.4K
    for i in range(1,N+1):
        #print(i)
        K_para=K_para+(2*J+1)*math.exp(-J*(J+1)*T_rot/T)
        J+=1
        K_ortho=K_ortho+(2*J+1)*math.exp(-J*(J+1)*T_rot/T)
        J+=1
    Yo=3*K_ortho/(K_para+3*K_ortho)
    return Yo


T=[]
Yo=[]
#Generate the data points:
for i in range(15,300+1):
    T.append(i)
    Yo.append(Yo_equilib(i))


#Set font properties before generating plot
plt.rcParams.update({'font.family':"Times New Roman"})
plt.rcParams.update({'font.size':12})


fig, YoPlot = plt.subplots()#Create teh plot

YoPlot.plot(T,Yo)#Plot the data we want
YoPlot.set_xlabel('Temperature [K]')
YoPlot.set_ylabel('Equilibrium Orthohydrogen Mass Fraction [-]')
YoPlot.set_xticks(np.arange(0, 300+20, 20))
YoPlot.set_yticks(np.arange(0, .75+.05, .05))
YoPlot.set_xlim([0,300])
YoPlot.set_ylim([0,0.751])
YoPlot.grid(True)

plt.show()

#https://matplotlib.org/stable/gallery/lines_bars_and_markers/psd_demo.html#sphx-glr-gallery-lines-bars-and-markers-psd-demo-py


