#Title: RelEnthDiffPlots2.py
#Author: Greg Wallace
#Organization: Washington State University, HYPER Lab
#Date Created: September 6, 2021
#Date Last Edited: September 6, 2021
#Purpose: Plot the difference in relative enthalpies of orthohydrogen, normal hydrogen, and parahydrogen compared to equilibrium.


import CoolProp.CoolProp as CP #Source of all thermodynamic data
import math #for "floor" function
import matplotlib.pyplot as plt #For plotting
import numpy as np #For creating plot axes


#Returns equilibrium ortho fraction based off temperature and number of iterations to use.  If no iterations are given, defaults to 7, which computs accurately from 0 to 300K up to 16 decimals
#At 20K, it should be about 0.00 (0%)
#At 77K, it shoul dbe about 0.49 (49%)
#At 300K and above, it should be 0.75 (75%)
#Source of Formula: Brandt Pedrow's thesis: https://wsuwp-uploads.s3.amazonaws.com/uploads/sites/44/2016/05/b_pedrow_011423571.pdf (pdf page 21, or numbered page 8)
#Citation: Brandt Patrick Pedrow, May 2016, "Parahydrogen-Orthohydrogen Conversion on Catalyst Loaded Scrim for Vapor Cooled Shielding of Cryogenic Storage Vessels", Master Thesis, Washington State University, Pullman.
def Yo_equilib(T,N=7):
    if N<1 or (type(N) != int):#returns if N is not an int of 1 or above
        print("Yo_equilib requires integer N value of 1 or above")
        return None
    if T<=0:#returns if T is negative or 0, because negative temps don't exist
        print("Yo_equilib requires temperature must be in positive Kelvin")
        return None
    J=0
    K_para=0.0
    K_ortho=0.0
    T_rot=85.4#characteristic rotational temperature of hydrogen, 85.4K
    for i in range(1,N+1):
        #print(i)
        K_para=K_para+(2*J+1)*math.exp(-J*(J+1)*T_rot/T)
        J+=1
        K_ortho=K_ortho+(2*J+1)*math.exp(-J*(J+1)*T_rot/T)
        J+=1
    Yo=3*K_ortho/(K_para+3*K_ortho)
    return Yo

#Returns the Relative Enthalpy of ortho-para mixture based off mass-weighted average
def h_mix(T,P,Yo):
    h_ortho=CP.PropsSI('H','T',T,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','T',T,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return (Yo*h_ortho+(1-Yo)*h_para)


#Setting Up Coolprop:
P_ref=101325#Pa
#Parahydrogen Reference State does not need to change:
h_ref=CP.PropsSI('H','P',P_ref,'Q',0,'parahydrogen')#should be about zero
s_ref=CP.PropsSI('S','P',P_ref,'Q',0,'parahydrogen')#should be roughly zero
#print(h_ref)#Should be about 0
#print(s_ref)#Should be about 0
#Orthohydrogen Reference State needs to be adjusted such that its enthalpy is 702.98 kJ/kg and entropy is 0.018269 kJ/kg-K at NBP
Dmolar_ref_ortho=CP.PropsSI('Dmolar','P',P_ref,'Q',0,'orthohydrogen')#reference state must be in molar density at liquid of normal boiling point
T_ref_ortho=CP.PropsSI('T','P',P_ref,'Q',0.5,'orthohydrogen')
CP.set_reference_state('orthohydrogen',T_ref_ortho,Dmolar_ref_ortho,702.98*(2*1.00784),0.018269*(2*1.00784)) #sets so that enthalpy and entropy difference correctly correspond to ~702.98kJ/kg and 0.018269 kJ/kg-K at 20K; but function inputs must be molar
#print(CP.PropsSI('H','P',P_ref,'Q',0,'orthohydrogen'))#should be 702980 J/kg
#print(CP.PropsSI('S','P',P_ref,'Q',0,'orthohydrogen'))#should be 18.269 J/kg-K


P=101325#Pa, pressure to test at
T=[]
Yo=[]
ReldH75=[]
ReldH100=[]
ReldH00=[]
for i in range(15,300+1):#Properties of hydrogen do not go below 14K, Room temp is about 300K
    T.append(i)
    Yo.append(Yo_equilib(i))
    ReldH100.append(h_mix(i,P,1.0)/1000-h_mix(i,P,Yo_equilib(i))/1000)#convert to kJ/kg, orthohydrogen - equilib
    ReldH75.append(h_mix(i,P,0.75)/1000-h_mix(i,P,Yo_equilib(i))/1000)#convert to kJ/kg, normal hydrogen - equilib
    ReldH00.append(h_mix(i,P,0.0)/1000-h_mix(i,P,Yo_equilib(i))/1000)#convert to kJ/kg, normal parahydrogen - equilib




#Set font properties before generating plot
plt.rcParams.update({'font.family':"Times New Roman"})
plt.rcParams.update({'font.size':12})

fig, ReldHPlot=plt.subplots()#start plot
ReldHPlot.plot(T,ReldH100,color='red')#plot enthalpy difference between orthohydrogen and equilibrium
ReldHPlot.plot(T,ReldH75,color='green')#plot enthalpy difference between normal hydrogen and equilibrium
ReldHPlot.plot(T,ReldH00,color='darkblue')#plot enthalpy difference between parahydrogen and equilibrium
ReldHPlot.set_xlabel('Temperature [K]')
ReldHPlot.set_ylabel('Relative Enthalpy Difference from Equilibrium [kJ/kg]')
ReldHPlot.set_xticks(np.arange(0, 300+20, 20))
ReldHPlot.set_yticks(np.arange(-400, 750+50, 50))
ReldHPlot.set_xlim([0,300])
ReldHPlot.set_ylim([-400,751])
ReldHPlot.grid(True)
ReldHPlot.legend(['Orthohydrogen','Normal Hydrogen','Parahydrogen'])
ReldHPlot.axes.axhline(y=0, color='black', linewidth=1)#Make horizontal line at y=0
plt.show()

#https://matplotlib.org/stable/gallery/lines_bars_and_markers/psd_demo.html#sphx-glr-gallery-lines-bars-and-markers-psd-demo-py
#https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html#matplotlib.axes.Axes.legend
#https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
#https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.axhline.html