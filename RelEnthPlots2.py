#Title: RelEnthPlots2.py
#Author: Greg Wallace
#Organization: Washington State University, HYPER Lab
#Date Created: December 3, 2021
#Date Last Edited: December 6, 2021
#Purpose: Plots (and saves a printable paper-sided pdf of) the relative enthalpies of parahydrogen, orthohydrogen, normal hydrogen, and equilibrium hydrogen using ortho increments.


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

# def h_ortho(T,P):
    # return CP.PropsSI('H','T',T,'P',P,'orthohydrogen') #Enthalpy in (J/kg)

#Returns the Relative Enthalpy of ortho-para mixture based off mass-weighted average
def h_mix(T,P,Yo):
    h_ortho=CP.PropsSI('H','T',T,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','T',T,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return (Yo*h_ortho+(1-Yo)*h_para)

#Returns the Orthohydrogen fraction of ortho-para mixture based off a relative enthalpy
def Yo_mix(T,P,h_mix):
    h_ortho=CP.PropsSI('H','T',T,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','T',T,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return ((h_mix-h_para)/(h_ortho-h_para))

#Finds final temperature when hydrogen is catalyzed from a known T_initial, P, and h_mix to an equilibrium T and Yo.  Isenthalpic, isobaric process.
#Temperature between 0K and 500K
def T_isenth(P,h_mixture):
    i=0
    T_guess=150.0
    T_max=500.0
    T_min=14.0
    h_guess=1000.0
    h_err=1.0
    Yo_guess=0.375
    while abs(h_err) > 0.000001:#Ends when enthalpy error is less than 0.001J/kg
        Yo_guess=Yo_equilib(T_guess)
        h_guess=h_mix(T_guess,P,Yo_guess)
        print("Temp Guess: " + str(T_guess))
        print("Enth Guess: " + str(h_guess))
        h_err=h_guess-h_mixture
        print("Enth err: " + str(h_err))
        if h_err>0:
            T_max=T_guess
            T_guess=(T_max+T_min)/2
        else:
            T_min=T_guess
            T_guess=(T_max+T_min)/2
        
        if i>1000:
            print("Did not converge in 1000 iterations!")
            return None
        i+=1
    print(i)
    return T_guess

#Returns the reverse carnot refrigeration efficiency when pulling heat out at T_L and pumping it up to T_H
def COPRefrig(T_H,T_L):
    if T_L>T_H:print("Warning: T_L is higher than T_H")
    return (T_L/(T_H-T_L))
#https://en.wikipedia.org/wiki/Heat_pump_and_refrigeration_cycle

#Finds the temperature where a certain equilibrium fraction occurs.  
#Temperature reasonably accurate between 14K and 400K, but less accurate at very high and very low temps because ortho fraction approaches a constant (0 or 0.75)
def T_equilib(Yo):
    if(Yo>=0.75 or Yo<=0.0):
        print("Orthohydrogen fraction must be above 0 and below 0.75, it does not exist as equilibrium outside this range.")
        return None
    i=0
    
    T_max=1000.0
    T_guess=T_max
    T_min=0
    Yo_err=0.25
    Yo_guess=0.375
    while abs(Yo_err) > 0.000001:#Ends when equilibrium error is less than 1/1000 of 0.1%
        if Yo_err>0:
            T_max=T_guess
            T_guess=(T_max+T_min)/2
        else:
            T_min=T_guess
            T_guess=(T_max+T_min)/2
        Yo_guess=Yo_equilib(T_guess)

        print("Temp Guess: " + str(T_guess))
        print("Yo Guess: " + str(Yo_guess))
        Yo_err=Yo_guess-Yo
        print("Yo err: " + str(Yo_err))
        
        
        if i>1000:
            print("Did not converge in 1000 iterations!")
            return None
        i+=1
    print(str(i) + " iterations")
    return T_guess




#Setting Up Coolprop:
P_ref=101325#Pa
#Parahydrogen Reference State does not need to change:
h_ref=CP.PropsSI('H','P',P_ref,'Q',0,'parahydrogen')#should be about zero
s_ref=CP.PropsSI('S','P',P_ref,'Q',0,'parahydrogen')#should be roughly zero
#print(h_ref)
#print(s_ref)
#Orthohydrogen Reference State needs to be adjusted such that its enthalpy is 702.98 kJ/kg and entropy is 0.018269 kJ/kg-K at NBP
Dmolar_ref_ortho=CP.PropsSI('Dmolar','P',P_ref,'Q',0,'orthohydrogen')#reference state must be in molar density at liquid of normal boiling point
T_ref_ortho=CP.PropsSI('T','P',P_ref,'Q',0.5,'orthohydrogen')
CP.set_reference_state('orthohydrogen',T_ref_ortho,Dmolar_ref_ortho,702.98*(2*1.00784),0.018269*(2*1.00784)) #sets so that enthalpy and entropy difference correctly correspond to ~702.98kJ/kg and 0.018269 kJ/kg-K at 20K; but function inputs must be molar
#print(CP.PropsSI('H','P',P_ref,'Q',0,'orthohydrogen'))#should be 702980 J/kg
#print(CP.PropsSI('S','P',P_ref,'Q',0,'orthohydrogen'))#should be 18.269 J/kg-K


P=101325#Pa, pressure to test at

#Temperature range to take enthalpy values at
T_low=15#K
T_high=150#K
dT=0.1#K, step size
N_T=round((T_high-T_low)/dT)+1#Number of temperatures to check
T_arr=np.linspace(T_low,T_high,N_T)#makes an array of temperatures to check between the T_low and T_high
#print(T_arr)

Yo_low=0.0#unitless
Yo_high=0.75
dYo=0.05#Step size
N_Yo=round((Yo_high-Yo_low)/dYo)+1#Number of ortho fractions to arrange
Yo_arr=np.linspace(Yo_low,Yo_high,N_Yo)#makes an array of ortho fractions to check between the Yo_low and Yo_high
#print(Yo_arr)

h_arrs=np.zeros([N_Yo,N_T])#initialize array of array of enthalpies.  Primary array is for that ortho fraction, secondary is for each temp: [Yo1[T1,T2,T3...],Yo2[T1,T2,T3...]]
h_eq=np.zeros(N_T)#initialize array for the equilibrium enthalpies.  Each value in the array is for a specific temp

Yo_Names=[]

for k in range(N_T):#For each temperature
    for i in range(N_Yo):#for each of the ortho-fraciton setpoints (ex: 0%, 5%, 10%, etc.)
        #Generate enthalpy at each temp and ortho fraction
        h_arrs[i][k]=h_mix(T_arr[k],P,Yo_arr[i])/1000#kJ/kg, do for each temp and each ortho fraction
        if 0==k:
            Yo_Names.append(str(round(Yo_arr[i]*100))+" %")#Create string name for that ortho fraction
    #Generate enthalpy at equilib at that temp
    h_eq[k]=h_mix(T_arr[k],P,Yo_equilib(T_arr[k]))/1000#kJ/kg, equilibrium,only do for each temp
#print(h_arrs)
#print(h_eq)

Yo_Names.append("Equilibrium")

#Flip the names and the enthalpy arrays (only flips the order of the ortho-fractions, the temperature relations will not change
#This is used for plotting (so that it will plot highest ortho fracitons first, and lowest last), to make the plot's key show the lines in order of greatest to smallest
Yo_Names=np.flip(Yo_Names)
h_arrs=np.flip(h_arrs,0)

#Now the Plot:
#Set font properties before generating plot
plt.rcParams.update({'font.family':"Times New Roman",'font.size':10})

fig, h_plot=plt.subplots(figsize=(7.5,10))#start plot, give size to be 7.5" wide, 10" tall
h_plot.plot(T_arr,h_eq,color='black',linewidth=1.0)#Prints equilib curve line because I want it to show up on top in the Key
for j in range(N_Yo):
    h_plot.plot(T_arr,h_arrs[j],linewidth=1.0)
h_plot.plot(T_arr,h_eq,color='black',linewidth=1.0)#Re-prints the equilib curve a second time because I want it to show up on top of the other curves
h_plot.set_xlabel('Temperature [K]')
h_plot.set_ylabel('Relative Enthalpy [kJ/kg]')
h_plot.set_xticks(np.arange(0, 300+5, 5))
h_plot.set_yticks(np.arange(-500, 2500+50, 50))
h_plot.set_xlim([15,150])
h_plot.set_ylim([-100,2501])
h_plot.grid(True)
h_plot.legend(Yo_Names,title='Ortho Fraction:',fontsize=10,loc="lower right")
plt.savefig("RelEnthPlots2.pdf",format="pdf")#Save the figure as a pdf (vector) which can be printed vertically
#File is saved in whatever the working directory is when you run the file (ex: the user folder, or the local directory, etc.)
plt.show()#Show the plot

#https://matplotlib.org/stable/gallery/lines_bars_and_markers/psd_demo.html#sphx-glr-gallery-lines-bars-and-markers-psd-demo-py
#https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html#matplotlib.axes.Axes.legend
#https://stackoverflow.com/questions/3899980/how-to-change-the-font-size-on-a-matplotlib-plot
#https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure.savefig
#https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linewidth
#https://matplotlib.org/stable/gallery/subplots_axes_and_figures/figure_size_units.html
#https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linestyle