#Title: "H2_Functions" 
#Author: Greg Wallace (11326742)
#Organization: Washington State Univeristy
#Last Edited: 2022.27.04
#Purpose: A bunch of useful functions I've written for hydrogen properties, liquefiers, ortho-para stuff, etc.  



#Libraries:
import CoolProp.CoolProp as CP #where we get all the thermofluid data
import math #for "floor" function
import matplotlib.pyplot as plt
import numpy as np


#Functions:

#Returns the equilibrium fraction of orthohydrogen (between 0 and 1) based off input temperature, and optionally the number of iterations to use and an optional rotational temperature)
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

#A couple that return relative enthalpy at the saturation point, based off pressure/temperature and ortho fraction
def h_satL_mixP(P,Yo):
    h_ortho=CP.PropsSI('H','Q',0.0,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','Q',0.0,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return (Yo*h_ortho+(1-Yo)*h_para)

def h_satG_mixP(P,Yo):
    h_ortho=CP.PropsSI('H','Q',1.0,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','Q',1.0,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return (Yo*h_ortho+(1-Yo)*h_para)

def h_satL_mixT(T,Yo):
    h_ortho=CP.PropsSI('H','Q',0.0,'T',T,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','Q',0.0,'T',T,'parahydrogen') #Enthalpy in (J/kg)
    return (Yo*h_ortho+(1-Yo)*h_para)

def h_satG_mixT(T,Yo):
    h_ortho=CP.PropsSI('H','Q',1.0,'T',T,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','Q',1.0,'T',T,'parahydrogen') #Enthalpy in (J/kg)
    return (Yo*h_ortho+(1-Yo)*h_para)

#Returns the Orthohydrogen fraction of ortho-para mixture based off a relative enthalpy
def Yo_mix(T,P,h_mix):
    h_ortho=CP.PropsSI('H','T',T,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','T',T,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return ((h_mix-h_para)/(h_ortho-h_para))

def Yo_SatL_mixP(P,h_mix):
    h_ortho=CP.PropsSI('H','Q',0.0,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','Q',0.0,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return ((h_mix-h_para)/(h_ortho-h_para))

def Yo_SatG_mixP(P,h_mix):
    h_ortho=CP.PropsSI('H','Q',1.0,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','Q',1.0,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return ((h_mix-h_para)/(h_ortho-h_para))

def Yo_SatL_mixT(T,h_mix):
    h_ortho=CP.PropsSI('H','T',T,'Q',0.0,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','T',T,'Q',0.0,'parahydrogen') #Enthalpy in (J/kg)
    return ((h_mix-h_para)/(h_ortho-h_para))

def Yo_SatG_mixT(T,h_mix):
    h_ortho=CP.PropsSI('H','T',T,'Q',1.0,'orthohydrogen') #Enthalpy in (J/kg)
    h_para=CP.PropsSI('H','T',T,'Q',1.0,'parahydrogen') #Enthalpy in (J/kg)
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

#Returns the Relative Internal Energy of ortho-para mixture based off mass-weighted average
def u_mix(T,P,Yo):
    u_ortho=CP.PropsSI('U','T',T,'P',P,'orthohydrogen') #Internal Energy in (J/kg)
    u_para=CP.PropsSI('U','T',T,'P',P,'parahydrogen') #Internal Energy in (J/kg)
    return (Yo*u_ortho+(1-Yo)*u_para)

#Returns the Orthohydrogen fraction of ortho-para mixture based off a relative enthalpy
def Yo_mixu(T,P,u_mix):
    u_ortho=CP.PropsSI('U','T',T,'P',P,'orthohydrogen') #Enthalpy in (J/kg)
    u_para=CP.PropsSI('U','T',T,'P',P,'parahydrogen') #Enthalpy in (J/kg)
    return ((u_mix-u_para)/(u_ortho-u_para))

#Finds final temperature when hydrogen is catalyzed from a known T_initial, P, and h_mix to an equilibrium T and Yo.  Isenthalpic, isobaric process.
#Temperature between 0K and 500K
def T_isenthu(P,u_mixture):
    i=0
    T_guess=150.0
    T_max=500.0
    T_min=14.0
    u_guess=1000.0
    u_err=1.0
    Yo_guess=0.375
    while abs(h_err) > 0.000001:#Ends when enthalpy error is less than 0.001J/kg
        Yo_guess=Yo_equilib(T_guess)
        u_guess=u_mix(T_guess,P,Yo_guess)
        print("Temp Guess: " + str(T_guess))
        print("Enth Guess: " + str(u_guess))
        u_err=u_guess-u_mixture
        print("Int E err: " + str(u_err))
        if u_err>0:
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

def s_mix_rough(T,P,Yo):
    s_ortho=CP.PropsSI('S','T',T,'P',P,'orthohydrogen') #Entropy in (J/kg-K)
    s_para=CP.PropsSI('S','T',T,'P',P,'parahydrogen') #Entropy in (J/kg-K)
    return (Yo*s_ortho+(1-Yo)*s_para)



#Setting Up Coolprop:
P_ref=101325#Pa
#Parahydrogen Reference State does not need to change:
h_ref=CP.PropsSI('H','P',P_ref,'Q',0,'parahydrogen')#should be about zero
s_ref=CP.PropsSI('S','P',P_ref,'Q',0,'parahydrogen')#should be roughly zero
#print(h_ref)#Should be about 0
#print(s_ref)#Should be about 0
#Orthohydrogen Reference State needs to be adjusted such that its enthalpy is 702.98 kJ/kg and entropy is 0.018269 kJ/kg-K at NBP
Dmolar_ref_ortho=CP.PropsSI('Dmolar','P',P_ref,'Q',0.0,'orthohydrogen')#reference state must be in molar density at liquid of normal boiling point
T_ref_ortho=CP.PropsSI('T','P',P_ref,'Q',0.0,'orthohydrogen')
CP.set_reference_state('orthohydrogen',T_ref_ortho,Dmolar_ref_ortho,702.98*(2.0*1.00784),0.018269*(2.0*1.00784)) #sets so that enthalpy and entropy difference correctly correspond to ~702.98kJ/kg and 0.018269 kJ/kg-K at 20K; but function inputs must be molar
#Gotten from Jacob Leachman's blog post: https://hydrogen.wsu.edu/2015/06/22/why-equilibrium-hydrogen-doesnt-exist/
#print(CP.PropsSI('H','P',P_ref,'Q',0,'orthohydrogen'))#should be 702980 J/kg
#print(CP.PropsSI('S','P',P_ref,'Q',0,'orthohydrogen'))#should be 18.269 J/kg-K

#Previous method, don't use anymore:
# CP.set_reference_state('orthohydrogen',20.3800689304,35150.6373702,1417.12332,0.036828) 
# fluid_thermo ='orthohydrogen'



#Main Code: 

Q=10#W, cooling power in watts of the cryocooler
P=14.7*6894.76#psi * 6894.76 = Pa
print("Pressure (Pa): "+str(P))
Ti=CP.PropsSI('T','Q',0.0,'P',101325,'nitrogen')
print("Intitial Temp [K]: "+str(Ti))

Yo_i=Yo_equilib(Ti)
print("Initial Ortho Fraction [-]: "+str(Yo_i))
h_initial=h_mix(Ti,P,Yo_i)
print("Initial Enthalpy [J/kg]: "+str(h_initial))


Tf=CP.PropsSI('T','Q',0.0,'P',P,'hydrogen')
print("Final Temp [K]: "+str(Tf))
Yo_f=Yo_equilib(Tf)
print("Final Ortho Fraction [-]: "+str(Yo_f))

h_final_cat=h_satL_mixP(P,Yo_f)
print("Final h [J/kg], cat: "+str(h_final_cat))
h_final_nocat=h_satL_mixP(P,Yo_i)
print("Final h [J/kg], nocat: "+str(h_final_nocat))

dh_cat=(h_initial-h_final_cat)#J/kg
dh_nocat=(h_initial-h_final_nocat)#J/kg

m_dot_cat=1000*Q/dh_cat#g/s, the multiple is converting kg to g
m_dot_nocat=1000*Q/dh_nocat#g/s
print("Liquefaction rate [g/s] with "+str(Q)+" W cooling, with catalysis: "+str(m_dot_cat))
print("Liquefaction rate [g/s] with "+str(Q)+" W cooling, without catalysis: "+str(m_dot_nocat))

