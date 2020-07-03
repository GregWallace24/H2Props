"""
Author: Greg Wallace
Company: Washington State University
Last Edited: July 2020

Adapted from: 
Carl Bunge
Washington State University
June 2018

Adapted by Carl from @author: Luka Denies from TU Delft.

"""

import CoolProp.CoolProp as CP #Source of all thermodynamic data
import math #for "floor" function



#****************************************************************************************

#Set a reference state so that enthalpy and entropy are related between ortho/para/normal hydrogen
#Not necessary if you don't care if ortho-para conversion energy matches

P_ref=101325#Pa

#For hydrogen, enthalpy and entropy are taken to be 0 in the liquid phase at the normal boiling point
#This is fine, unless you want the conversion between ortho and para hydrogen already included
#By setting the reference state of othohydrogen at it's normal boiling point (liquid) to be H=702.89kJ/kg, 
#S=0.018269kJ/kg-K, then the difference in enthalpy and entropy will be taken into account automatically.  
#This way, at an orthohydrogen fraction other than 0 or 1, enthalpy and entropy can taken as a mass average between
#the two states at the same temperature and pressure.  
#See: https://hydrogen.wsu.edu/2015/06/22/why-equilibrium-hydrogen-doesnt-exist/
h_ref=CP.PropsSI('H','P',P_ref,'Q',0,'parahydrogen')#should be about zero
s_ref=CP.PropsSI('S','P',P_ref,'Q',0,'parahydrogen')#should be roughly zero
print(h_ref)
print(s_ref)

Dmolar_ref_ortho=CP.PropsSI('Dmolar','P',P_ref,'Q',0,'orthohydrogen')#reference state must be in molar density at liquid of normal boiling point
T_ref_ortho=CP.PropsSI('T','P',P_ref,'Q',0.5,'orthohydrogen')
CP.set_reference_state('orthohydrogen',T_ref_ortho,Dmolar_ref_ortho,702.98*(2*1.00784),0.018269*(2*1.00784)) #sets so that enthalpy and entropy difference correctly correspond to ~702.98kJ/kg and 0.018269 kJ/kg-Kat 20K; but function inputs must be molar
print(CP.PropsSI('H','P',P_ref,'Q',0,'orthohydrogen'))#should be 702980 J/kg
print(CP.PropsSI('S','P',P_ref,'Q',0,'orthohydrogen'))#should be 18.269 J/kg-K

#We can do the same for normal hydrogen, which is 75% ortho and 25% para, but it is unlikely you will ever use this
Dmolar_ref_normal=CP.PropsSI('Dmolar','P',P_ref,'Q',0,'hydrogen')#reference state must be in molar density at liquid of normal boiling point
T_ref_normal=CP.PropsSI('T','P',P_ref,'Q',0.5,'hydrogen')
CP.set_reference_state('hydrogen',T_ref_normal,Dmolar_ref_normal,702.98*(2*1.00784)*0.75,0.018269*(2*1.00784)*0.75) #sets so that enthalpy and entropy difference correctly correspond to ~702.98kJ/kg and 0.018269 kJ/kg-Kat 20K; but function inputs must be molar
print(CP.PropsSI('H','P',P_ref,'Q',0,'hydrogen'))#should be 702980 J/kg
print(CP.PropsSI('S','P',P_ref,'Q',0,'hydrogen'))#should be 18.269 J/kg-K



#****************************************************************************************
#Choose the type of hydrogen: parahydrogen, orthohydrogen, or normal hydrogen
fluid_thermo ='orthohydrogen'
#fluid_thermo = 'parahydrogen'
#fluid_thermo ='hydrogen'


#Fluid for transport model (viscosity will not always be defined, otherwise)
fluid_transport = 'hydrogen'


#Property (our symbol, coolprop symbol): Units
#Temperature (T,T): K
#Pressure (P,P): Pa
#Enthalpy (H,H): J/kg
#Density (rho,D): kg/m^3
#Specific Heat (Cp,C): J/kgK (constant pressure)
#Entropy (S,S): J/kgK
#Speed of Sound (c,A): m/s
#Internal Energy (E,U): J/kg
#Dynamic Viscosity (mu,V): Pa*s
#Thermal Conductivity (kappa,L): J/kgK




#****************************************************************************************

#Temperature limits

#Start, End, and Iterator of first temp range
T0=100 #K
T1=180 #K
dT1=1 #K

#Start, End, and Iterator of second temp range, if desired
# T2=math.ceil((CP.PropsSI(fluid_thermo,'Tcrit'))/10)*10 #K, automatically chooses the first multiple of 10 above the critical temp
# T3=90 #K
# dT2=10 #K

#List of all Temps to read at:
TRange = []
for i in range(T0,T1+1,dT1):#This adds all the first temperatrues starting at T0 and going until one iteration before the second number (T1+1 means it will include on the T1 number)
    TRange.append(i)
# for i in range(T2,T3+1,dT2):#The second temperature range, if used
    # TRange.append(i)
print(TRange)



#****************************************************************************************

#Pressure limits

P0 = 5e4 #Pa
P1 = 50e4 #Pa
dP=1e4 #Pa

#List all the Pressures:
pRange = []
# while i<=P1:
    # pRange.append(i)
    # i=i+dP
# print(pRange)
for i in range(int(P0),int(P1+dP),int(dP)):#Can't do the same for pressure as done for temp because by using "e" as a power of 10, it treats the number as a float and not an integer.  We could use the for-loop if we wrote out all the numbers: 100000,500000,50000.
    pRange.append(i)
    i=i+dP
print(pRange)   


#Saturation Temperature Limits:
# Tsat0=T0
# Tsat1=math.floor(CP.PropsSI(fluid_thermo,'Tcrit'))
# dTsat=1

# TsatRange=[]
# for i in range(Tsat0,Tsat1+1,dTsat):
    # TsatRange.append(i)
# print(TsatRange)



#****************************************************************************************
#Initiallize Arrays: 

#Critical Point (to determine phase)
T_crit=CP.PropsSI(fluid_thermo,'Tcrit')
P_crit=CP.PropsSI(fluid_thermo,'Pcrit')

props = []
props.append(["Temperature","Pressure","Density","Speed","Conductivity","Enthalpy","Entropy","Viscosity","InternalE","Cp","Cv","CpMCv"])
#print(props)

rho = []    #Density (kg/m^3)
mu = []     #Dynamic viscosity (Pa*s)
kappa = []  #Thermal conductivity (J/kgK)
Cp = []     #Heat capacity (constant pressure) (J/kgK)
H = []      #Enthalpy (J/kg)
CpMCv = []  #Specific Gas Constant (R_sp) (J/kgK)
E = []      #Internal Energy (J/kg)
S = []      #Entropy (J/kgK)
c = []      #Speed of Sound (m/s)


#Saturation Properties at L and G, and saturation pressure:
# pSat = []

# satL_rho = []
# satL_mu = []
# satL_kappa = []
# satL_Cp = []
# satL_H = []
# satL_CpMCv = []
# satL_E = []
# satL_S = []
# satL_c = []

# satG_rho = []
# satG_mu = []
# satG_kappa = []
# satG_Cp = []
# satG_H = []
# satG_CpMCv = []
# satG_E = []
# satG_S = []
# satG_c = []

#also do sat properties in another set of arrays so we don't need another file?




#****************************************************************************************

#obtain properties (populate property arrays)

for p in pRange:#Iterates through all pressures
    for T in TRange:#iterates through all temepratures (for each pressure)
        #get all the data and put into it's array... Do everything in this second level of the loop
        
        #Order: ["Temperature","Pressure","Density","Speed","Conductivity","Enthalpy","Entropy","Viscosity","InternalE","Cp","Cv","CpMCv"]
        array=[None]*len(props[0])#initializes empty array of the same length as the header row
        
        array[0]=T #adds temperature to the array in K
        array[1]=p #adds pressure to the array in Pa
        array[2] = CP.PropsSI('D','T',T,'P',p,fluid_thermo)#adds density in kg/m^3
        rhoCur = CP.PropsSI('D','T',T,'P',p,fluid_thermo)#current density, used for further properties because it is more robust than T,p if close to liquid/gas mixture.  Probably not necessary
        array[3] = CP.PropsSI('A','D',rhoCur,'P',p,fluid_thermo) #speed of sound (m/s)
        #array[4] = CP.PropsSI('L','D',rhoCur,'P',p,fluid_thermo)
        #Thermal conductivity is difficult 
        #It isn't formulated for ortho, so I take a mass average of normal and para, which probably isn't accurate at all
        #Also, for para, it throws an error below 49.407K because the formulas are inaccurate.  So for that, we will assume normal hydrogen conductivity, which probably isn't very accurate
        #the CoolProp devs did this intentionally, see explanation here: https://github.com/CoolProp/CoolProp/blob/master/FAQ.md
        if "orthohydrogen" == fluid_thermo:
            #array[4] = ((CP.PropsSI('L','D',rhoCur,'T',T,fluid_transport)-(0.25*CP.PropsSI('L','D', rhoCur,'T',T,'REFPROP::parahydrogen')))*(1.3333333))
            if 50<=T:
                array[4] = (1/0.75)*CP.PropsSI('L','D',rhoCur,'T',T,'hydrogen')-(1/0.75)*(0.25/1)*CP.PropsSI('L','D',rhoCur,'T',T,'parahydrogen')
            else:
                array[4] = CP.PropsSI('L','D',rhoCur,'T',T,'hydrogen')
            #"Ortho"
            #No thermal conductitity models are available for orthohydrogen, 
        elif "parahydrogen" == fluid_thermo:
            #array[4] = CP.PropsSI('L','D',rhoCur,'T',T,'REFPROP::parahydrogen')
            if 50<=T:
                array[4] = CP.PropsSI('L','D',rhoCur,'T',T,fluid_thermo)
            else:
                array[4] = CP.PropsSI('L','D',rhoCur,'T',T,'hydrogen')
            #"Para"
        elif "hydrogen" ==fluid_thermo:
            array[4] = CP.PropsSI('L','D',rhoCur,'T',T,fluid_thermo)
            #"Normal"
        #Thermal conductivitgy (W/m-K)
        #otherwise, thermal conductivity is left empty
        array[5] = CP.PropsSI('H','D',rhoCur,'P',p,fluid_thermo) #Enthalpy in (J/kg)
        array[6] = CP.PropsSI('S','D',rhoCur,'P',p,fluid_thermo) #Entropy in (J/kgK)
        array[7] = CP.PropsSI('V','D',rhoCur,'P',p,fluid_transport) #Dynamic Viscosity (Pa-s)
        #No viscosity models are available for parahydrogen or orthohydrogen, but it should be exactly the same as normal hydrogen
        array[8] = CP.PropsSI('U','D',rhoCur,'P',p,fluid_thermo) #Internal Energy in J/kg
        array[9] = CP.PropsSI('C','D',rhoCur,'P',p,fluid_thermo) #Heat Capacity (Cp) in J/kgK (const pressure)
        array[10] = CP.PropsSI('O','D',rhoCur,'P',p,fluid_thermo) #Heat Capacity (Cv) in J/kgK (const volume)
        array[11] = array[9]-array[10]#Cp-Cv=R, specific gas constnat (J/kgK)
        props.append(array)

#print(props)

print("ranges complete")



#****************************************************************************************

#write to file

#don't really need to make a function here, but it is a residual of a previous version of this code and it works just fine
def writetofile(list):#this is a function definition.  It has to be defined BEFORE the function is ever called.  
    
    #decides what fluid name to give the files
    fluid=fluid_thermo
    if 'orthohydrogen' == fluid_thermo:
        fluid="ortho"
    elif 'parahydrogen' == fluid_thermo:
        fluid="para"
    elif 'hydrogen' == fluid_thermo:
        fluid="normal"
    print(fluid)
    
    #Create files for liquid, gas, and supercritical, which include the temperature and pressure ranges in their filenames
    propFileSuper=open("S" + fluid + "_" + str(T0) + "-" + str(T1) + "K" + "_" + str(math.floor(P0)) + "-" + str(math.floor(P1)) + "Pa" + ".csv","w")
    propFileLiquid=open("L" + fluid + "_" + str(T0) + "-" + str(T1) + "K" + "_" + str(math.floor(P0)) + "-" + str(math.floor(P1)) + "Pa" + ".csv","w")
    propFileGas=open("G" + fluid + "_" + str(T0) + "-" + str(T1) + "K" + "_" + str(math.floor(P0)) + "-" + str(math.floor(P1)) + "Pa" + ".csv","w")
    
    
    for i in list:
        if i == list[0]: #if we are looking at the row, then write in the headers for all of the files
            for n in i:
                propFileSuper.write(n+",")
                propFileLiquid.write(n+",")
                propFileGas.write(n+",")
            propFileSuper.write("\n")
            propFileLiquid.write("\n")
            propFileGas.write("\n")
            
        else: #add the row's data to it's proper file by properties
            if i[0] >= T_crit and i[1] >= P_crit: #if above the critical point, then a supercritical fluid, so add it's values and then newline
                for n in i:
                    propFileSuper.write(str(n)+",")
                propFileSuper.write("\n")
            elif i[0] >= T_crit: #if not supercritical but above the critical temperature, then a gas
                for n in i:
                    propFileGas.write(str(n)+",")
                propFileGas.write("\n")
            elif i[1] >= P_crit: #if not supercritical but above critical pressure, then a liquid
                for n in i:
                    propFileLiquid.write(str(n)+",")
                propFileLiquid.write("\n")
            elif i[0] > CP.PropsSI('T','P',i[1],'Q',0.5,fluid_thermo): #checks if temp is greater than the saturation temp at that pressure, if so, then it is a gas
                for n in i:
                    propFileGas.write(str(n)+",")
                propFileGas.write("\n")
            elif i[0] < CP.PropsSI('T','P',i[1],'Q',0.5,fluid_thermo): #checks if temp is less than the saturation temp at that pressure, if so, then it is a liquid
                for n in i:
                    propFileLiquid.write(str(n)+",")
                propFileLiquid.write("\n")
            else:
                print("Value Ignored")
            #Ignores values that are in liquid-gas-super phase    
        
        #Note: Star CCM+ does not care if there is an extra comma before the newline key, so this should still upload just fine with the extra commas

    propFileSuper.close()
    propFileLiquid.close()
    propFileGas.close()
    
writetofile(props)

