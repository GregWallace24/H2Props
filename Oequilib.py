#Author: Greg Wallace (11326742)
#Date Started: July 3, 2020
#Date Last Edited: July 3, 2020
#Organization: Washington State Univeristy
#Purpose: Calcualte the fraction of hydrogen in the orthohydrogen spin state based off temperature

#import CoolProp.CoolProp as CP #not actually needed for this code
import math #For exp() function



T=77#In Kelvin, this is the input



#Takes temperature between 0 and 500K (higher is probably okay) and outputs the ortho fraction (between 0 at low temps and 0.75 at high temps)
n=7#iterator: how many times it shouold iterate.  Number shoul be odd.  Number should be above 5, but higher is more accurate
if 0 == (n%2):
        print("CAUTION: Iterator should be an odd number for best results")
        #n++

K_para=0.0#these need will be summations
K_ortho=0.0

for J in range(0,n+1):
    if 0 == (J%2):
        #print("even")
        K_para=K_para+(2*J+1)*math.exp(-J*(J+1)*85.4/T)
    elif 1 == (J%2):
        #print("odd")
        K_ortho=K_ortho+(2*J+1)*math.exp(-J*(J+1)*85.4/T)
    else:
        print("ERROR: iterator must be whole number")

X=3*K_ortho/(3*K_ortho+K_para)#Fraction of ortho as opposed to all (ortho + para)

print("Ortho Fraction: ")
print(X)

#At 20K, it should be about 0.00 (0%)
#At 77K, it shoul dbe about 0.49 (49%)
#At 300K and above, it should be 0.75 (75%)

#Source of Formula: Brandt Pedrow's thesis: https://wsuwp-uploads.s3.amazonaws.com/uploads/sites/44/2016/05/b_pedrow_011423571.pdf (pdf page 21, or numbered page 8)
#Citation: Brandt Patrick Pedrow, May 2016, "Parahydrogen-Orthohydrogen Conversion on Catalyst Loaded Scrim for Vapor Cooled Shielding of Cryogenic Storage Vessels", Master Thesis, Washington State University, Pullman.
