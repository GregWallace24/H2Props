# H2Props

File: "'H2PropTableGenerator.py"

Tabulates properties of hydrogen (normal, para, ortho) above the triple point

The tables are put into a file, with the name depending on whether it is gas/liquid/supercrical, para/ortho/normal, and the temp and pressure ranges.  Tables for all 3 phases will always be generated, but some might be empty depending on the temp and pressure range used (as seen in the example tables).

Caution: this code is intended to output properties, even when they aren't very accurately known.  Viscosity is only known for normal hydrogen, and so para and ortho hydrogen output normal viscosity.  Thermal conductivity is known for para above 50K, but not below.  No conductivity data is known for orthohydrogen.  So for para, below 50K, it adopts normal hydrogen conductivity.  For orhto, it assumes mass-averaged between para and normal.  You will need to decide for yourself whether the numbers look accurate enough for your applications.  

I'm not much of a coder, so much of the code is written very linearly, with the intent that a non-coder (most engineers) can understand with only the basic syntax known.



File: "Oequilb.py"

Prints out the equilibrium orthohydrogen fraction of hydrogen, based off the input temperature (in Kelvin).



File: "OrthoFractionPlot.py"

Creates a table of equilibrium orthohydrogen fraction in bulk hydrogen as a function of temperature (in Kelvin).



File: "OrthoFractionPlotDiffSpins.py"

Compares the output of "OrthoFractionPlot.py" for two different values of hydrogen's characteristic rotational temperature.  The one used is for "low temperatures", but a second value is given for "high temperatures", as described in the source of the formula.



File: "RelEnthPlots.py"

Tabulates an plots the relative enthalpy between normal hydrogen, orthohydrogen, parahydrogen, and equilibrium hydrogen for used in energy balances requiring hydrogen conversion
~703 kJ/kg is released when converting from 100% orthohydrogen to 100% parahydrogen at the normal boiling point.



File: "RelEnthDiffPlots1.py" and "RelEnthDiffPlots2.py"

These tabulate and plot the difference in relative enthalpy of normal hydrogen, parahydrogen, orthohydrogen, and equilibrium hydrogen.  The first uses parahydrogen as the base (it is the lowest energy).  The second uses equilibrium hydrogen.  



File: "AvailableCooling.py"

Calculates the amount of cooling power available if hydrogen is converted from primarily parahydrogen to having more orthohydrogen.  Takes input temperature and pressure.  

Assumes initial orthohydrogen fraction is equal to the liquid-temperature equilibrium (based off the input pressure).  Assumes relative enthalpy is constant during conversion (heat is absorbed at constant temperature to the new equilibrium, at the input temperature).  


