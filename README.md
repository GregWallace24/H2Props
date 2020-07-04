# H2Props
Tabulates properties of hydrogen (normal, para, ortho) above the triple point

The tables are put into a file, with the name depending on whether it is gas/liquid/supercrical, para/ortho/normal, and the temp and pressure ranges.  Tables for all 3 phases will always be generated, but some might be empty depending on the temp and pressure range used (as seen in the example tables).

Caution: this code is intended to output properties, even when they aren't very accurately known.  Viscosity is only known for normal hydrogen, and so para and ortho hydrogen output normal viscosity.  Thermal conductivity is known for para above 50K, but not below.  No conductivity data is known for orthohydrogen.  So for para, below 50K, it adopts normal hydrogen conductivity.  For orhto, it assumes mass-averaged between para and normal.  You will need to decide for yourself whether the numbers look accurate enough for your applications.  

I'm not much of a coder, so much of the code is written very linearly, with the intent that a non-coder (most engineers) can understand with only the basic syntax known.
