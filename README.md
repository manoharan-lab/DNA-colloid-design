# DNA-colloid-design

The code for generating random sequences is in randSeqGenBest.m
This code take in (len,seq,max,num), which correspond to:
len: the length of the sequence that you want to generate
seq: the sequence of other DNA strands that you do not want to bind to your generated sequence
max: the number of repeated 3 base regions that you will allow in your new sequence
num: the number of times you want the code to try to generate a sequence for you

The code for generating melting curves for particles coated in complementary DNA is in twoSequenceMeltingCurvePlot.m
The code for generating melting curves for particles coated in complementary DNA with a single displacing strand is in threeSeqReentrantMeltingCurvePlot.m
The code for generating melting curves for particles coated in complementary DNA with two displacing strands is in fourSeqReentrantMeltingCurvePlot.m
These packages take in some subset of (seqA, seqD1, seqD2, Cd1, Cd2, saltConc), which correspond to:
seqA: the DNA sequence that will hold your two particles together when bound
seqD1: the sequence of the first displacer
seqD2: the sequence of the second displacer
Cd1: the concentration of the first displacer in solution (M)
Cd2: the concentration of the second displacer in solution (M)
saltConc: the salt concentration of the solution (M)

The rest of the packages involved are merely those needed to make the code mentioned above function.
The only thing in these packages that ever needs to be changed are in the BindingStrength packages where:
Ts = 65; % number of bases in the poly T chain
should be changed to reflect the total length (in bases) of your ssDNA strand from the surface of the particle to the end of the sticky end
N_DNA = 5000; %strands per particle
should be changed to reflect the DNA coverage of your particle (note: this code works for relatively low coverage <~20,000 strands per particle)
