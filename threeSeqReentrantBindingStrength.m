function [bindingStrength,info]=threeSeqReentrantBindingStrength(seqA, seqD1, Cd1, TempC, saltConc)

%1 Liter = 1 dm^3

Temp = TempC + 273.15;

R = 1.9872/1000; %kcal/(K mol)
Nav = 6.02*10^23; %Avogadro's number
a = 5*10^(-6); %radius of the particles (dm)
Ts = 65; % number of bases in the poly T chain
Kuhn = 3*10^(-8); %Kuhn length (dm)
ssBaseLength = 6.5 *10^(-9); %(dm)
n = (65)*ssBaseLength/Kuhn;
L = sqrt(n)*Kuhn; %Length of the DNA (dm) THIS NEEDS DEPENDENCE ON SALT CONCENTRATION
RT = R*Temp;
f1 = 1.0; %1.2
f2 = 1.0; %1.0

[del_G_AB,Tmelt_AB] = twoSequenceBindingEnergy(seqA, TempC, saltConc); %kcal/mol
[del_G_AD1,Tmelt_AD1] = twoSequenceBindingEnergy(seqD1, TempC, saltConc); %kcal/mol

Kappa1 = exp(-f1*del_G_AB/RT);
Kappa2 = exp((f1*del_G_AB-f2*del_G_AD1)/RT);
N_DNA = 50000; %strands per particle 6500
V_shell = (4/3)*pi * ((a + L)^3 - a^3); %Volume of DNA shell around particle (dm^3)
CAo = (1)*(N_DNA/V_shell) * (1/Nav); %Concentration of DNA around each particle (mol/dm^3 -> Molar)
CBo = (1)*(N_DNA/V_shell) * (1/Nav); %Concentration of B DNA around each particle (mol/dm^3 -> Molar)
b = CBo/CAo; %b<1

Co = 1; %Reference concetration from NN model (Molar)
%Vi = (pi/12) * (4*(a+L) + (2*a+L)) * (2*(a+L) - (2*a+L))^2; %Overlap Volume of two particle's DNA shells (dm^3)
Vi = (pi/12) * (6*a*L^2 + 5*L^3); %Overlap Volume of two particle's DNA shells (dm^3)
NeffA = CAo*Vi*Nav; %Avg number of A DNA strands in the gap 
NeffB = CBo*Vi*Nav; %Avg number of B DNA strands in the gap 
Neff = min(NeffA,NeffB); %minimum of Neff1 and Neff2

Kt = (Kappa1*CAo/Co)/((1+Kappa1*Kappa2*Cd1/Co));
%p = 1 - (sqrt(1+4*Kt) - 1)/(2*Kt);
p = ((b+1)*Kt + 1 - sqrt((b-1)^2*(Kt)^2 + 2*(b+1)*Kt + 1))/(2*b*Kt); %This b on the bottom only appears when b<1 otherwise just drop it
bindingStrength = log((1-p)^Neff); %kT


info=[del_G_AB,Tmelt_AB,del_G_AD1,Tmelt_AD1];

end