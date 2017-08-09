% Copyright 2017, Vinothan N. Manoharan, Emily W. Gehrels
%
% This file is part of DNA-colloid-design.
%
% DNA-colloid-design is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% DNA-colloid-design is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with DNA-colloid-design.  If not, see <http://www.gnu.org/licenses/>.

function [bindingStrength,info]=fourSeqReentrantBindingStrength(seqA, seqD1, seqD2, Cd1, Cd2, TempC, saltConc)

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
f1 = 1;
f2 = 1;

[del_G_AB,Tmelt_AB] = twoSequenceBindingEnergy(seqA, TempC, saltConc); %kcal/mol
[del_G_AD1,Tmelt_AD1] = twoSequenceBindingEnergy(seqD1, TempC, saltConc); %kcal/mol
[del_G_BD2,Tmelt_BD2] = twoSequenceBindingEnergy(seqD2, TempC, saltConc); %kcal/mol

Kappa1 = exp(-f1*del_G_AB/RT);
Kappa2 = exp((f1*del_G_AB-f2*del_G_AD1)/RT);
Kappa3 = exp((f1*del_G_AB-f2*del_G_BD2)/RT);

N_DNA = 5000; %strands per particle
V_shell = (4/3)*pi * ((a + L)^3 - a^3); %Volume of DNA shell around particle (dm^3)
CAo = (N_DNA/V_shell) * (1/Nav) %Concentration of A DNA around each particle (mol/dm^3 -> Molar)
CBo = (1)*(N_DNA/V_shell) * (1/Nav); %Concentration of B DNA around each particle (mol/dm^3 -> Molar)
b = CBo/CAo;

Co = 1; %Reference concentration from NN model (Molar)
%Vi = (pi/12) * (4*(a+L) + (2*a+L)) * (2*(a+L) - (2*a+L))^2; %Overlap Volume of two particle's DNA shells (dm^3)
Vi = (pi/12) * (6*a*L^2 + 5*L^3); %Overlap Volume of two particle's DNA shells (dm^3)
Neff1 = CAo*Vi*Nav; %Avg number of A DNA strands in the gap 
Neff2 = CBo*Vi*Nav; %Avg number of B DNA strands in the gap 
Neff = min(Neff1,Neff2); %minimum of Neff1 and Neff2

Kt = Kappa1*CAo*Co/((1+Kappa1*Kappa2*Cd1/Co)*(1+Kappa1*Kappa3*Cd2/Co));
%p = 1 - (sqrt(1+4*Kt) - 1)/(2*Kt);
p = ((b+1)*Kt + 1 - sqrt((b-1)^2*(Kt)^2 + 2*(b+1)*Kt + 1))/(2*b*Kt);
bindingStrength = log((1-p)^Neff); %kT


info=[del_G_AB,Tmelt_AB,del_G_AD1,Tmelt_AD1,del_G_BD2,Tmelt_BD2];

end
