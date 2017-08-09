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

function [bindingStrength,info]=twoSequenceBindingStrength(sequence1, TempC, saltConc)

%1 Liter = 1 dm^3

Temp = TempC + 273.15;

R = 1.9872/1000; %kcal/(K mol)
Nav = 6.02*10^23; %Avogadro's number
a = 5*10^(-6); %radius of the particles (dm)
Ts = 65; % number of bases in the DNA chain
Kuhn = 3*10^(-8); %Kuhn length (dm) (SHOULD THIS BE MORE LIKE 4.5?)
ssBaseLength = 6.5 *10^(-9); %(dm)
n = Ts*ssBaseLength/Kuhn;
L = sqrt(n)*Kuhn; %Length of the DNA (dm) THIS NEEDS DEPENDENCE ON SALT CONCENTRATION
RT = R*Temp;
f=1;

[del_G,Tmelt] = twoSequenceBindingEnergy(sequence1, TempC, saltConc); %kcal/mol
Kappa = exp(-f*del_G/RT);

N_DNA = 6500; %strands per particle
V_shell = (4/3)*pi * ((a + L)^3 - a^3); %Volume of DNA shell around particle (dm^3)
CAo = (1)*(N_DNA/V_shell) * (1/Nav); %Concentration of A DNA around each particle (mol/dm^3 -> Molar)
CBo = (1)*(N_DNA/V_shell) * (1/Nav); %Concentration of B DNA around each particle (mol/dm^3 -> Molar)
b = CBo/CAo; %which way?!?

Co = 1; %Reference concentration from NN model (Molar)
%Vi = (pi/12) * (4*(a+L) + (2*a+L)) * (2*(a+L) - (2*a+L))^2; %Overlap Volume of two particle's DNA shells (dm^3)
Vi = (pi/12) * (6*a*L^2 + 5*L^3); %Overlap Volume of two particle's DNA shells (dm^3)
Neff1 = CAo*Vi*Nav; %Avg number of A DNA strands in the gap 
Neff2 = CBo*Vi*Nav; %Avg number of B DNA strands in the gap 
Neff = min(Neff1,Neff2); %minimum of Neff1 and Neff2

%p2 = 1 + (Co - sqrt(Co^2 + 4*Kappa*CAo*Co))/(2*Kappa*CAo);
Kt = Kappa*CAo/Co;
p = ((b+1)*Kt + 1 - sqrt((b-1)^2*(Kt)^2 + 2*(b+1)*Kt+1))/(2*b*Kt);
bindingStrength = log((1-p)^Neff); %kT

%bindingStrength2 = N_DNA*(Vi/V_shell) * log((-Co + sqrt(Co^2+4*(Kappa*N_DNA*Co)/(V_shell*Nav))) / (2*Kappa*N_DNA/(V_shell*Nav))); %sanity check
%bindingStrength3 = log((1-p2)^Neff) %kT
info=[del_G,Tmelt];

end
