function [H_pair,S_pair,G_pair]=basePairBindingEnergies(TempC,saltConc)

% Temperature in Kelvin
% These equations work for salt concentrations between 50mM and 1.1M
Temp = TempC + 273.15; 

H_pair = [0.2, -7.6, -7.2, -7.2, -8.5, -8.4, -7.8, -8.2, -10.6, -9.8, -8.0, 2.2]; %in kcal/mol
S_pair_1M = [-5.7, -21.3, -20.4, -21.3, -22.7, -22.4, -21.0, -22.2, -27.2, -24.4, -19.9, 6.9]; %in e.u.
S_pair = S_pair_1M + 0.368*(1/2)*log(saltConc); % S&H pg 422
%G_pair=[1.96, -1, -0.88, -0.58, -1.45, -1.44, -1.28, -1.30, -2.17, -2.24, -1.84, 0.05]; in kcal/mol

G_pair = H_pair - Temp*S_pair/1000;


% 2 aa,tt- -7.6, -21.3, -1
% 3 at- -7.2, -20.4, -0.88
% 4 ta- -7.2, -21.3, -0.58
% 5 ca,tg- -8.5, -22.7, -1.45
% 6 gt,ac- -8.4, -22.4 -1.44
% 7 ct,ag- -7.8, -21.0, -1.28
% 8 ga,tc- -8.2, -22.2, -1.30
% 9 cg- -10.6, -27.2, -2.17
% 10 gc- -9.8, -24.4, -2.24
% 11 gg,cc- -8.0, -19.9, -1.84
% 12 initiation- 0.2, -5.7, 1.96
% 1 terminal at- 2.2, 6.9, 0.05


% These numbers from Santalucia, Hicks 2004 Annu. Rev. Biophys. Biomol. Struct. 
end