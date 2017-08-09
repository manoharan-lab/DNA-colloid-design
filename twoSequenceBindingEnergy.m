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

function [sequenceEnergy,meltingTemp]=twoSequenceBindingEnergy(sequence, TempC, saltConc)

%MAKE SURE SEQUENCE IS ENTERED 5' TO 3'

%sequence1
sequence1=lower(sequence);
if isempty(sequence1)
    sequenceEnergy=0;
    meltingTemp=0;
else
    Temp = TempC + 273.15;
    R = 1.9872; %cal/(K mol)
    Ct = 30*150*10^(-6); %This is the approximate actual concentration of linkers in overlap region (Molar)
    [H_pair,S_pair,G_pair] = basePairBindingEnergies(TempC, saltConc);
    
    %initiation energy
    sequenceEnergy = G_pair(1);
    del_H = H_pair(1);
    del_S = S_pair(1);
    
    RepeatGC = 0; %keeps track of number of gg and cc pairs
    ii=0;
    for ii=1:(length(sequence1)-1)
        if (sequence1(ii)=='a'&&sequence1(ii+1)=='a')||(sequence1(ii)=='t'&&sequence1(ii+1)=='t')
            bob='aa or tt';
            sequenceEnergy = sequenceEnergy + G_pair(2);
            del_H = del_H + H_pair(2);
            del_S = del_S + S_pair(2);
        elseif (sequence1(ii)=='a'&&sequence1(ii+1)=='t')
            bob='at';
            sequenceEnergy = sequenceEnergy + G_pair(3);
            del_H = del_H + H_pair(3);
            del_S = del_S + S_pair(3);
        elseif (sequence1(ii)=='t'&&sequence1(ii+1)=='a')
            bob='ta';
            sequenceEnergy = sequenceEnergy + G_pair(4);
            del_H = del_H + H_pair(4);
            del_S = del_S + S_pair(4);
        elseif (sequence1(ii)=='c'&&sequence1(ii+1)=='a')||(sequence1(ii)=='t'&&sequence1(ii+1)=='g')
            bob='ca or tg';
            sequenceEnergy = sequenceEnergy + G_pair(5);
            del_H = del_H + H_pair(5);
            del_S = del_S + S_pair(5);
        elseif (sequence1(ii)=='g'&&sequence1(ii+1)=='t')||(sequence1(ii)=='a'&&sequence1(ii+1)=='c')
            bob='gt or ac';
            sequenceEnergy = sequenceEnergy + G_pair(6);
            del_H = del_H + H_pair(6);
            del_S = del_S + S_pair(6);
        elseif (sequence1(ii)=='c'&&sequence1(ii+1)=='t')||(sequence1(ii)=='a'&&sequence1(ii+1)=='g')
            bob='ct or ag';
            sequenceEnergy = sequenceEnergy + G_pair(7);
            del_H = del_H + H_pair(7);
            del_S = del_S + S_pair(7);
        elseif (sequence1(ii)=='g'&&sequence1(ii+1)=='a')||(sequence1(ii)=='t'&&sequence1(ii+1)=='c')
            bob='ga or tc';
            sequenceEnergy = sequenceEnergy + G_pair(8);
            del_H = del_H + H_pair(8);
            del_S = del_S + S_pair(8);
        elseif (sequence1(ii)=='c'&&sequence1(ii+1)=='g')
            bob='cg';
            sequenceEnergy = sequenceEnergy + G_pair(9);
            del_H = del_H + H_pair(9);
            del_S = del_S + S_pair(9);
        elseif (sequence1(ii)=='g'&&sequence1(ii+1)=='c')
            bob='gc';
            sequenceEnergy = sequenceEnergy + G_pair(10);
            del_H = del_H + H_pair(10);
            del_S = del_S + S_pair(10);
        elseif (sequence1(ii)=='g'&&sequence1(ii+1)=='g')||(sequence1(ii)=='c'&&sequence1(ii+1)=='c')
            bob='gg or cc';
            RepeatGC = RepeatGC + 1;
            sequenceEnergy = sequenceEnergy + G_pair(11);
            del_H = del_H + H_pair(11);
            del_S = del_S + S_pair(11);
        else
            bob='something is wrong';
        end
    end
    
    
    %Terminal AT penalties for start and end of sequence
    if (sequence1(1)=='a')||(sequence1(1)=='t')
        bob='beginning at';
        sequenceEnergy=sequenceEnergy+G_pair(12);
        del_H = del_H + H_pair(12);
        del_S = del_S + S_pair(12);
    end
    
    if (sequence1(length(sequence1))=='a')||(sequence1(length(sequence1))=='t')
        bob='terminal at';
        sequenceEnergy=sequenceEnergy+G_pair(12);
        del_H = del_H + H_pair(12);
        del_S = del_S + S_pair(12);
    end
    
    meltingTemp=del_H*1000/(del_S+R*log(Ct/4))-273.15;
    finalH=del_H;
    finalS=del_S;
    
    RepeatGC;
end

end
