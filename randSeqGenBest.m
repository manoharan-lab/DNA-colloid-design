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

function output=randSeqGenBest(len,seq,max,num)

timei = cputime;
randomTemp = 273+30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: to put in multiple sequences use the format:
% sequences={'atcg'; 'tggcc'; 'atataatatatatat'}
% This allows for sequences of different lengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

possibles=[''];
melties=[];
space=[''];
countOut=[];
for ii=1:num
    [test_seq,counter] = randSeqGenExcMax(len,seq,max);
    [junk2,test_melt] = twoSequenceBindingEnergy(test_seq, randomTemp, 0.25);
    if isempty(test_seq)
        messageOut ='empty sequence';
    else
        messageOut ='not empty';
        possibles=[possibles;test_seq];
        melties = [melties;test_melt];
        space = [space;'   '];
        countOut=[countOut;counter];
    end
end

[out_len,junk]=size(possibles);
if out_len>0
   out1 = [possibles  space num2str(countOut) space num2str(melties)];
   output = out1;
else
   output='sorry no available sequences'; 
end

%space num2str(melties)

timef=cputime-timei
end
