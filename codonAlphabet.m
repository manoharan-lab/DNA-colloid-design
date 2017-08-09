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

function Codons=codonAlphabet()

    Codons=zeros(60,3);
    bases=['a', 't', 'c', 'g'];
    hold=npermutek(bases,3);
    threepeat=zeros(1,64);

    for ii=1:4^3
        if hold(ii,1)==hold(ii,2)
           if hold(ii,1)==hold(ii,3)
               threepeat(ii)=1;
           end
        end
    end
    
    

    Codons=hold(threepeat(:)==0,:);
end
