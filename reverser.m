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

function [ revSeq ] = reverser( seq )
%this takes a sequence and finds it's reversed complement


leng=length(seq);
ii=0;
jj=0;
for ii=1:leng
    kk=[leng:-1:1];
    if seq(ii)=='a'
        revSeq(kk(ii))='t';
    elseif seq(ii)=='t'
        revSeq(kk(ii))='a';
    elseif seq(ii)=='c'
        revSeq(kk(ii))='g';
    elseif seq(ii)=='g'
        revSeq(kk(ii))='c';
    end
end

end

