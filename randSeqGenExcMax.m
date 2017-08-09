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

function [randomSequence,tester]=randSeqGenExcMax(len,sequences, max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: to put in multiple sequences use the format:
% sequences={'atcg'; 'tggcc'; 'atataatatatatat'}
% This allows for sequences of different lengths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate codonAlphabet and holder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Codon=codonAlphabet();
codonHold=Codon;
capCodonHold=capCodonAlphabet();
usedCodons=[''];

counter=0;

%generate a random sequence and test each codon against used ones
%if there are more than 'max' reused codons reject the sequence 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Break sequences into codons and put them and their
% reverse complement in the used codon library
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[numOfSeq,junk]=size(sequences);
for ii=1:numOfSeq
    seq_Hold=sequences(ii,:);
    seq = char(seq_Hold);
    for xx=1:(length(seq)-2)
        segments(xx,1:3) = seq(xx:(xx+2));
        compSegments(xx,(1:3))=reverser(seq(xx:(xx+2)));
    end
    
    markAsUsed=[];
    for yy=1:(length(seq)-2)
        [~,indx]=ismember(segments(yy,:),codonHold,'rows');
        if indx==0
        else
            markAsUsed = [markAsUsed,indx];
        end
        [~,indx]=ismember(compSegments(yy,:),codonHold,'rows');
        if indx==0
        else
            markAsUsed = [markAsUsed,indx];
        end
    end
    
    usedCodons(markAsUsed,:)=codonHold(markAsUsed,:);
    %codonHold(markAsUsed,:)=[];
end
%codonHold
%usedCodons


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First codon randomly generated, put in randomSequence and it (w/ comp)
% added to usedCodons and replaced with caps. Last two bases held.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomSequence=[''];

randNum=randi(length(codonHold),1,1);
randomSequence(1:3)=codonHold(randNum,:);
revCodon(1:3)=reverser(codonHold(randNum,:));

numberOfBases=3;
partialHolder= codonHold(randNum,[2,3]);

markAsUsed = [];
[~,indx]=ismember(randomSequence,codonHold,'rows');
[~,indx2]=ismember(randomSequence,usedCodons,'rows');
[~,indx3]=ismember(revCodon,codonHold,'rows');
[~,indx4]=ismember(revCodon,usedCodons,'rows');

if indx==0
else
    markAsUsed=[markAsUsed,indx];
end

if (indx2==0||indx4==0)
else
    randomSequence(1:3)=capCodonHold(randNum,:);
    counter=counter+1;
end

if indx3==0
else
    markAsUsed=[markAsUsed,indx3];
end


usedCodons(markAsUsed,:)=codonHold(markAsUsed,:);
%randomSequence
%usedCodons


while numberOfBases<len
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find all possible next codons (which have same first two bases as
    % last two bases of previous) if none available exit routine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    twoBases=[];
    for jj=1:size(codonHold,1)
        if codonHold(jj,[1,2])==partialHolder
            twoBases=[twoBases,jj];
        end
    end
    
    if isempty(twoBases)==1
        %%exit routine, return message, and clear randomSequence
        messageOut ='unable to complete sequence';
        numberOfBases=len+1;
        randomSequence=[];
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Randomly choose next codon from list of possibles, put final base
        % in randomSequence and removed (w/ complement) from library.
        % Last two bases held.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        randNum=randi(length(twoBases),1,1);
        nextBaseNumber=twoBases(randNum);
        randomSequence(numberOfBases+1)=codonHold(nextBaseNumber,3);
        revCodon(1:3)=reverser(codonHold(nextBaseNumber,:));
        
        numberOfBases=numberOfBases+1;
        partialHolder= codonHold(nextBaseNumber,[2,3]);
        
        
        markAsUsed = [];
%         randomSequence(numberOfBases-2:numberOfBases)
        [~,indx]=ismember(lower(randomSequence(numberOfBases-2:numberOfBases)),codonHold,'rows');
        [~,indx2]=ismember(revCodon,codonHold,'rows');
        [~,indx3]=ismember(lower(randomSequence(numberOfBases-2:numberOfBases)),usedCodons,'rows');
        [~,indx4]=ismember(revCodon,usedCodons,'rows');
        
        if indx==0
        else
            markAsUsed=[markAsUsed,indx];
        end
        
        if indx2==0
        else
            markAsUsed=[markAsUsed,indx2];
        end
        
        if (indx3==0||indx4==0)
        else
            randomSequence(numberOfBases-2:numberOfBases)=capCodonHold(indx,:);
            counter=counter+1;
        end

        usedCodons(markAsUsed,:)=codonHold(markAsUsed,:);
        
    end    
    
end

if counter>max
    randomSequence=[''];
    tester=counter;
else
    tester=counter;
end

end

