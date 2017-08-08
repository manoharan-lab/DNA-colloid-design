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