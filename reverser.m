function [ revSeq ] = reverser( seq )
%UNTITLED3 Summary of this function goes here

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

