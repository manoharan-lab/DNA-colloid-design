function Codons=capCodonAlphabet()

    Codons=zeros(60,3);
    bases=['A', 'T', 'C', 'G'];
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