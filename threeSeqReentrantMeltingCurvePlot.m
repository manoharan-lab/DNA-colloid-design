function f=threeSeqReentrantMeltingCurvePlot(seqA, seqD1, Cd1, saltConc)

    tempStep = 0.25;
    Tmin = 25;
    Tmax = 70;

    del_F = zeros(1,(Tmax-Tmin)/tempStep);
    singletFraction = zeros(1,(Tmax-Tmin)/tempStep);
    Xtemp = (Tmin:tempStep:Tmax);
    R = 1.9872/1000; %kcal/(K mol)
    Cp = 0.1; %areal concentration of particles (um)^-2 (This is a normal concentration to use in lab)
    
    %things to understand later
    l = 0.015;%range of interaction (um) (according to Ben)
    z = 3; %coordination of each particle (according to Ben)
    
    check=0;
    for tempi=1:((Tmax-Tmin)/tempStep+1)
        [del_F(tempi),info] = threeSeqReentrantBindingStrength(seqA, seqD1, Cd1, Xtemp(tempi), saltConc); %kcal/mol
        %del_F(tempi)
        RT = R*(Xtemp(tempi)+273);
        K = (l/2)^2*exp(-(z/2)*del_F(tempi)/RT);
        singletFraction(tempi) = (1/(2*(K*Cp)^2)) * (1+2*K*Cp-sqrt(1+4*K*Cp));
        if (check==0 && singletFraction(tempi)>=0.5)
            MeltTemp = Xtemp(tempi)
            check = 1;
        end
        %singletFraction(tempi)
    end

    %Xtemp
    %delG
    
%      figure
%      plot(Xtemp,del_F)
%      set(gca,'fontsize',16,'fontname','Helvetica','fontweight','b')
%      xlabel('Temperature (C)')
%      ylabel('Binding Strength (kT)')
    

    %Xtemp
    %singletFraction
%    dataX = 45;
%    dataY = 0.5;
    dataX = 44;%31;
    dataY = 0.5;
    figure('Color','white')
    hold all
    plot(Xtemp,(singletFraction),'b','LineWidth',2)
    plot(dataX,dataY,'db','MarkerSize',12,'MarkerFaceColor','b')
    set(gca,'fontsize',16,'fontname','Helvetica','fontweight','b')
    axis([Tmin Tmax 0 1])
    xlabel('Temperature (C)')
%    ylabel('binding probability')
%    title(['salt: ' num2str(saltConc) ', A/B length: ' num2str(length(seqA)) ', D1 length: ' num2str(length(seqD1))])
    hold off
    printer = [' TmeltAB= ', num2str(info(2)), ', TmeltAD1 = ', num2str(info(4))]

    
    

end