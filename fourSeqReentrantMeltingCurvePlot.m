function f=fourSeqReentrantMeltingCurvePlot(seqA, seqD1, seqD2, Cd1, Cd2, saltConc)

    tempStep = 0.25;
    Tmin = 0;
    Tmax = 80;

    del_F = zeros(1,(Tmax-Tmin)/tempStep);
    singletFraction = zeros(1,(Tmax-Tmin)/tempStep);
    Xtemp = (Tmin:tempStep:Tmax);
    R = 1.9872/1000; %kcal/(K mol)
    Cp = 0.1; %areal concentration of particles (um)^-2 (This is a normal concentration to use in lab)
    
    %things to understand later
    l = 0.015;%range of interaction (um) (according to Ben)
    z = 3; %coordination of each particle (according to Ben)
    
check1=0;
check2=0;
    for tempi=1:((Tmax-Tmin)/tempStep+1)
        [del_F(tempi),info] = fourSeqReentrantBindingStrength(seqA, seqD1, seqD2, Cd1, Cd2, Xtemp(tempi), saltConc); %kcal/mol
        %del_F(tempi)
        RT = R*(Xtemp(tempi)+273);
        K = (l/2)^2*exp(-(z/2)*del_F(tempi)/RT);
        singletFraction(tempi) = (1/(2*(K*Cp)^2)) * (1+2*K*Cp-sqrt(1+4*K*Cp));
        if (check1==0 && check2==0 && singletFraction(tempi)<=0.5)
            MeltTemp1 = Xtemp(tempi)
            check1 = 1;
        end
        if (check1==1 && check2==0 && singletFraction(tempi)>=0.5)
            MeltTemp2 = Xtemp(tempi)
            check2 = 1;
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
    
    %figure
    %Xtemp
    %singletFraction
    dataX = [31,50];
    dataY = [0.5,0.5];
    
    figure('Color','white')
    hold on
    plot(Xtemp,(1-singletFraction),'g','LineWidth',2)
    plot(dataX,dataY,'dg','MarkerSize',12,'MarkerFaceColor','g')
    set(gca,'fontsize',16,'fontname','Helvetica','fontweight','b')
    xlabel('Temperature (C)')
    ylabel('binding probability')
    %title(['salt: ' num2str(saltConc) ', A/B length: ' num2str(length(seqA)) ', D1/D2 length: ' num2str(length(seqD1))])
    hold off
    
    
    printer = [' TmeltAB= ', num2str(info(2)), ', TmeltAD1 = ', num2str(info(4)),', TmeltBD2 = ', num2str(info(6))]


end