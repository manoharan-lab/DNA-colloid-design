function f=twoSequenceMeltingCurvePlot(seqA, saltConc)

    tempStep = 0.25;
    Tmin = 40;
    Tmax = 48;

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
        [del_F(tempi),info] = twoSequenceBindingStrength(seqA, Xtemp(tempi), saltConc); %kcal/mol
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
    
    figure
    %Xtemp
    %singletFraction
    plot(Xtemp,singletFraction,'LineWidth',2)
    set(gca,'fontsize',16,'fontname','Helvetica','fontweight','b')
    xlabel('Temperature (C)')
    ylabel('singlet fraction')
%    title(['salt: ' num2str(saltConc) ', length: ' num2str(length(seqA))])
    
    printer = [' TmeltAB= ', num2str(info(2))];


end