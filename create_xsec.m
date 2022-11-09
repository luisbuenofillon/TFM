function create_xsec(ci)


    T_fueli = 5.592504e+02;
    T_modi = 5.592097e+02;
    density_i = 7.523201e+02/1000;
     %T_fueli = 618.3;
      %T_modi = 306.6;
    %  density_i = 712.5;
    fileID = fopen('../Prueba_2/NEACRP_th_wr.xsec','w');
    intro;
    for i=1:length(intro1)
        fprintf(fileID,intro1{i});
        fprintf(fileID,'\n');
    end
    %% Boron Concentration
    format long
    c0=1200.2;
    T_mod0 = 306.6+273.15;
    density_0=0.7125;
    T_fuel0 = 618.3+273.15;
    initial_data_OCR;
    for k=1:length(sigma)
        sigmak = sigma{k};
        if isempty(boron{k})==0
            sigmak_c = boron{k}*(ci-c0);
            sigmak_Tm = T_mod{k}*(T_modi-T_mod0);
            sigmak_rho = density{k}*(density_i-density_0);
            sigmak_Tf = T_fuel{k}*(sqrt(T_fueli)-sqrt(T_fuel0));
            sigma{k}=sigmak + sigmak_c + sigmak_Tm + sigmak_rho + sigmak_Tf;  
        end
    end
    sigma{12} = sigma{1} + sigma_CA1_18 + sigma_driver;
    sigma{13} = sigma{1} + sigma_CA2_18+ sigma_driver;
    sigma{14} = sigma{4} + sigma_CA1;
    sigma{15} = sigma{6} + sigma_CA2;

    for i=1:15
        s = sigma{i};
        s1=s(1,:);
        s2=s(2,:);
        for j=1:length(s1)
            if s1(j)< 0.0
                s1(j) = 0.0;
            end
        end
        for j=1:length(s2)
            if s2(j)< 0.0
                s2(j) = 0.0;
            end
        end
        fprintf(fileID,'# Material %d\n%d 	%.5e 	%.5e 	%.5e 	%.5e 	%.5e\n 	%.5e 	%.5e 	%.5e 	%.5e\n',i,i,s1,s2(1:end-1));
    end
    
    fclose(fileID);
fprintf('ci:  %f    ',ci)
end
