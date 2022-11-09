%% CA1 CA2 and driver
sigma_CA1 = [0.37322e-2  0.24777e-2 -0.102786e-3 -0.377989e-4 -0.319253e-2;
            -0.219926e-1 0.255875e-1 -0.282319e-2 -0.115483e-2 1];
        
sigma_CA2 = [0.37409e-2  0.242926e-2 -0.122634e-3 -0.459250e-4 -0.314239e-2;
            -0.167503e-1 0.256478e-1 -0.328086e-2 0.134262e-2 1];

sigma_CA1_18 = [0.37322e-2  0.24777e-2 -0.102786e-3 0 0;
            -0.219926e-1 0.255875e-1 0 0 1];
        
sigma_CA2_18 = [0.37409e-2  0.242926e-2 -0.122634e-3  0 0;
            -0.167503e-1 0.256478e-1 0 0 1];
sigma_driver =[-0.697102e-02  0.879034e-04 0.000000e+00 0.000000e+00 -0.119034e-02;
-0.113498E-01 0.170043e-02 0.000000e+00 0.000000e+00 1];
%% Components 1 to 11
file= table2cell(readtable('initial_data.txt'));
data=cell(11,1);
for i=1:11
    data0=zeros(1,50);
    for j = 1:50
        k=(i-1)*51 + j;
        data0(j)=file{k};
    end
    data{i} = reshape(data0,[10,5]);
end

%% Create Variables
sigma = cell(15,1);
boron = cell(15,1);
T_mod = cell(15,1);
density = cell(15,1);
T_fuel = cell(15,1);
for i=1:11
    datai = data{i};
    sigma{i} = [datai(1,[1 3 4 5 2]);datai(2,:)];
    boron{i} = [datai(3,[1 3 4 5 2]);datai(4,:)];
    T_mod{i} = [datai(5,[1 3 4 5 2]);datai(6,:)];
    density{i} =[datai(7,[1 3 4 5 2]);datai(8,:)];
    T_fuel{i} =[datai(9,[1 3 4 5 2]);datai(10,:)];
end
sigma{12} = sigma{1} + sigma_CA1_18+sigma_driver;
sigma{13} = sigma{1} + sigma_CA2_18 + sigma_driver;
sigma{14} = sigma{4} + sigma_CA1;
sigma{15} = sigma{6} + sigma_CA2;