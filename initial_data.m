%initial data

%% Secciones 1 - 11

%comp_num   1    !interior reflector
%------------------------------------------------------------------------------
 sigma1=[ 5.320580e-02  3.732790e-04  0.000000e+00  0.000000e+00  2.645540e-02;
            3.864060e-01  1.772150e-02  0.000000e+00  0.000000e+00 1];
 c1  =[6.118330e-08  1.877310e-07  0.000000e+00  0.000000e+00  7.914570e-10;
            5.175350e-06  1.026350e-05  0.000000e+00  0.000000e+00 1];
 dxs_ddm1   =[7.457560e-02  2.076880e-04  0.000000e+00  0.000000e+00  3.713100e-02;
            5.336340e-01  7.584210e-03  0.000000e+00  0.000000e+00 1];
 %comp_num   2    !exterior reflector
%------------------------------------------------------------------------------
 sigma2= [2.956090e-01  1.187820e-03  0.000000e+00  0.000000e+00  2.316130e-02;
            2.459310e+00  2.526180e-01  0.000000e+00  0.000000e+00 1];
c2=[   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00;
            7.761840e-04  8.446950e-05  0.000000e+00  0.000000e+00 1];

 %comp_num   3    !corner reflector
%------------------------------------------------------------------------------
 sigma3=[ 2.956090e-01  1.187820e-03  0.000000e+00  0.000000e+00  2.008080e-02;
            2.459310e+00  2.526180e-01  0.000000e+00  0.000000e+00 1];
c3=   [0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00;
            7.761840e-04  8.446950e-05  0.000000e+00  0.000000e+00 1];

 %comp_num   4    !fuel 1
%------------------------------------------------------------------------------
 sigma4 =[2.221170e-01  8.717740e-03  4.982770e-03  0.190224e-2  1.824980e-02;
            8.031400e-01  6.525500e-02  8.390260e-02  0.343581e-1 1];
c4   =[3.478090e-08  1.285050e-07 -1.120990e-09 -0.548361e-9 -1.085900e-07;
           -9.765100e-06  7.088070e-06 -2.430450e-06 -0.995273e-6 1];
 dxs_dtm4   =[-2.033100e-06  2.121910e-07  1.247090e-07  1.430354e-18  8.096760e-07;
           -1.086740e-04 -3.155970e-05 -4.164390e-05 -5.467221e-16 1];
 dxs_ddm4    =[1.356650e-01  1.551850e-03  9.206940e-04  1.023919e-14  2.931950e-02;
            9.926280e-01  2.526620e-02  2.477460e-02  3.252554e-13 1];
 dxs_dtf4   =[-3.091970e-05  3.497090e-05  6.401340e-07  7.154124e-18 -2.755360e-05;
           -1.372920e-04 -3.718060e-05 -5.630370e-05 -7.391879e-16 1];
 cdf4      =[  1.0069   0.9307   1.0034   0.9646   1.1040   1.4493   1.0096   1.1580];

 %comp_num   5    !fuel 2
%------------------------------------------------------------------------------
 sigma5=[ 2.219140e-01  9.061330e-03  5.576590e-03  0.214498e-2  1.800400e-02;
            7.955380e-01  7.233540e-02  9.986290e-02  0.408938e-1 1];
c5 =[  3.538260e-08  1.267090e-07 -1.678800e-09 -0.777980e-9 -1.069510e-07;
           -8.501690e-06  6.823110e-06 -2.724450e-06 -0.111566e-5 1];
 dxs_dtm5 =[  -1.980800e-06  2.260000e-07  1.351450e-07  1.568962e-18  8.584740e-07;
           -9.061500e-05 -3.214350e-05 -4.531020e-05 -5.948574e-16 1];
 dxs_ddm5  =[  1.357480e-01  1.614910e-03  9.641600e-04  1.081410e-14  2.926960e-02;
            9.819850e-01  2.866670e-02  3.149930e-02  4.135416e-13 1];
 dxs_dtf5 =[  -3.086070e-05  3.517980e-05  9.974310e-07  1.186847e-17 -2.767660e-05;
           -1.174810e-04 -3.770390e-05 -6.041550e-05 -7.931704e-16 1];
 cdf5    =[    1.0069   0.9307   1.0034   0.9646   1.1040   1.4493   1.0096   1.1580];

 %comp_num   6    !fuel 3
%------------------------------------------------------------------------------
sigma6=[ 2.217150e-01  9.384960e-03  6.150470e-03  0.237972e-2  1.776700e-02;
            7.892530e-01  7.892030e-02  1.146670e-01  0.469561e-1 1];
c6=[   3.598380e-08  1.249860e-07 -2.210380e-09 -0.996653e-9 -1.053740e-07;
           -7.462510e-06  6.597980e-06 -2.958830e-06 -0.121164e-5 1];
 dxs_dtm6 =[  -1.924340e-06  2.399390e-07  1.490840e-07  1.754220e-18  9.034940e-07;
           -7.627860e-05 -3.237760e-05 -4.784750e-05 -6.281741e-16 1];
 dxs_ddm6 =[   1.358270e-01  1.680150e-03  1.014100e-03  1.147707e-14  2.921540e-02;
            9.722670e-01  3.195710e-02  3.810970e-02  5.003282e-13 1];
 dxs_dtf6 =[  -3.091650e-05  3.538410e-05  1.418470e-06  1.742686e-17 -2.783900e-05;
           -1.013370e-04 -3.775580e-05 -6.309600e-05 -8.283631e-16 1];
 cdf6  =[      1.0069   0.9307   1.0034   0.9646   1.1040   1.4493   1.0096   1.1580];

 %comp_num   7    !fuel 4
%------------------------------------------------------------------------------
 sigma7=[ 2.220390e-01  9.316920e-03  5.550100e-03  0.213629e-2  1.713810e-02;
            7.762300e-01  7.963280e-02  9.855760e-02  0.403596e-1 1];
c7=[   3.378060e-08  1.198690e-07 -1.713230e-09 -0.777980e-9 -1.008730e-07;
           -6.737440e-06  6.293100e-06 -2.553590e-06 -0.104561e-5 1];
 dxs_dtm7 =[  -2.696340e-06  2.485300e-07  1.407730e-07  1.568962e-18  7.013110e-07;
           -7.624350e-05 -3.001190e-05 -4.202020e-05 -5.516691e-16 1];
 dxs_ddm7 =[   1.310330e-01  1.683970e-03  9.819510e-04  1.081410e-14  2.824890e-02;
            9.346970e-01  3.142400e-02  3.515880e-02  4.617153e-13 1];
 dxs_dtf7  =[ -3.137460e-05  3.486990e-05  9.454310e-07  1.186847e-17 -2.735500e-05;
           -1.082710e-04 -3.727480e-05 -5.796620e-05 -7.608485e-16 1];
 cdf7     =[   1.0069   0.9307   1.0034   0.9646   1.1040   1.4493   1.0096   1.1580];

 %comp_num   8    !fuel 5
%------------------------------------------------------------------------------
 sigma8=[ 2.220830e-01  9.400320e-03  5.540830e-03  0.213318e-2  1.685010e-02;
            7.699690e-01  8.210870e-02  9.800590e-02  0.401338e-1 1];
c8 =[  3.324950e-08  1.175850e-07 -1.724210e-09 -0.793325e-9 -9.885780e-08;
           -6.197250e-06  6.119040e-06 -2.488800e-06 -0.101904e-5 1];
 dxs_dtm8 =[  -3.079050e-06  2.618540e-07  1.432350e-07  1.678974e-18  6.173800e-07;
           -7.333970e-05 -2.919290e-05 -4.077010e-05 -5.352614e-16 1];
 dxs_ddm8 =[   1.293790e-01  1.719720e-03  9.884370e-04  1.113221e-14  2.788950e-02;
            9.181710e-01  3.247150e-02  3.632510e-02  4.770784e-13 1];
 dxs_dtf8 =[  -3.155030e-05  3.472740e-05  9.260780e-07  1.089348e-17 -2.723810e-05;
           -1.055210e-04 -3.718080e-05 -5.711080e-05 -7.495752e-16 1];
 cdf8   =[     1.0069   0.9307   1.0034   0.9646   1.1040   1.4493   1.0096   1.1580];

 %comp_num   9    !fuel 6
%------------------------------------------------------------------------------
 sigma9=[ 2.221270e-01  9.482860e-03  5.531370e-03  0.213003e-2  1.656260e-02;
            7.638130e-01  8.459120e-02  9.741090e-02  0.398902e-1 1];
c9 =[  3.272010e-08  1.153190e-07 -1.735020e-09 -0.796917e-9 -9.684890e-08;
           -5.682200e-06  5.947110e-06 -2.422400e-06 -0.991817e-6 1];
 dxs_dtm9 =[  -3.538770e-06  2.743130e-07  1.460190e-07  1.716653e-18  5.165470e-07;
           -7.137110e-05 -2.830410e-05 -3.943190e-05 -5.176891e-16 1];
 dxs_ddm9 =[   1.276820e-01  1.749890e-03  9.951750e-04  1.122094e-14  2.752020e-02;
            9.012930e-01  3.359450e-02  3.744990e-02  4.919001e-13 1];
 dxs_dtf9 =[  -3.172810e-05  3.460260e-05  9.058020e-07  1.061656e-17 -2.711690e-05;
           -1.025250e-04 -3.702010e-05 -5.615430e-05 -7.369691e-16 1];
 cdf9 =[       1.0069   0.9307   1.0034   0.9646   1.1040   1.4493   1.0096   1.1580];

%comp_num  10    !fuel 7
%------------------------------------------------------------------------------
sigma10 =[ 2.218360e-01  9.637200e-03  6.123820e-03  0.237097e-2  1.690430e-02;
            7.707050e-01  8.611870e-02  1.132410e-01  0.463724e-1 1];
c10 =[   3.438590e-08  1.181860e-07 -2.243350e-09 -0.996653e-9 -9.933120e-08;
           -5.868980e-06  6.084430e-06 -2.776570e-06 -0.113696e-5 1];
 dxs_dtm10= [   -2.639070e-06  2.642890e-07  1.558580e-07  1.754220e-18  7.443200e-07;
           -6.395540e-05 -3.035090e-05 -4.444310e-05 -5.834832e-16 1];
 dxs_ddm10=[    1.311160e-01  1.755280e-03  1.035220e-03  1.147707e-14  2.818770e-02;
            9.249250e-01  3.498530e-02  4.206930e-02  5.523868e-13 1];
 dxs_dtf10 =[  -3.141920e-05  3.506370e-05  1.356420e-06  1.742686e-17 -2.750490e-05;
           -9.388860e-05 -3.714030e-05 -6.050520e-05 -7.942518e-16 1];
 cdf10  =[      1.0069   0.9307   1.0034   0.9646   1.1040   1.4493   1.0096   1.1580];

 %comp_num  11    !fuel 8
%------------------------------------------------------------------------------
sigma11=[ 2.218780e-01  9.719370e-03  6.114440e-03  0.236781e-2 1.661750e-02;
            7.647040e-01  8.854880e-02  1.126350e-01  0.461246e-1 1];
c11=[   3.385590e-08  1.159170e-07 -2.253690e-09 -0.101112e-8 -9.732910e-08;
           -5.383450e-06  5.916970e-06 -2.707800e-06 -0.110878e-5 1];
 dxs_dtm11=[   -3.021470e-06  2.790600e-07  1.588140e-07  1.885280e-18  6.595210e-07;
           -6.169840e-05 -2.956260e-05 -4.315880e-05 -5.666222e-16 1];
 dxs_ddm11=[    1.294630e-01  1.794990e-03  1.042910e-03  1.185343e-14  2.782590e-02;
            9.084560e-01  3.610320e-02  4.332150e-02  5.688569e-13 1];
 dxs_dtf11=[   -3.159080e-05  3.491190e-05  1.333360e-06  1.627694e-17 -2.738350e-05;
           -9.171260e-05 -3.699090e-05 -5.962840e-05 -7.827158e-16 1];
 cdf11 =[        1.0069   0.9307   1.0034   0.9646   1.1040   1.4493   1.0096   1.1580];

 
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
 %% Calculo 12 y 13


% sigma12=[ 2.7396e-03 0.0000e+00 0.0000e+00 0 2.4796e-02;
%  3.3753e-02 0.0000e+00 0.0000e+00 0 1];
% c12=[];
% sigma13=[ 2.4190e-03 0.0000e+00 0.0000e+00 0 2.5209e-02;
%  3.3753e-02 0.0000e+00 0.0000e+00 0 1];
% c13=[];
sigma12 = sigma1 + sigma_CA1_18+sigma_driver;
c12=[];
sigma13 = sigma1 + sigma_CA2_18 + sigma_driver;
c13=[];
%% Calculo secciones 14 y 15



sigma14 = sigma4 + sigma_CA1;
c14=[];
sigma15 = sigma6 + sigma_CA2;
c15=[];

sigma = {sigma1,sigma2,sigma3,sigma4,sigma5,sigma6,sigma7,sigma8,sigma9,sigma10,sigma11,sigma12,sigma13,sigma14,sigma15};
c = {c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15};
% for i=1:15
%     s = sigma{i};
%     s1=s(1,:);
%     s2=s(2,:);
%     fprintf('# Material %d\n%d 	%.5e 	%.5e 	%.5e 	%.5e 	%.5e\n 	%.5e 	%.5e 	%.5e 	%.5e\n',i,i,s1,s2(1:end-1))
% end

