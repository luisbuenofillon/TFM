function  k0 = run_femffusion()

[~, out ]= system('cd ../../; ./femffusion.exe -f 3D_NEACRP_A1/Prueba_2/NEACRP_A1_20.prm');
index = strfind(out,'K0');
k0 = str2double(out(index+5:index+5+6));
if isnan(k0)
    k0=str2double(out(index+5:index+5+2));
end
fprintf('Femffusion run completed %f\n',k0)
end