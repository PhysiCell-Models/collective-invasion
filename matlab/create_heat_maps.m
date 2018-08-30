%Move this and read_microenvironment.m into the directory where your
%output...ECM.mat files are

clear all
for i = 0:1440
    filename1 = sprintf('output%08u_ECM.mat', i) ;
    DATA1 = read_ecm_data(filename1);
    PLOT = heatmap(DATA1.X, DATA1.Y, DATA1.data{1,1});
    title({'ECM Anisotropy' ; }) ;
    filename3 = sprintf('output%08u.png',i) ;
    saveas(PLOT,filename3) ;
end
    