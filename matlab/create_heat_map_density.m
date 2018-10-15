clear all
for i = 0:288
    filename1 = sprintf('output%08u_ECM.mat', i) ;
    DATA1 = read_ecm_data(filename1);
    PLOT = heatmap(DATA1.X, DATA1.Y, DATA1.data{1,2});
    title({'ECM Density' ; }) ;
    print(sprintf('output%08u.svg',i), '-dsvg', '-r300', '-painters') ;
end
    