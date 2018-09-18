%Move this and read_microenvironment.m into the directory where your
%output...ECM.mat files are

clear all;
for i = 0:144
    filename1 = sprintf('output%08u_ECM.mat', i) ;
   % filename2 = sprintf('output%08u.xml',i) ;
    ECM = read_ecm_data(filename1);
   % DATA2 = read_MultiCellDS_xml(filename2);
    
    fig = figure;
    set(fig, 'Visible', 'off');
    scaled_X_fiber = ECM.data{1,1} .* ECM.data{1,3};
    scaled_Y_fiber = ECM.data{1,1} .* ECM.data{1,4};
    quiver(ECM.X, ECM.Y, scaled_X_fiber, scaled_Y_fiber);

    axis ([-1000 1000 -1000 1000]);
    time = 0;
   %string(sprintf('%d',DATA2.metadata.current_time/60)) ;
    % title({'ECM Fiber Alignment' ; 'Time = ' + time + ' hrs'}) 
    print(sprintf('output%08u.svg',i), '-dsvg', '-r300') ;
end