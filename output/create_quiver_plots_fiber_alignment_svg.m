%IMPORTANT NOTE: Before running this script for a data set larger than
%about 200 of our current time steps, go in Matlab to Home, then
%preferences, then general, then Java Heap Memory, and increase it
%I'm not sure where it needs to be set exactly, but I made mine to
%~1000 MB


%Move this and read_microenvironment.m into the directory where your
%output...ECM.mat files are

clear all;
for i = 0:1440
    filename1 = sprintf('output%08u_ECM.mat', i) ;
   % filename2 = sprintf('output%08u.xml',i) ;
    ECM = read_ecm_data(filename1);
   % DATA2 = read_MultiCellDS_xml(filename2);
    
    fig = figure;
    set(fig, 'Visible', 'off');
    scaled_X_fiber = ECM.data{1,1} .* ECM.data{1,3};
    scaled_Y_fiber = ECM.data{1,1} .* ECM.data{1,4};
    quiver(ECM.X, ECM.Y, scaled_X_fiber', scaled_Y_fiber');

    axis image;
    time = 0;
   %string(sprintf('%d',DATA2.metadata.current_time/60)) ;
    % title({'ECM Fiber Alignment' ; 'Time = ' + time + ' hrs'}) 
    %using -painters argument makes sure we're outputting an actual 
    %vector file type and not just a bitmap image embedded in an svg
    print(sprintf('output%08u.svg',i), '-dsvg', '-r300', '-painters') ;
end