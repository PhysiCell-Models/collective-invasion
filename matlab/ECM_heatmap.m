%Move this and read_microenvironment.m into the directory where your
%output...ECM.mat files are

clear all

rawColors = [1 0 1; 0 0 1; 0 1 1; 0 1 0; 1 1 1; 1 1 0; 1 0.5 0; 1 0 0; 1 0 1]; %Raw colors used (magenta, blue, cyan, green, white, yellow, orange, red, magenta)
depth = 10; %How many colors between two of the raw colors
customColorMap = []; %Matrix that holds the color map
for i = 1:size(rawColors, 1)-1
    customColorMap = [customColorMap; reshape(linspace(rawColors(i, 1), rawColors(i+1,1), depth), [depth, 1]) reshape(linspace(rawColors(i, 2), rawColors(i+1,2), depth), [depth, 1]) reshape(linspace(rawColors(i, 3), rawColors(i+1,3), depth), [depth, 1])];
end


for i = 0:0
    filename1 = sprintf('output%08u_ECM.mat', i) ;
    ECM = read_ecm_data(filename1);
%     filename2 = sprintf('output%08u.xml',i) ;
%     DATA2 = read_MultiCellDS_xml(filename2);
    % LINEWIDTH
%     scaled_X_fibers = DATA1.data{1,1}*DATA1.data{1,3};
%     scaled_Y_fibers = DATA1.data{1,1}*DATA1.data{1,4};
%     PLOT = quiver(DATA1.X, DATA1.Y,scaled_X_fibers , scaled_Y_fibers);
    %set(figure(1), 'Visible', 'off');
    %figure(1)
    scaled_X_fiber = ECM.data{1,1} .* ECM.data{1,3};
    %scaled_Y_fiber = ECM.data{1,1} .* ECM.data{1,4};

    %h = heatmap(ECM.X, ECM.Y, ECM.data{1,3});
    h = heatmap(ECM.X, ECM.Y, scaled_X_fiber);
    h.ColorMethod = 'none'; %Sets the coloring method to use the value of z (x or y component)
    h.ColorLimits = [-1,1]; %Sets the scale to be constant with min=-1 and max=1
    %h.Colormap = hsv; %This color scale wraps around such that min = [1 0 0] and max = [1 0 0].  
                      %This is good for looking at the fibers as a value near -1 is equivilent to a value of 1 (as the model treats fibers as lines not vectors).
    %h.Colormap = parula; %This is a more gradual color scale going from
    %[0 0 1] to [1 1 0] linearly
    
    h.Colormap = customColorMap;
    
    grid off
    %axis image
    time = 0;%string(sprintf('%d',DATA2.metadata.current_time/60)) ;
%     title({'ECM Fiber Alignment' ; 'Time = ' + time + ' hrs'}) ;
    filename3 = sprintf('output%08u.png',i) ;
    saveas(figure(1),filename3);
    images{i + 1} = imread(filename3);
end

% create the video writer with 1 fps
 writerObj = VideoWriter('heatmap.mp4');
 writerObj.FrameRate = 60;

 % open the video writer
 open(writerObj);
 % write the frames to the video
 for u=1:length(images)
     % convert the image to a frame
     frame = im2frame(images{u});
     
     writeVideo(writerObj, frame);
     
 end    
 % close the writer object
 close(writerObj);
