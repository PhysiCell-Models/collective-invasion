%Move this and read_microenvironment.m into the directory where your
%output...ECM.mat files are

clear all
numOutput = 1152;
images = cell(numOutput, 1);

for i = 0:numOutput
    filename1 = sprintf('output%08u_ECM.mat', i) ;
    ECM = read_ecm_data(filename1);
%     filename2 = sprintf('output%08u.xml',i) ;
%     DATA2 = read_MultiCellDS_xml(filename2);
    % LINEWIDTH
%     scaled_X_fibers = DATA1.data{1,1}*DATA1.data{1,3};
%     scaled_Y_fibers = DATA1.data{1,1}*DATA1.data{1,4};
%     PLOT = quiver(DATA1.X, DATA1.Y,scaled_X_fibers , scaled_Y_fibers);
    
    scaled_X_fiber = ECM.data{1,1} .* ECM.data{1,3};
    scaled_Y_fiber = ECM.data{1,1} .* ECM.data{1,4};

    PLOT = figure(1);
    set(PLOT, 'Visible', 'off');
    title('ECM Fiber Alignment') ;
    
    
    quiver(ECM.X, ECM.Y, scaled_X_fiber, scaled_Y_fiber)

    axis image
%     time = 0;%string(sprintf('%d',DATA2.metadata.current_time/60)) ;
    
    filename3 = sprintf('fiber_alignment%08u.png',i) ;
    saveas(PLOT,filename3) ;
    
    images{i + 1} = imread(filename3);
end

% create the video writer with 1 fps
 writerObj = VideoWriter('fiber_alignment_quiver.mp4');
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
