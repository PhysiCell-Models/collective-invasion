clear all
numOutput = 31;
images = cell(numOutput, 1);

for i = 0:numOutput
    filename1 = sprintf('output%08u_ECM.mat', i) ;
    filename2 = sprintf('output%08u.xml',i) ;
    DATA1 = read_ecm_data(filename1);
    DATA2 = read_MultiCellDS_xml(filename2);
    PLOT = figure ;
    set(PLOT, 'Visible', 'off');
%     PLOT = heatmap(DATA1.X, DATA1.Y, DATA1.data{1,1});
    title({'ECM Anisotropy' ; }) ;
    contourf(DATA1.X,DATA1.Y,DATA1.data{1,1}) ;
    time = string(sprintf('%d',DATA2.metadata.current_time/60)) ;
    title({'ECM' ; 'Time = ' + time + ' hrs'}) ;
    colorbar ;
   set(colorbar,'Ylim',[0,1])
   axis image ;
   hold on ;
   plot( DATA2.discrete_cells.state.position(:,1) , DATA2.discrete_cells.state.position(:,2) ,'ko' ) ;
   filename3 = sprintf('output%08u.png',i) ;
   saveas(PLOT,filename3) ;

   images{i + 1} = imread(filename3);
end

% create the video writer with 1 fps
 writerObj = VideoWriter('anistropyContour.mp4');
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
    