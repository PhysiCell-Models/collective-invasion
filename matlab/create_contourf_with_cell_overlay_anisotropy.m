clear all

for i = 0:1000
    
    filename1 = sprintf('./Output/output%08u_ECM.mat', i) ;
    DATA1 = read_ecm_data(filename1);
    PLOT = contourf(DATA1.X, -DATA1.Y, -DATA1.data{1,1});
    hold on
    title({'ECM Density with Cell Overlay' ; }) ;
    MCDS = read_MultiCellDS_xml(sprintf('output%08u.xml', i), './Output');
    plot(MCDS.discrete_cells.state.position(:,2), -MCDS.discrete_cells.state.position(:,1) ,'o', 'MarkerSize', 3, 'MarkerEdgeColor', [1,0,0],'MarkerFaceColor', [1,1,1] ); 
    axis image;
    print(sprintf('./Output/output%08u.png',i), '-dpng', '-r0') ;
end