clear all
for i = 0:1440
    filename1 = sprintf('./Output/output%08u_ECM.mat', i) ;
    DATA1 = read_ecm_data(filename1);
    PLOT = contourf(DATA1.X, -DATA1.Y, DATA1.data{1,2});
    hold on
    title({'ECM Density' ; }) ;
    MCDS = read_MultiCellDS_xml(sprintf('output%08u.xml', i), './Output');
    plot(MCDS.discrete_cells.state.position(:,2), -MCDS.discrete_cells.state.position(:,1) ,'bo' ); 
    axis image;
    print(sprintf('./Output/output%08u.png',i), '-dpng', '-r0') ;
end
    