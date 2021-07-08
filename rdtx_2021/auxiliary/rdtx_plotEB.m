% This function plots the components of the electromagnetic field in 6 windows
%
% function rdtx_plotEB(dir,number)


function rdtx_plotEB(dir,number)
[Ex,Ey,Ez,Bx,By,Bz,zgrid,xgrid,time] = rdtx_openEB(dir,number);   
    
subplot(2,3,1); imagesc(zgrid,xgrid,Ex); title(['E_x at t=' num2str(time)]);  colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
subplot(2,3,2); imagesc(zgrid,xgrid,Ey); title(['E_y at t=' num2str(time)]);  colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
subplot(2,3,3); imagesc(zgrid,xgrid,Ez); title(['E_z at t=' num2str(time)]); colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
subplot(2,3,4); imagesc(zgrid,xgrid,Bx); title(['B_x at t=' num2str(time)]);  colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
subplot(2,3,5); imagesc(zgrid,xgrid,By); title(['B_y at t=' num2str(time)]);  colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
subplot(2,3,6); imagesc(zgrid,xgrid,Bz); title(['B_z at t=' num2str(time)]); colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
load('rdtxcolormap','cmap')
set(gcf,'Colormap',cmap);