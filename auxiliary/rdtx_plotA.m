% This function plots the components of the four potential in 4 windows
%
% function rdtx_plotA(dir,number)



function rdtx_plotA(dir,number)

[Ax,Ay,Az,phi,zgrid,xgrid,time]=rdtx_openA(dir,number);

subplot(2,2,1); imagesc(zgrid,xgrid,Ax); title(['A_x at t=' num2str(time)]);  colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
subplot(2,2,2); imagesc(zgrid,xgrid,Ay); title(['A_y at t=' num2str(time)]);  colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
subplot(2,2,3); imagesc(zgrid,xgrid,Az); title(['A_z at t=' num2str(time)]); colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
subplot(2,2,4); imagesc(zgrid,xgrid,phi); title(['\phi at t=' num2str(time)]); colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight

load('rdtxcolormap','cmap')
set(gcf,'Colormap',cmap);