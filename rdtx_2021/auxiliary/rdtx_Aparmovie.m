% plots particle trajectories (indicated by par_range) overlaid on 
% potential data for various time slices (indicated by time_range)
%
% function rdtx_Aparmovie(dir,par_range,time_range)

function rdtx_Aparmovie(dir,par_range,time_range)


ncontours=8;


screen = get(0,'ScreenSize');
figure('Position',[0.1*screen(3) 0.1*screen(4) 0.7*screen(3) 0.7*screen(4)]);
[xc,vc,ac,Pavec,namec] = rdtx_openpar(dir,par_range);
%[Ax,Ay,Az,phi,xgrid,ygrid,time]=OpenFLD(dir,time_range(1));
[Ax,Ay,Az,phi,zgrid,xgrid,time]=rdtx_openA(dir,0);

[Z,X]=colrowgrid(zgrid,xgrid);
for ii = time_range
    
%clf reset;
[Ax,Ay,Az,phi,zgrid,xgrid,time]=rdtx_openA(dir,ii);   
x=xc{1};
tlim=(x(:,1)<time);
inttime=sum(tlim);
if (~inttime)
    inttime=1;
end


modA=(Ax.^2+Ay.^2+Az.^2);
imagesc(zgrid,xgrid,modA); title(['A and \phi at t=' num2str(time)]); colorbar;xlabel('z-v_{frame}t'); ylabel('x'); axis tight

hold on
%[C,h] =contour(Z,X,phi,ncontours,'k'); %title(['\phi at t=' num2str(time)]); colorbar;xlabel('z-v_{ph}t'); ylabel('x'); axis tight
%clabel(C,h);
caxis([min(min(modA)),max(max(modA))]);
iii=1;
for num=par_range
    x=xc{iii};
    v=vc{iii};
    a=ac{iii};
    mylinecolor = rand(1,3);
plot(x(1:inttime,4),x(1:inttime,2),'k','LineWidth',2);% xlabel('z'); ylabel('x');
scatter(x(inttime,4),x(inttime,2),'MarkerFaceColor','b','LineWidth',2);%xlabel('z'); ylabel('x');
iii=iii+1;
end
    pause(0.1);
end
