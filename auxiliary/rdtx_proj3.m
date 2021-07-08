% rdtx analysis 2011
% phasespace projection from scattered data on 3D volume
% on x against y
% gridx,y,z define the grid -> 1D grids
% function grid=rdtx_proj3(x,y,gridx,gridy,gridz)
function grid=rdtx_proj3(x,y,gridx,gridy,gridz)

dx=gridx(2)-gridx(1);
dy=gridy(2)-gridy(1);
dz=gridz(2)-gridz(1);

xmin=min(gridx);
ymin=min(gridy);
zmin=min(gridz);
%get rid of offset

normx=(x-xmin)/dx;
normy=(y-ymin)/dy;
normgridx=(gridx-xmin)/dx;
normgridy=(gridy-ymin)/dy;

Nx=max(size(gridx));
Ny=max(size(gridy));
Npar=size(x);

if (size(y)~=Npar)
    warning('x and y grids not same size');
end

grid=zeros(Ny,Nx); % rows, columns!

for ii=1:Npar
    xpos=int16(normx(ii))+1;
    ypos=int16(normy(ii))+1;
 
    if (xpos>0) && (xpos<=Nx) && (ypos>0) && (ypos<=Ny)
        grid(xpos,ypos)=grid(xpos,ypos)+1;
    end
end