% rdtx analysis 2011
% phasespace projection from scattered data
% x against y
% gridx,y define the grid -> 1D grids
function grid=rdtx_proj(x,y,gridx,gridy)

dx=abs(gridx(2)-gridx(1));
dy=abs(gridy(2)-gridy(1));

xmin=min(gridx);
ymin=min(gridy);

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
        grid(ypos,xpos)=grid(ypos,xpos)+1;
    end
end

grid = grid/dx/dy;