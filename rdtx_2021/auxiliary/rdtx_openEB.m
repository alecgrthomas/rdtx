% function to open 2D field data from rdtx. directory dir and
% filenumber number are input, and return is the 3 electric and 
% 3 magnetic field components, the 2 spatial grid vectors z,x and 
% the time corresponding to the file number
%
% function [Ex,Ey,Ez,Bx,By,Bz,zgrid,xgrid,time] = rdtx_openEB(dir,number) 

function [Ex,Ey,Ez,Bx,By,Bz,zgrid,xgrid,time] = rdtx_openEB(dir,number) 
filename = [dir '/fields/fields_' int2str(number) '.dat'];
fid = fopen(filename);
tempdata = textscan(fid, '%s', 1);
while strcmp(tempdata{:},'0')==0;
tempdata = textscan(fid, '%s', 1);
end
name = textscan(fid, '%s', 1);
name = name{1};
coords = textscan(fid, '%s', 1);
coords=coords{1};
nx = textscan(fid, '%n', 1);
nx = nx{1};
nz = textscan(fid, '%n', 1);
nz = nz{1};
x = textscan(fid, '%n', 1);
x = x{1};
z = textscan(fid, '%n', 1);
z = z{1};
time = textscan(fid, '%n', 1);
time = time{1};
zgrid = linspace(-z*0.5,z*0.5,nz);
xgrid = linspace(-x*0.5,x*0.5,nx);
for i=1:nx
    
        tempdata=textscan(fid, '%n', nz);
        Ex(i,:)=tempdata{1};
        tempdata=textscan(fid, '%n', nz);
        Bx(i,:)=tempdata{1};
        tempdata=textscan(fid, '%n', nz);
        Ey(i,:)=tempdata{1};
        tempdata=textscan(fid, '%n', nz);
        By(i,:)=tempdata{1};
        tempdata=textscan(fid, '%n', nz);
        Ez(i,:)=tempdata{1};
        tempdata=textscan(fid, '%n', nz);
        Bz(i,:)=tempdata{1};

end
fclose(fid);

clear tempdata ans fid filename i j x y 