% function to open 2D potential data from rdtx. directory dir and
% filenumber number are input, and return is the four potential 
% components (phi,A), the 2 spatial grid vectors z,x and the time 
% corresponding to the file
%
% function [Ax,Ay,Az,phi,zgrid,xgrid,time] = rdtx_openA(dir,number) 

function [Ax,Ay,Az,phi,zgrid,xgrid,time] = rdtx_openA(dir,number) 
filename = [dir '/fields/potentials_' int2str(number) '.dat'];
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
        phi(i,:)=tempdata{1};
        tempdata=textscan(fid, '%n', nz);
        Ax(i,:)=tempdata{1};
        tempdata=textscan(fid, '%n', nz);
        Ay(i,:)=tempdata{1};
        tempdata=textscan(fid, '%n', nz);
        Az(i,:)=tempdata{1};

end
fclose(fid);

clear tempdata ans fid filename i j x y Ex Ey 