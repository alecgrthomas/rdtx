%
% 
% This function opens particle track data, returning fourvectors for 
% spacetime x, four velocity v, four potential seen by particle a, average
% power emitted by particle Pave, the QED paramter chi and the name of
% particle (which is just `Par_num_N'). The inputs are the root directory
% and the range of particles as a vector e.g. 0:2:20. Data are returned as
% an array of cells each containing an Nt x 4 matrix.
%
% [x,v,a,Pave,chi,name,damping] = rdtx_openpar(directory,par_range)
%

function [x,v,a,Pave,chi,S,name] = rdtx_openpar(directory,par_range)  
ii=1;
for num = par_range
filename=[directory '/output_par_' int2str(num) '.dat'];
fid = fopen(filename);
tempdata = textscan(fid, '%s', 1);
while strcmp(tempdata{:},'0')==0
tempdata = textscan(fid, '%s', 1);
end
data=[];
name = textscan(fid, '%s', 1);
name = name{1};
coords = textscan(fid, '%s', 1);
coords=coords{1};
nx=textscan(fid, '%d', 1);
nx = nx{1};
ny = 18; % added chi+S
vphase = textscan(fid, '%n', 1);
vphase=vphase{1};
textscan(fid, '%d', 1);
for j=1:ny  
    textscan(fid, '%s', 1);
end

for i=1:nx
    tempdata=textscan(fid, '%n', ny);
    if tempdata{1}~=1.2345
    data(i,:)=tempdata{1};
    end
end
fclose(fid);

x{ii}=data(:,1:4); v{ii}=data(:,5:8); a{ii}=data(:,10:13); 
Pave{ii}=data(:,9); chi{ii} = data(:,14); S{ii}=data(:,15:18);
ii=ii+1;
end

clear tempdata ans fid filename i j data