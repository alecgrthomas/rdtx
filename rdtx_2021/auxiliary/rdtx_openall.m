%
% This function opens x and v positions of all particles at the start and
% end of the simulation.
%
%  [x0,v0,x,v,Erad] = rdtx_openall(dir) 
%
function [x0,v0,x,v,Erad] = rdtx_openall(dir)  %maximum energy in MeV
filename=[dir '/output_all.dat'];
fid = fopen(filename);
tempdata = textscan(fid, '%s', 1);
while strcmp(tempdata{:},'0')==0
tempdata = textscan(fid, '%s', 1);
end
nx = textscan(fid, '%d', 1);
nx = nx{1};
ny = 19;

%for j=1:ny  
    textscan(fid, '%s', ny);
%end

for i=1:nx
        tempdata=textscan(fid, '%n', ny);
        data(i,:)=tempdata{1};     
end
fclose(fid);
x0=data(:,2:5); v0=data(:,6:9); x=data(:,10:13); v=data(:,14:17); Erad=data(:,18);