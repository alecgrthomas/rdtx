% This function opens a spectrum from base directory basedir. component is
% either '' or 'para' or 'perp', for unseparated or separated (parallel and
% perpendicular) polarizations for the radiation. range is the range of
% theta bins, e.g. 0:2:19. Output spectrum is f
%
% [omega,theta,f] = rdtx_openspectrum(basedir,component,range)
%
function [omega,theta,f] = rdtx_openspectrum(basedir,component,range)  

switch component
    case ''
        
    case 'para'
        component = 'para_';
        
    case 'perp'
        component = 'perp_';
        
    otherwise
        component = '';
end

iii=1;
for number=range
num = num2str(number);

filename=[basedir '/spectrum_' component 'theta_' num '.dat'];
fid = fopen(filename);
for ii=1:4
tempdata = textscan(fid, '%s', 1);
end
thetatemp = textscan(fid, '%n', 1);
thetatemp = thetatemp{1};
for ii=1:4
tempdata = textscan(fid, '%s', 1);
end
while strcmp(tempdata{:},'E')==0
tempdata = textscan(fid, '%s', 1);
end
nx = textscan(fid, '%d', 1);
nx = nx{1};
ny = 2;
textscan(fid, '%d', 1);
%for j=1:ny  
    textscan(fid, '%s', ny);
%end

for ii=1:nx
 %   for j=1:ny 
        tempdata=textscan(fid, '%n', ny);
        data(ii,:)=tempdata{1};        
 %   end
end
fclose(fid);
omega=data(:,1); f(:,iii)=data(:,2); theta(:,iii) = thetatemp;
iii=iii+1;
end
f=f';
end % of function
