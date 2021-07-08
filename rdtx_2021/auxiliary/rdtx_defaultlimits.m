% define default limits for phase space
% based on input grids X0 and Y0
%
% function llimits=rdtx_defaultlimits(X0,Y0)
function llimits=rdtx_defaultlimits(X0,Y0)

defaultN=30; % default grid size
defaultfac=1; % default max/min rage factor (Xcenter+-Xrange*defaultfac)

centerx=(max(X0)+min(X0))*0.5;
rangex=(max(X0)-min(X0));
centery=(max(Y0)+min(Y0))*0.5;
rangey=(max(Y0)-min(Y0));

%fix zero limits
if rangex==0
    rangex = 1;
end
if rangey==0
    rangey = 1;
end
llimits=[centerx-rangex*defaultfac,centerx+rangex*defaultfac,defaultN,centery-rangey*defaultfac,centery+rangey*defaultfac,defaultN];

