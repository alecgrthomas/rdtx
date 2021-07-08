% Function that plots electron spectra before and after
% interaction. 
% 
% function rdtx_espec(dir)

function rdtx_parspec(dir)

[x0,v0,x,v,Erad] = rdtx_openall(dir);

%Nsteps = ceil(length(v0(:,1)));
Nsteps=2000;
E0 = v0(:,1);
Ef = v(:,1);

Emax=(max(max([E0 Ef]))-1.0)*1.1; %add 10%
E = linspace(0,Emax,Nsteps);
dE = E(2)-E(1);
dNdE = histc(E0-1.0,E)/dE;

dNdEf = histc(Ef-1.0,E)/dE;

subplot(2,1,1);
plot(E,dNdE);
xlabel('\gamma - 1');
ylabel('dN/dE');
title('Initial Energy Spectrum');

subplot(2,1,2);
plot(E,dNdEf);
xlabel('\gamma - 1');
ylabel('dN/dE');
title('Final Energy Spectrum');