% Function that plots electron spectra before and after
% interaction. (non relativistic version )
% 
% function rdtx_espec(dir)

function rdtx_parspec(dir,rest_mass_energy)

[x0,v0,x,v,Erad] = rdtx_openall(dir);

%Nsteps = ceil(length(v0(:,1)));
Nsteps=200;
E0 = 0.5*(v0(:,2).^2+v0(:,3).^2+v0(:,4).^2);
Ef = 0.5*(v(:,2).^2+v(:,3).^2+v(:,4).^2);

Emax=(max(max([E0 Ef])))*1.1; %add 10%
E = linspace(0,Emax,Nsteps);
dE = E(2)-E(1);
dNdE = histc(E0,E)/dE;

dNdEf = histc(Ef,E)/dE;

subplot(2,1,1);
semilogy(E*rest_mass_energy,dNdE);
xlabel('E /keV');
ylabel('dN/dE');
title('Initial Energy Spectrum');

subplot(2,1,2);
semilogy(E*rest_mass_energy,dNdEf);
xlabel('E / keV');
ylabel('dN/dE');
title('Final Energy Spectrum');