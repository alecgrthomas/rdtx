% This function displays a spectrum from base directory basedir. component 
% is either '' or 'para' or 'perp', for unseparated or separated (parallel 
% and perpendicular) polarizations for the radiation. range is the range of
% theta bins, e.g. 0:2:19. smoothing is the smoothing level, less than 1 
%is no smoothing 
%
% rdtx_plotspectrum(basedir,component,range,smoothing)

function rdtx_plotspectrum(basedir,component,range,smoothing)

[omega,theta,f] = rdtx_openspectrum(basedir,component,range) ;

if (smoothing >1)
    f=rdtx_smooth(f,smoothing);
end
[N,M] = size(f);

% need to check if logarithmic
if (omega(2) - omega(1)) == (omega(3) - omega(2)) % linear
    

if (N>1)

    imagesc(omega,theta,f);
    colorbar;
    load('rdtxcolormap','mycmap');
    set(gcf,'Colormap',mycmap);
    xlabel('\omega / \omega_0');
    ylabel('\theta / \rm{rad}');
else
    plot(omega,f);
    xlabel('\omega / \omega_0');
    ylabel('d^2I/{d\omega}d\Omega (eV s sr^{-1})'); axis tight;
end

else
    if (N>1)

    imagesc(log10(omega),theta,f);
    colorbar;
    load('rdtxcolormap','mycmap');
    set(gcf,'Colormap',mycmap);
    xlabel('log_{10}(\omega / \omega_0)');
    ylabel('\theta / \rm{rad}');
else
    semilogx(omega,f);
    xlabel('\omega / \omega_0');
    ylabel('d^2I/{d\omega}d\Omega (eV s sr^{-1})'); axis tight;
    end
end

end % of function