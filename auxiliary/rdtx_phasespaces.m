% This function generates phasespaces from all particles at the start and
% end of the simulation. fields are phasespaces desired, e.g. {'x0p0',
% 'x2x1'}. limits are an {} array with same size as fields of an array of
% [x0min,x0max,xmin,xmax,Nx,y0min,y0max,ymin,ymax,Ny] where x and y are the two spaces
% 0 x0,y0 is at start and x,y is end
% 'default' will give limits at + and - 2* actual limits 
%
% rdtx_phasespaces(dir,fields,limits)
%

function rdtx_phasespaces(dir,fields,limits)

% otherwise surf crashes!
%opengl software
fontname='times';
fontsize=16;
[x0,v0,x,v,Erad] = rdtx_openall(dir);

Nfields = length(fields);
for ii=1:Nfields
    
    [X0,Y0]=rdtx_phase2d(x0,v0,fields{ii});
    [X,Y]=rdtx_phase2d(x,v,fields{ii});
    llimits=limits{ii};
    
    if strcmpi(llimits,'default')
         llimits0=rdtx_defaultlimits(X0,Y0);
         llimits=rdtx_defaultlimits(X,Y);
    else     
        llimits0=[llimits(1) llimits(2) llimits(5) llimits(6) llimits(7) llimits(10) ];
        llimits=[llimits(3) llimits(4) llimits(5) llimits(8) llimits(9) llimits(10) ];
    end
    gridx0{ii}=linspace(llimits0(1),llimits0(2),llimits0(3));
    gridy0{ii}=linspace(llimits0(4),llimits0(5),llimits0(6));
    gridd0{ii}=rdtx_proj(X0,Y0,gridx0{ii},gridy0{ii});
 
    gridx{ii}=linspace(llimits(1),llimits(2),llimits(3));
    gridy{ii}=linspace(llimits(4),llimits(5),llimits(6));
    gridd{ii}=rdtx_proj(X,Y,gridx{ii},gridy{ii});
end

    

for ii=1:Nfields
    maxgrdd=max(max(max(gridd{ii})), max(max(gridd0{ii})));
    subplot(2,Nfields,ii);
    imagesc(gridx0{ii},gridy0{ii},gridd0{ii}); axis xy  
   
  hold on;
   surfc(gridx0{ii},gridy0{ii},gridd0{ii}); axis xy  
 set(gca,'fontname',fontname,'fontsize',fontsize);
     name=fields{ii};
    xlabel([name(1) '_' name(2) rdtx_units(name(1))],'fontname',fontname,'fontsize',fontsize);
    ylabel([name(3) '_' name(4) rdtx_units(name(1))],'fontname',fontname,'fontsize',fontsize);
    %title([name(1) '_' name(2) name(3) '_' name(4) ' initial']);
    axis tight square;
    c=colorbar;  
    initpos = get(c,'Position');
set(c, ...
   'Position',[initpos(1)*1.05+initpos(3)*0.25 initpos(2)+initpos(4)*0.25 ...
      initpos(3)*0.5 initpos(4)*0.5])
     set(c,'fontname',fontname,'fontsize',fontsize);
    caxis([0 maxgrdd]);
    
    myname=['$\frac{d^2N}{' 'd' name(1) '_' name(2) 'd' name(3) '_' name(4) '}$'];
	title(c,myname,'interpreter','latex','fontname',fontname,'fontsize',fontsize);

    subplot(2,Nfields,ii+Nfields);
    imagesc(gridx{ii},gridy{ii},gridd{ii}); axis xy
     hold on;
   surfc(gridx{ii},gridy{ii},gridd{ii}); axis xy  
   set(gca,'fontname',fontname,'fontsize',fontsize);
    xlabel([name(1) '_' name(2) rdtx_units(name(1))],'fontname',fontname,'fontsize',fontsize);
    ylabel([name(3) '_' name(4) rdtx_units(name(1))],'fontname',fontname,'fontsize',fontsize);
    %title([name(1) '_' name(2) name(3) '_' name(4) ' final']);
   axis tight square;
   c=colorbar;%('NorthOutside'); 
   initpos = get(c,'Position');
set(c, ...
   'Position',[initpos(1)*1.05+initpos(3)*0.25 initpos(2)+initpos(4)*0.25 ...
      initpos(3)*0.5 initpos(4)*0.5])
     set(c,'fontname',fontname,'fontsize',fontsize);
    caxis([0 maxgrdd]);
    title(c,myname,'interpreter','latex','fontname',fontname,'fontsize',fontsize);


end
    load('rdtxcolormap','mycmap');
    set(gcf,'Colormap',mycmap);
clear