% function that plots various particle trajectories of 3 coordinates
% par range is a vector inticated the particle labels to be displayed
% more than 20 is not recommended
%
% function  rdtx_plotpar3(directory,par_range)  

function  rdtx_plotpar3(directory,par_range)  

 [xc,vc,ac,Pavec,chic,namec] = rdtx_openpar(directory,par_range);  
 
figure;
ii=1;
for num=par_range
    x=xc{ii};
    v=vc{ii};
    a=ac{ii};
    mylinecolor = rand(1,3);
subplot(1,3,1); plot3(v(:,2),v(:,3),v(:,4),'Color',mylinecolor); xlabel('p_x/mc'); ylabel('p_y/mc');  zlabel('p_z/mc'); 
hold on;
subplot(1,3,2); plot3(x(:,2),x(:,3),x(:,4),'Color',mylinecolor); xlabel('x\omega_p/c'); ylabel('y\omega_p/c');  zlabel('z\omega_p/c'); 
hold on;
subplot(1,3,3); plot3(a(:,2),a(:,3),a(:,4),'Color',mylinecolor); xlabel('eA_x/mc'); ylabel('eA_y/mc');  zlabel('eA_z/mc'); 
hold on;
ii=ii+1;
end
% 
% figure;
% subplot(2,3,1); plot(x(:,1),v(:,1)); title('t-\gamma'); xlabel('t'); ylabel('\gamma'); 
% subplot(2,3,2); plot(v(:,4),v(:,2)); title('vz-vx');xlabel('vz'); ylabel('vx'); 
% subplot(2,3,3); plot(x(:,4)-vphase*x(:,1),x(:,2)); title('(z-v_{ph}t)-x');xlabel('z-v_{ph}t'); ylabel('x'); 
% subplot(2,3,4); plot(x(:,4),x(:,2)); title('z-x');xlabel('z'); ylabel('x'); 
% subplot(2,3,5); plot(x(:,4),x(:,3)); title('z-y');xlabel('z'); ylabel('y'); 
% subplot(2,3,6); plot(x(:,2),x(:,3)); title('x-y');xlabel('x'); ylabel('y'); 
% 
% figure;
% subplot(1,2,1); plot(x(:,1),-(max(v(:,1))-v(:,1))/max(v(:,1))*100); title('t-{\Delta}E/E'); xlabel('t'); ylabel('{\Delta}E/E (%)'); 
% %subplot(1,2,2); plot(x(:,1),v(:,2)); title('t-v_x'); xlabel('t'); ylabel('v_x'); 
% 
% vdot = v(:,2);
% dt=(x(:,1) - circshift(x(:,1),1));
% vdot = (vdot - circshift(vdot,1))./dt;
% 
% %
% asq = (a(:,1).^2+a(:,2).^2+a(:,3).^2+a(:,4).^2);
% drag = -4*v(:,1).^2.*asq;
% subplot(1,2,2); plot(x(:,1),drag); title('t-drag'); xlabel('t');  ylabel('drag');
% 
% deltagamma = (v(1,1)- 6.4e-24*2.36e15*sum(drag.*dt))/v(1,1)-1;
%  
%  %deltagamma*100
%  
%  %figure; plot(x(:,1),asq./max(asq),x(:,1),exp(-2*((x(:,1)-65).^2/(65/2)^2))); xlabel('t');  ylabel('');
%  
%  figure; subplot(1,2,1); plot(x(:,1),chi); title('t-\chi'); xlabel('t');  ylabel('\chi');
%  subplot(1,2,2); plot(x(:,1),(3.7*chi.^3+31.0*chi.^2+12.0*chi+1.0).^(-4.0/9.0)); title('t-g(\chi)'); xlabel('t');  ylabel('g(\chi)');