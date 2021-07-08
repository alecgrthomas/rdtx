% function that plots various particle trajectories of 2 coordinates
% par range is a vector inticated the particle labels to be displayed
% more than 20 is not recommended
%
% function  rdtx_plotpar(directory,par_range)  

function  rdtx_plotpar(directory,par_range)  

 [xc,vc,ac,Pavec,chic,Sc,namec] = rdtx_openpar(directory,par_range);  

%figure;
ii=1;

for num=par_range
    x=xc{ii};
    v=vc{ii};
    a=ac{ii};
    S=Sc{ii};

    mylinecolor = rand(1,3);
subplot(2,4,1); plot(x(:,1),v(:,1),'Color',mylinecolor); title('t-\gamma'); xlabel('t'); ylabel('\gamma'); 
hold on;
subplot(2,4,2); plot(v(:,4),v(:,2),'Color',mylinecolor); title('vz-vx');xlabel('vz'); ylabel('vx'); 
hold on;
subplot(2,4,3); plot(v(:,2),v(:,3),'Color',mylinecolor); title('v_x-v_y');xlabel('v_x'); ylabel('v_y'); 
hold on;
subplot(2,4,4); plot(x(:,1),x(:,2),'Color',mylinecolor); title('t-x'); xlabel('t'); ylabel('x'); 
hold on;
subplot(2,4,5); 
plot(x(:,4),x(:,2),'Color',mylinecolor); title('z-x, z-y');xlabel('z'); ylabel('x'); 
hold on;
plot(x(:,4),x(:,3),'Color',(1-mylinecolor.^2));  ylabel('x (y)');
legend('x','y');
%subplot(2,3,5); plot(x(:,4),x(:,3),'Color',mylinecolor); title('z-y');xlabel('z'); ylabel('y'); 
%hold on;
subplot(2,4,6); plot(x(:,2),x(:,3),'Color',mylinecolor); title('x-y');xlabel('x'); ylabel('y'); 
hold on;

subplot(2,4,7); plot(x(:,4),S(:,2),'Color',mylinecolor); xlabel('z'); ylabel('S'); title('S');
hold on;
plot(x(:,4),S(:,3));
plot(x(:,4),S(:,4));
Smag = sqrt(S(:,2).^2+S(:,3).^2+S(:,4).^2);
plot(x(:,4),Smag);
legend('S_x','S_y','S_z','S');


% Rest frame spin
Sx = S(:,2);
Sy = S(:,3);
Sz = S(:,4);
betax = v(:,2)./v(:,1);
betay = v(:,3)./v(:,1);
betaz = v(:,4)./v(:,1);
gamma = v(:,1);
SdotB = betax.*Sx+betay.*Sy+betaz.*Sz;

sx = Sx - gamma./(gamma+1.0).*SdotB.*betax;
sy = Sy - gamma./(gamma+1.0).*SdotB.*betay;
sz = Sz - gamma./(gamma+1.0).*SdotB.*betaz;


subplot(2,4,8); plot(x(:,4),sx,'Color',mylinecolor); xlabel('z'); ylabel('S'); title('Rest frame spin components');
hold on;
plot(x(:,4),sy);
plot(x(:,4),sz);
s = sqrt(sx.^2+sy.^2+sz.^2);
plot(x(:,4),s);
legend('s_x','s_y','s_z','s');
ii=ii+1;
end
