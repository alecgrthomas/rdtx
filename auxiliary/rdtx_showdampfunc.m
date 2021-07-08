% Function that shows the damping function
% returns tau and damping
% function [tau,damping]=rdtx_showdampfunc(directory)

function [tau,damping]=rdtx_showdampfunc(directory)

filename=[directory '/dampinginfo.dat'];
fid = fopen(filename);
tempdata=textscan(fid, '%n', 1);
        nt=tempdata{1};
        tempdamp=[];
for i=1:nt
        tempdata=textscan(fid, '%n', 2);
        tempdamp=[tempdamp; tempdata{1}'];   
end

tau=tempdamp(:,1);
damping=tempdamp(:,2);

figure;
plot(tau,damping); xlabel ('\omega_0\tau'); ylabel('damping factor'); title('Spectrum damping factor');