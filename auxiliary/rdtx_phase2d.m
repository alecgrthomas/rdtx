% this function returns 2 vectors depending on the phasespace selectroin
% for example 'x1p1' returns x(2), p(2) (1 offset due to indexing) 
%
% function [X,Y] = rdtx_phase2d(x,v,phasespace) 

function [X,Y] = rdtx_phase2d(x,v,phasespace) 

p1=phasespace(1);
p2=phasespace(3);
n1=str2num(phasespace(2));
n2=str2num(phasespace(4));

switch p1
    case 'x'
        X=x(:,n1+1);
    case 'p'
        X=v(:,n1+1);
end

switch p2
    case 'x'
        Y=x(:,n2+1);
    case 'p'
        Y=v(:,n2+1);
end
