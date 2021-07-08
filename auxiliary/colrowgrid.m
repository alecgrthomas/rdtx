%
% This function returns mesh similar to meshgrid in 2D but with the number
% of gridpoints in the direction x the same as in the vector
%
% [ X,Y ] = colrowgrid( x,y )
%
function [ X,Y ] = colrowgrid( x,y )

dummyx=ones(size(x));
dummyy=ones(size(y));
X=dummyy'*x;
Y=y'*dummyx;
end

