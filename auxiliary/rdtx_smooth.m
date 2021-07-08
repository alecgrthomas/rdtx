% This function performs a <numpoint> point smooth of a function f
% returns fsmooth
% function[fsmooth]= rdtx_smooth(f,numpoints) 
function[fsmooth]= rdtx_smooth(f,numpoints) % Size of the averaging window

window = ones(numpoints,numpoints)/numpoints/numpoints; 
fsmooth = convn(f,window,'same');