% This function plots the components of the electromagnetic fields in 6 windows
% and steps through ith a slight pause to make a movie. 
% User specifies range, which is a vector of the form start:step:end
%
% function rdtx_EBmovie(dir,range)

function rdtx_EBmovie(dir,range)
screen = get(0,'ScreenSize');
figure('Position',[0.1*screen(3) 0.1*screen(4) 0.7*screen(3) 0.7*screen(4)]);
for ii = range
  rdtx_plotEB(dir,ii);
  pause(0.1);
end
