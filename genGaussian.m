function [acc ,vel, pos]=genGaussian(mag,sigma,duration)
%this is from Jing and Jian - how they generate the actual trajectories
%Updated by Adam. (+ updated quickly at the bottom take to into account the constant needed to multiply e.g. acc X div to get good magnitude)
offset=0.0;
tdist = mag/2;
xval = sqrt(2.0) * (duration - duration/2.0)/(duration/sigma);
namp = sqrt(2.0) / (sqrt(pi)*duration/2.0/sigma*erf(xval));
amplitude = tdist*namp;
area_half = -amplitude/2.0 * sqrt(pi) * sqrt(2.0) * duration/2.0/sigma * erf((-duration/2.0) * sqrt(2.0) / (duration/sigma));
%div=60;
div=1000; %AZ can't change since used in in POST_monkey_regression
trajectory_length = floor(duration*div+0.5) + 1;
pcalc = sqrt(pi)*duration/2.0/sigma;

g0=zeros(trajectory_length,1);
g1=zeros(trajectory_length,1);
g2=zeros(trajectory_length,1);

for i=1:trajectory_length
    timeStep = (i-1) / div;
   
    xval = sqrt(2.0) * (timeStep - duration/2.0)/(duration/sigma);
	namp = sqrt(2.0) / (pcalc*erf(xval));
	ti = amplitude/namp;
    
    g0(i)=ti + area_half + offset;
    
    g1(i)=amplitude*exp(-xval^2.0);
end

g0=g0/max(g0)*mag;
g2=diff(g0)*div;
g3=diff(g2)*div/9.8;
velmax=max(g2);

% plot(g1,'-b'); hold on;
% plot(g2,'-r');
% %plot(g0,'-g');
% plot(g3,'-g');

%by AZ
acc=diff(g1);
vel=cumtrapz(acc);
pos=cumtrapz(vel);

%normalize (added quickly...)
acc=acc*div;
vel=vel;
pos = pos/div;
% pos = abs(mag)*pos(1:end)/max(abs(pos));


