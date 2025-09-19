%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method 1
x = 0:.1:20;
y1 = 5*sin(x) + 2*x - x.^2 +.3*x.^3 - .2*(x-15).^4 - 10*x.^2.*cos(x./3+12).^3 + .5*(x-12).^4;

% add lots of noise:
r = randi(1000,1,201) - 500;
y2 = y1+r;

%Now make a 1D Gaussian filter, normalize it and convolve it with our function:
g = gausswin(20); % <-- this value determines the width of the smoothing window
g = g/sum(g);
y3 = conv(y2, g, 'same')

figure;
hold on; 
plot(y1, 'r', 'linewidth', 3); 
plot(y2, 'b'); 
plot(y3, 'g', 'linewidth', 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%method 2
x = 1:100;
A = cos(2*pi*0.05*x+2*pi*rand) + 0.5*randn(1,100);
[B, window] = smoothdata(A,'gaussian');
window    % window = 4

% Smooth the original data with a larger window of length 20. Plot the smoothed data for both window lengths.
C = smoothdata(A,'gaussian',20);
plot(x,B,'-o',x,C,'-x')
legend('Small Window','Large Window')