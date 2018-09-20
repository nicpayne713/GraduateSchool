% Nicholas Payne
% Math 561 HW 5 #2
%
%lose all
clc
clear all
load co2.dat
t = (1:216)';
p = polyfit(t,co2,1);
co2_1 = co2 - polyval(p,t);

% PART a

figure
plot(t,co2, 'x-') 
title('Plot of data containing linear progression')
xlabel('Data points taken once per month for 216 months')
ylabel('Co2 levels measured at each month')

figure
plot(t, co2_1, 'x-')
title('Plot of data eliminating the linear progression')
xlabel('Data points taken once per month')
ylabel('Change in Co2 levels at points of measurement')

Y = fft(co2_1); % computes the DFT using FFT algorithm
Y1 = abs(Y(1:108)); % creates vector with absolute values
%                   of the first half of the entires in the FFT vector
t1 = t(1:108); % time vector to plot Y1

figure
plot(t1, Y1, 'r-')
title('Plot of absolute value of DFT compared to the data with which the FFT was performed on')
xlabel('Measurement points')
ylabel('FFT outputs')
% PART b
%
% These for loops just tell us at which data points in the FFT the peaks
% occur - I found the peaks by looking at the vector and the loops just
% tell us where the entires occur

for j = 1:50;
    if abs(Y1(j) - 274.9881) < 1;
        w = j;
    end
end
for k = 1:50;
    if abs(Y1(k) - 74.6986) < 1;
        r = k;
    end
end
% these are the entires at which the peaks occur in the FFT
%First peak = 2 
w; % = 19
r; % = 37

