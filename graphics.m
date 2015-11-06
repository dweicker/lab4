%This scripts plots the variables figures needed for part f) of LAB4
%We use : hx = 0.1
%         ht = 0.005
%         tend = 2
close all;
clear all;

[u,x,t] = tempEE(10,0.005,2);

figure;
surf(x(length(x):-1:1),t,u');title('Evolution of the temperature in the rod');xlabel('Length [m/L]');
ylabel('Time [s/tp]');zlabel('Temperature [°C/T0]');

figure;
T = [0.5 1 1.5 2];
ti = floor(T/0.005)+1;
for i = 1:4
    subplot(2,2,i);
    plot(x(length(x):-1:1),u(:,ti(i)));title(sprintf('Temperature at tau = %f',T(i)));
    xlabel('Length [m/L]'); ylabel('Temperature [°C/T0]');axis([0 1 0 1]);
end

%We can also plot an unstable solution
%We use : hx = 0.1
%         ht = 0.0051
%         tend = 2
[uUnstable,x,t] = tempEE(10,0.0051,2);

figure;
surf(x(length(x):-1:1),t,uUnstable');title('Unstable solution');xlabel('Length [m/L]');
ylabel('Time [s/tp]');zlabel('Temperature [°C/T0]');



