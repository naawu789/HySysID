clc; clear all; close all;

t = 0.01:0.01:40;
x1 = @(t) ones(1,t); x2 = @(t) t; x3 = @(t) t.^2; x4 = @(t) t.^3;
%a1 = 1; a2 = 2; a3 = 3; a4 = 4;
%y = a1*x1(t) + a2*x2(t) + a3*x3(t) + a4*x4(t);
%plot(t,y)

theta = [1 2 3; 4 5 6];
phi = [x2(t); x3(t); x4(t)];
y = theta*phi;

gamma = 1/10;
hat_theta = zeros(size(theta));
error = zeros(length(t),1);

for i = 1:length(t)
    hat_theta = hat_theta - gamma*(  phi(:,i)* ( hat_theta*phi(:,i) - y(:,i) )'  )';

    error(i) = norm(hat_theta - theta);
end
hat_theta - theta

figure(1);
plot(t,y);
figure(2);
xlabel('t'); ylabel('$e$','interpreter','latex');
