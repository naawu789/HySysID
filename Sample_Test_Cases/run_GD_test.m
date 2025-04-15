clc; clear all; close all;

t = 0.01:0.01:40;
%theta = 2; phi = 2*sin(t);
%theta = 100*rand(2,2); phi = [rand(1,length(t)); rand(1,length(t))];
%theta = 100*rand(2,3); phi = [rand(1,length(t)); rand(1,length(t)); rand(1,length(t))];
theta = [0 1 0; 0 0 -9.8]; phi = [abs(sin(t*pi/2)); -mod(t,2)+1; ones(1,length(t))];

y = theta*phi;
figure(1);
subplot(2,1,1); hold on
for i = 1:size(y,1)
    plot(t,y(i,:))
end
ylabel('y');
subplot(2,1,2); hold on
for i = 1:size(phi,1)
    plot(t,phi(i,:))
end
ylabel('$\phi$','interpreter','latex'); xlabel('t');
%}

gamma = 1/10;
hat_theta = zeros(size(theta));
error = zeros(length(t),1);

for i = 1:length(t)
    hat_theta = hat_theta - gamma*(  phi(:,i)* ( hat_theta*phi(:,i) - y(:,i) )'  )';

    error(i) = norm(hat_theta - theta);
end
hat_theta - theta

figure(2);
plot(t,error)
xlabel('t'); ylabel('$e$','interpreter','latex');
