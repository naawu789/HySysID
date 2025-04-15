clc; clear all; close all;

global AB HJ wtrue btrue r exnoise tx jx

% Bouncing ball example
AB = [0 1 0; 0 0 -9.81];
HJ = [0 0 0; 0 -0.8 0];
wtrue = [1 0];
btrue = [0];

% Numerical example
AB = [-0.7 3 1; -2 -0.5 0];
HJ = [-1 0.5 0; -0.5 -1 1];
wtrue = [3 2];
btrue = [-0.01];

% Prerecorded data
load('x_bb_SVM.mat'); load('t_bb_SVM.mat'); load('j_bb_SVM.mat');
%load('x_nu_SVM.mat'); load('t_nu_SVM.mat'); load('j_nu_SVM.mat');

% Simulated data
%{
x1_0 = 0.5;
x2_0 = 0;
x0 = [x1_0,x2_0,0,0]; %last 2 stars are t and j
r = 1; tx = 0; jx = 0;

%rng(1);
%exnoise = 0.0001*rand(4,10000);
exnoise = 0.0001*[sin(1:10000); sin(5001:15000); cos(1:10000); cos(5001:15000)];

% simulation horizon
TSPAN=[0 10];
JSPAN = [0 20];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate
[t,j,x] = HyEQsolver( @f_ex1,@g_ex1,@C_ex1,@D_ex1,...
    x0,TSPAN,JSPAN,rule,options,'ode23t');
%}

% plot solution
figure(1) % position
clf
subplot(2,1,1), plotHarc(t,j,x(:,1));
grid on
ylabel('x_1 position')
subplot(2,1,2), plotHarc(t,j,x(:,2));
grid on
ylabel('x_2 velocity')

% plot phase plane
figure(2) % position
clf
plotHarcColor(x(:,1),j,x(:,2),t);
xlabel('x_1')
ylabel('x_2')
grid on

%% Example data
%s1 = [3 3 6 6; 1 -1 1 -1]'; % -1 category
%s2 = [1 0 0 -1; 0 1 -1 0]'; % +1 category
%s1 = [2 2 -2 -2; 2 -2 2 -2]';
%s2 = [1 1 1.5 1.5; 1 1.5 1 1.5]';
%[dH D] = HausdorffDist(s1,s2);

%maxx = max(max(abs(x(:,1))),max(abs(x(:,2))));
maxx = 1;
syms xc yc;
[sol_x_true, sol_y_true] = vpasolve([wtrue(1)*xc+wtrue(2)*yc+btrue==0,xc^2+yc^2==maxx^2],[xc,yc]);
sol_x_true = double(sol_x_true);
sol_y_true = double(sol_y_true);
xctrue = linspace(sol_x_true(1,end),sol_x_true(2,end),100);
yctrue = linspace(sol_y_true(1,end),sol_y_true(2,end),100);

k = 1;
for i = 1:length(j)-2
    if j(i) ~= j(i+1)
       jp(k) = i+2;
       k = k + 1;
    end
end

%m_start = 1; m_end = 500:500:4000; %m_end = 50:10:200;
%for i = m_start:length(m_end)
for i = 2:length(jp)
    x2 = x;
    if size(x2,1) > size(x2,2)
        x2 = x2';
    end
    %x = x(:,m_start:m_start+m_end(i)-1);
    x2 = x2(:,1:jp(i));

    %[w(:,i) b(i), dlim] = svm_v4(x,t,j);
    [w(:,i) b(i)] = svm_v5(x2,t,j);

    %% Hausdorff Distance
    %%{
    clear xc; clear yc;
    syms xc yc
    [sol_x, sol_y] = vpasolve([w(1,end)*xc+w(2,end)*yc+b(end)==0,xc^2+yc^2==maxx^2],[xc,yc]);
    sol_x = double(sol_x);
    sol_y = double(sol_y);
    if sign(sol_x(1)) ~= sign(sol_x_true(1))
       sol_x = -sol_x;
       sol_y = -sol_y;
    end

    xc = linspace(sol_x(1,end),sol_x(2,end),100);
    yc = linspace(sol_y(1,end),sol_y(2,end),100);
    %xd = linspace(dlim(1,1),dlim(1,2),100);
    %yd = linspace(dlim(2,1),dlim(2,2),100);

    %{
    figure(3);
    hold on
    scatter([xc(1) xc(end)],[yc(1) yc(end)],10*i);
    scatter([xd(1) xd(end)],[yd(1) yd(end)],10*i);
    %}
    dH(i) = HausdorffDist([xc;yc],[xctrue;yctrue]);
    %dH(i) = HausdorffDist([xd;yd],[xdtrue;ydtrue]);

end

%{
figure(3);
hold on
plot(xc,yc);
plot(xctrue, yctrue,'-o');
%plot(xd,yd);
%}

%{
figure(5);
Line = @(x1,x2) w(1,end)*x1+w(2,end)*x2+b(end);
x1min = min(x(1,:));x1max = max(x(1,:));
x2min = min(x(2,:));x2max = max(x(2,:));
ezplot(Line,[x1min x1max x2min x2max]);
%}

figure(6);
plot(jp(2:i),real(dH(2:end)),'b--o');
xlabel('Number of Measurements n');
ylabel('Distance');
%title('Hausdorff Distance between $$D \cap \bar{A}$$ and $$\hat{D}$$','Interpreter','latex');
title('Hausdorff Distance between $$C$$ and $$\hat{C}$$','Interpreter','latex');

