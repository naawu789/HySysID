clc; clear all; close all;

%Numerical example
%{
AB = [-1 2 0; 3 4 0];
HJ = [0.5 -0.3 0; -0.4 0.5 0];
load('x_nu_DMD.mat'); load('t_nu_DMD.mat'); load('j_nu_DMD.mat');
%}

%Bouncing ball
%%{
AB = [0 1 0; 0 0 -9.81];
HJ = [0 0 0; 0 -0.8 0];
load('x_bb_08_DMD.mat'); load('t_bb_08_DMD.mat'); load('j_bb_08_DMD.mat');
%}

%Thermostat
%{
AB = [-1 50 50; 0 0 0]; HJ = [0 4 73; 0 -1 1];
load('x_th_DMD.mat'); load('t_th_DMD.mat'); load('j_th_DMD.mat');
%}

if size(x,1) > size(x,2)
    x = x';
end
%x = x(:,1:200);

%% Split continuous and discrete data
x_cont = []; t_cont = []; xc = 0;
x_disc = []; t_disc = []; xd = 0;

%Algorithm 1 (j counter)
%%{
for i = 1:size(x,2)-1
    if j(i+1) == j(i)
        xc = xc + 1;
        x_cont(:,xc) = x(:,i);
        t_cont(xc) = t(i);
        j_cont(xc) = i;
    elseif j(i+1) ~= j(i)
        xd = xd + 1;
        x_disc(:,xd) = x(:,i);
        x2_disc(:,xd) = x(:,i+1); %discrete X'
        t_disc(xd) = t(i);
        j_disc(xd) = i;
    end
end
%}

%Algorithm 2 (check for infinite derivative)
%{
for i = 1:size(x,2)-1
    xgrad(:,i) = (x(:,i+1) - x(:,i))/(t(i+1)-t(i));
    %if ~isinf(norm(xgrad(:,i)))
    if j(i+1) == j(i)
        xc = xc + 1;
        x_cont(:,xc) = x(:,i);
        t_cont(xc) = t(i);
        j_cont(xc) = i;
    %elseif isinf(norm(xgrad(:,i)))
    elseif j(i+1) ~= j(i)
        xd = xd + 1;
        x_disc(:,xd) = x(:,i);
        x2_disc(:,xd) = x(:,i+1); %discrete X'
        t_disc(xd) = t(i);
        j_disc(xd) = i;
    end
end
%{
figure;
subplot(212);
plot(t(1:end-1),xgrad);
hold on;
subplot(211);
plot(t(1:end-1),x(:,1:end-1));
%}
%}

%% Apply bias trick to X
x_cont_bias = [x_cont(:,1:end-1); ones(1,xc-1)]; %v1
%x_cont_bias = [x_cont(:,1:end); ones(1,xc)]; %v2
x_disc_bias = [x_disc(:,1:end-1); ones(1,xd-1)];

%% Compute or generate X'
x2_cont = [AB; 0 0 1]*[x_cont_bias]; %v1
%{
x2_cont = []; %v2
for i = 2:length(x_cont_bias)-1
   x2_cont(:,i) = [(x_cont_bias(1:end-1,i+1) - x_cont_bias(1:end-1,i))/(t_cont(i+1)-t_cont(i)); 1];
end
figure;
subplot(211); plot(t_cont(1:end-1), x2_cont(1,:));
subplot(212); plot(t_cont(1:end-1), x2_cont(2,:));
%}
x2_disc = [x2_disc; ones(1,xd)];
%{
%x2_disc computed earlier above
%x2_cont2 = gradient(x_cont);
%}

%% CONT Optionally, truncate data
for r1_cont = 1:xc-1
    x2_cont_trun = x2_cont(:,1:r1_cont);
    x_cont_trun = x_cont_bias(:,1:r1_cont);

%% CONT Compute SVD and pick reduced order r
    [s_cont, v_cont, d_cont] = svd(x_cont_trun);

%% CONT DMD to compute [A B]
    ABbar = x2_cont_trun*d_cont*pinv(v_cont)*s_cont';
    error_cont(r1_cont) = norm([AB; 0 0 1]-ABbar);
end

%% DISC Optionally, truncate data
for r1_disc = 1:xd-1
    x2_disc_trun = x2_disc(:,1:r1_disc);
    x_disc_trun = x_disc_bias(:,1:r1_disc);

%% DISC Compute SVD and pick reduced order r
    [s_disc, v_disc, d_disc] = svd(x_disc_trun);

%% DISC DMD to compute [H J]
    HJbar = x2_disc_trun*d_disc*pinv(v_disc)*s_disc';
    error_disc(r1_disc) = norm([HJ; 0 0 1]-HJbar);
end

%% Plotting
ABbar
HJbar
eig(1/r1_cont*x_cont_trun*x_cont_trun')
eig(1/r1_disc*x_disc_trun*x_disc_trun')
figure;
subplot(2,1,1); plot(error_cont);
title('e_{c} vs n');
ylabel('Magnitude');
subplot(2,1,2); plot(error_disc);
title('e_{d} vs n');
xlabel('Number of measurements');
ylabel('Magnitude');
