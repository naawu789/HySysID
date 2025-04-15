clc; clear all; close all;

global s gamma_c gamma_d dimpsi dimt A B H J u C D

%pick example
[s, gamma_c, gamma_d, psi0, A, B, H, J, u, C, D, TSPAN, JSPAN, vmin, vmax] = pickexample(1);

plotbounds = 0;

dimpsi = size(A,1); dimu = size(u,1); dimt = dimpsi + dimu;

theta_c = [A B];                      theta_d = [H J];                     % parameter
hat_theta_c = zeros(dimpsi,dimt);     hat_theta_d = zeros(dimpsi,dimt);    % estimate

for i = 1:dimpsi
    hat_theta_c_rows(1+(i-1)*dimt:dimt+(i-1)*dimt) = hat_theta_c(i,:);
    hat_theta_d_rows(1+(i-1)*dimt:dimt+(i-1)*dimt) = hat_theta_d(i,:);
end

% initial conditions
t0 = 0; k0 = 0; j0 = 0;
x0 = [psi0, hat_theta_c_rows, hat_theta_d_rows, t0, k0, j0];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;
options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate
[t,j,x] = HyEQsolver( @f_est,@g_est,@C_est,@D_est,...
    x0,TSPAN,JSPAN,rule,options,'ode23t');

for i = 1:dimpsi
    hat_theta_c(i,:) = x(end,dimpsi+1+(i-1)*dimt               : dimpsi+dimt+(i-1)*dimt);
    hat_theta_d(i,:) = x(end,dimpsi+1+(i-1)*dimt + dimpsi*dimt : dimpsi+dimt+(i-1)*dimt + dimpsi*dimt);
end

hat_theta_c
hat_theta_d
for i = 1:size(t,1)
    error_norm_c(i) = norm(x(i,1+dimpsi:dimpsi+dimpsi*dimt) - reshape(theta_c',[1,dimpsi*dimt]));
    error_norm_d(i) = norm(x(i,1+dimpsi+dimpsi*dimt:dimpsi+2*dimpsi*dimt) - reshape(theta_d',[1,dimpsi*dimt]));
end
error_norm_c(end)
error_norm_d(end)

%% PE, learning constants test
%psi_eig = eig(x(1:end-2,1:dimpsi)'*x(1:end-2,1:dimpsi));
%psi_m2 = max(x(1:end-2,1:dimpsi).*x(1:end-2,1:dimpsi));

%% Plots
% plot solution
%%{
figure(1) % position
clf
varm = {'x_1 position', 'x_2 velocity', 'x_3 state'};
for i = 1:dimpsi
    subplot(dimpsi,1,i), plotHarc(t,j,x(:,i));
    grid on
    ylabel(varm(i));
end

figure(2)
subplot(2,1,1);
plotHarc(t,j,error_norm_c);
%title('Norm of error $\bar\theta_c$ vs Time','interpreter','latex');
xlabel('Time'); ylabel('$\hat\theta_c$','interpreter','latex');
if plotbounds == 1
    hold on; plotHarc(t,j,boundc,[j(1),j(end)],[],modificatorJ);
end
subplot(2,1,2);
plotHarc(t,j,error_norm_d);
%title('Norm of error $\bar\theta_d$ vs Time','interpreter','latex');
xlabel('Time'); ylabel('$\hat\theta_d$','interpreter','latex');
if plotbounds == 1
    hold on; plotHarc(t,j,boundd,[j(1),j(end)],[],modificatorJ);
end
%}
