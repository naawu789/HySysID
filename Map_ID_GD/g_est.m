function xplus = g_est(x)

global s gamma_c gamma_d dimpsi dimt u A B H J D

% state
hat_theta_c_rows = x(dimpsi+1             : dimpsi+    dimpsi*dimt);
hat_theta_d_rows = x(dimpsi+1+dimpsi*dimt : dimpsi+ 2*(dimpsi*dimt));

t = x(end-2); k = x(end-1); j = x(end);
psi = x(1:dimpsi);

for i = 1:dimpsi
    hat_theta_c(i,:) = hat_theta_c_rows(1+(i-1)*dimt:dimt+(i-1)*dimt);
    hat_theta_d(i,:) = hat_theta_d_rows(1+(i-1)*dimt:dimt+(i-1)*dimt);
end

if t >= k+s
    bar_theta_c = hat_theta_c - [A B];
    hat_theta_c = hat_theta_c - s*gamma_c*psi*(bar_theta_c'*psi)';
end
if D(psi)
    bar_theta_c = hat_theta_c - [A B];
    bar_theta_d = hat_theta_d - [H J];
    hat_theta_c = hat_theta_c - gamma_d*psi*(bar_theta_c'*psi)'/(1+gamma_d*psi'*psi);
    hat_theta_d = hat_theta_d - gamma_d*psi*(bar_theta_d'*psi)'/(1+gamma_d*psi'*psi);
end

for i = 1:dimpsi
    hat_theta_c_rows(1+(i-1)*dimt:dimt+(i-1)*dimt) = hat_theta_c(i,:);
    hat_theta_d_rows(1+(i-1)*dimt:dimt+(i-1)*dimt) = hat_theta_d(i,:);
end

if t >= k+s
    xplus = [psi; hat_theta_c_rows; hat_theta_d_rows; t; k+s; j];
end
if D(psi)
    xplus = [H*psi+J*u(t); hat_theta_c_rows; hat_theta_d_rows; t; k; j+1];
end
end
