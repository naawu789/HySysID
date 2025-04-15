function xdot = f_est(x)

global s dimpsi dimt u A B
% state
psi = x(1:dimpsi);
thetarows = x(dimpsi+1:end-3);
t = x(end-2);
k = x(end-1);
j = x(end);

for i = 1:dimpsi
    theta(i,:) = thetarows(1+(i-1)*dimt:dimt+(i-1)*dimt);
end
for i = 1:dimpsi
    theta2rows(1+(i-1)*dimt:dimt+(i-1)*dimt) = theta(i,:);
end

% differential equations
xdot = [A*psi + B*u(t); zeros(size(theta2rows))'; zeros(size(theta2rows))'; 1; 0; 0];
end
