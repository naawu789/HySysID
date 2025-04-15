function [value] = C_est(x)

global s gamma dimpsi u C
% state
psi = x(1:dimpsi);
thetarows = x(dimpsi+1:end-3);
t = x(end-2);
k = x(end-1);
j = x(end);

if C(psi)
    value = 1;
else
    value = 0;
end
end
