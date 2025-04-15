function inside = D_est(x)

global s gamma dimpsi u D
% state
psi = x(1:dimpsi);
thetarows = x(dimpsi+1:end-3);
t = x(end-2);
k = x(end-1);
j = x(end);

if D(psi) || t >= k+s
    inside = 1;
else
    inside = 0;
end
end
