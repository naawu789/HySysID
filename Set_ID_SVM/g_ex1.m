function xplus = g_ex1(x)

% state
x1 = x(1);
x2 = x(2);
x3 = x(3); %t
x4 = x(4); %j
global HJ r jx
%xplus = HJ*[x1;x2;1]; %+ 0.01*rand(size(HJ,1),1);
xplus = [HJ*[x1;x2;1] + noise(3,r);x3;x4+1];
if x4 ~= jx
    r = r + 1;
end
%xplus = [0.5*x1-0.3*x2; -0.4*x1+0.5*x2];
%xplus = [-x1+2*x2; 3*x1+4*x2]
%xplus = [-1 0 0; 0 -1 0]*[x1;x2;1];
end
