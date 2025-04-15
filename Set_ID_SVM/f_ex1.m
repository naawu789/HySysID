function xdot = f_ex1(x)

% state
x1 = x(1);
x2 = x(2);
x3 = x(3); % t

global AB r tx
%xdot = AB*[x1;x2;1]; % + 0.01*rand(size(AB,1),1);
xdot = [AB*[x1;x2;1]+ noise(1,r);1;0];% ;
if x3 ~= tx
    r = r + 1;
end
tx = x3;
% differential equations
%xdot = [x2 ; gamma];
%xdot = [-3*x1+2*x2; 7*x1-5*x2];
%xdot = [-x1+2*x2; 3*x1+4*x2];
%xdot = [-0.7 3 1; -3 -0.7 0]*[x1; x2; 1];
%xdot = [-0.2 3 0; -3 -0.2 0]*[x1; x2; 1];
end
