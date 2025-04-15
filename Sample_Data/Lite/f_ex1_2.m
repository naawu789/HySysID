function xdot = f_ex1_2(x)

% state
x2 = x(2);

% differential equations
gamma = -9.81;
xdot = [x2 ; gamma];
end
