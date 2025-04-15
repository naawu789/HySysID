function xplus = g_ex1_2(x)

% state
x1 = x
x2 = x(2);
lambda = 0.8;
xplus = [0 ; -lambda*x2];
end
