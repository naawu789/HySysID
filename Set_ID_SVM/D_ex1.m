function inside = D_ex1(x)

x1 = x(1);
x2 = x(2);

global wtrue btrue

if wtrue(1)*x1 + wtrue(2)*x2 < btrue %(x1 > 0.25 || x2 < -0.25)
    inside = 1;
else
    inside = 0;
end
end
