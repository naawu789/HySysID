function [value] = C_ex1(x)

x1 = x(1);
x2 = x(2);

global wtrue btrue

if wtrue(1)*x1 + wtrue(2)*x2 >= btrue %x1 <= 0.25 && x2 >= -0.25
    value = 1;
else
    value = 0;
end
end
