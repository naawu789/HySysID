function output = noise(x,y)
global exnoise r
output = exnoise(x:x+1,r);
if r >= 9000
    r = 1;
end
%r = r + 1;
end