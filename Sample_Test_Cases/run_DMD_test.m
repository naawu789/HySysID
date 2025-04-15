clc; clear all; close all;

%{
A = [1 2; 3 4];
x = [7 1; 6 2];
xbar = A*x;
%[xs xv xd] = svd(x);
Abar = xbar*pinv(x)
%}

%% measurements >= dim(A)
%{
A = [1 2 3; 4 5 6];
x = [6 8 5; 8 2 6; 1 1 1]; %x = rand(3,3);
m = 3;
x2 = x(:,1:m-1);
xbar = A*x2;
[s v d] = svd(x2);
Abar = xbar*pinv(x2)
%}

%% Simple bouncing ball example
%{
A = [0 0 0; 0 -9.81 0];
x = [0 0; -8 -4; 1 1];
xbar = A*x;
Abar = xbar*pinv(x)
%}

%% Large number of states
%{
n = 10;
A = rand(n,n);
x = [rand(n-1,n+1); ones(1,n+1)];
error = zeros(1,n+1);
for m = 1:n+1
    x2 = x(:,1:m);
    xbar = A*x2;
    Abar = xbar*pinv(x2);
    error(m) = sum(sum(A-Abar));
end
plot(error); %error(end-2:end)
%}

%% Rank deficient case
%{
A = [0 0 0; 0 -9.81 0];
x = [0 0; -8 -4; 1 1];
xbar = A*x;
Abar = xbar*pinv(x)
%}

%% Thermometer
%{
x = [5 2 3; 1 0 1; 1 1 1];
AB = [-1 50 50; 0 0 0];
A = [0 4 73; 0 -1 1];
xbar = A*x;
Abar = xbar*pinv(x)
%}
