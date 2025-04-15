function [s, gamma_c, gamma_d, psi0, A, B, H, J, u, C, D, TSPAN, JSPAN, vmin, vmax] = pickexample(x)

switch x
    case 1 %bouncing ball
        s = 0.005; psi0 = [1 0];
        gamma_c = 1/(s*100); gamma_d = 1/(5000*s);
        TSPAN=[0 10]; JSPAN = [0 2000];
        %TSPAN = [0 30]; JSPAN = [0 5000]; s = 0.05;
        A = [0 1; 0 0]; B = [0; -9.8];
        H = [0 0; 0 -1]; J = [0; 0];
        u = @(t) 1;
        C = @(psi) psi(1) >= 0;
        D = @(psi) (psi(1) <= 0 && psi(2) <= 0);
        vmin = 546; vmax = 1633;
    case 2 %temperature control TBD
        s = 0.05; psi0 = [1 0];
        gamma_c = 1/(s*1); gamma_d = 1/(1*s);
        TSPAN=[0 10]; JSPAN = [0 2000];
        a = 1; Td = 90; Tr = 0;
        A = [-a -a*Tr+a*Td; 0 0]; B = [a*Tr; 0];
        H = [1 0; 0 -1]; J = [0; 1];
        u = @(t) 1;
        C = @(psi) 1;
        D = @(psi) (psi(1) <= 70 && psi(2) == 0) || (psi(1) >= 80 && psi(2) == 1);
        vmin = 1; vmax = 1;
    case 3 %academic example 1
        s = 0.005; psi0 = [0 1];
        gamma_c = 1/(50*s); gamma_d = 1/(0.05*s);
        TSPAN = [0 10]; JSPAN = [0 2000];
        A = [0 1; 0 0]; B = [1; 0];
        H = [0 0; 0 -1]; J = [0; 1];
        u = @(t) 1;
        C = @(psi) (1);
        D = @(psi) (psi(1) >= 1);
        vmin = 1; vmax = 603;
    case 4 %clock skew TBD
        s = 0.005; psi0 = [0 0];
        gamma_c = 1/(1*s); gamma_d = 1/(1*s);
        TSPAN = [0 10]; JSPAN = [0 2000];
        epsilon = 0.05;
        A = [0 1+epsilon; 0 0]; B = [0 0; 0 0];
        H = [0 0; 0 -1]; J = [0; 1];
        u = @(t) [2*floor(t)-floor(2*t)+1; 1];
        C = @(psi) 1;
        D = @(psi) psi(2) ~= u;
        vmin = 1; vmax = 1;
    case 5 %juggilng bouncing ball TBD
        s = 0.005;
        gamma_c = 1/(s*1); gamma_d = 1/(1*s);
        TSPAN = [0 10]; JSPAN = [0 2000];
        psi0 = [3 0 0];
        A = [0 1 0; 0 0 0; 0 0 0]; B = [0; -9.8; 0];
        H = [0 0 1; 0 -1 0; 0 0 -1]; J = [0; 0; 1];
        u = @(t) 1;
        C = @(psi) (1);
        D = @(psi) (psi(1) <= psi(3) && psi(2) <= 0);
        vmin = 1; vmax = 1;
    case 6 %academic system v2 TBD
        s = 0.01;
        gamma_c = 1/(1*s); gamma_d = 1/(1*s);
        TSPAN = [0 20]; JSPAN = [0 2000];
        psi0 = [0 0.0025];
        A = [-0.7 3; -2 -0.5]; B = [1; 0];
        H = [-1 0; 0 -1]; J = [0; 1];
        u = @(t) 1;
        C = @(psi) 1;
        D = @(psi) 3*psi(1) + 6*psi(2) < 0;
        vmin = 1; vmax = 1;
    case 7 %pressure mounter TBD
        k = 25; m = 0.5; a = 0.5; c = 1;
        A = [0 1; -k/m -c/m]; B = [0; 1/m];
        H = [1 0; 0 -a/c]; J = [0; 0];
        u = @(t) 1;
end
