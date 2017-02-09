J = 50;  K = 50;  D = ones(J-1,K-1);
[x,y] = ndgrid(-1:2/J:1,-1:2/K:1);
T0 = exp(-30*(x.*x + y.*y));
T = diffusion(1.0,1.0,J,K,D,D,D,D,T0,0.05);
