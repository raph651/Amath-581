%% HW 4 Raphael Liu
% computing Vorticity Streamfunction Equations
% using the method of i) A\b  ii) LU decomposition iii) FFT
% a)
L = 5;
n = 32;

xspan = linspace(-L, L, n+1);
yspan = linspace(-L, L, n+1);

xx = xspan(1:n);
yy = yspan(1:n);
[X, Y] = meshgrid(xx, yy);
w0 = exp(-X.^2-Y.^2/20);
w_initial = reshape(w0, [n^2, 1]);

tspan = 0:0.5:4;
h = xx(2) - xx(1);

e1 = ones(n^2, 1);
B = [-4 * e1, e1, e1, e1, e1, e1, e1];
d = [0, 1, -1, n, -n, (n - 1) * n, -((n - 1) * n)];
A = spdiags(B, d, n^2, n^2);
for i = 1:n
    A((i - 1)*n+1, i*n) = 1;
    A(i*n, (i - 1)*n+1) = 1;
    if i <= n - 1
        A(i*n, i*n+1) = 0;
        A(i*n+1, i*n) = 0;
    end
end
A(1, 1) = 2;


B_dx = [e1, -e1, -e1, e1];
d_dx = [n, -n, (n - 1) * n, -((n - 1) * n)];
B = spdiags(B_dx, d_dx, n^2, n^2);

B_dy = [e1, -e1];
d_dy = [1, -1];
C = spdiags(B_dy, d_dy, n^2, n^2);
for i = 1:n
    C((i - 1)*n+1, i*n) = -1;
    C(i*n, (i - 1)*n+1) = 1;
    if i <= n - 1
        C(i*n, i*n+1) = 0;
        C(i*n+1, i*n) = 0;
    end
end
v = 0.001;
A = A ./ h^2;
B = B ./ (2 * h);
C = C ./ (2 * h);

% i) using A\b
%tic
[T, W1] = ode45(@(t, w) rhspde1(t, w, v, A, B, C), tspan, w_initial);
%toc
A1 = W1;

% temp1=W1(end,:);
% temp1=reshape(temp1,n,n);
% pcolor(temp1);
%%
% ii) using LU
% tic
[T, W2] = ode45(@(t, w) rhspde2(t, w, v, A, B, C), tspan, w_initial);
% toc
A2 = W2;

%%
%b)
[T, W3] = ode45(@(t, w) rhspde3(t, w, v, A, B, C), tspan, w_initial);
A3 = W3;

% temp3=W3(end,:);
% temp3=reshape(temp3,n,n);
% pcolor(temp3);

[T, W4] = ode45(@(t, w) rhspde4(t, w, v, A, B, C), tspan, w_initial);
A4 = W4;

%c) 
%%
[T, W5] = ode45(@(t, w) rhspde5(t, w, v, A, B, C), tspan, w0);
A5 = W5;
% temp5=W5(end,:);
% temp5=reshape(temp5,n,n);
% pcolor(temp5);