%% HW 4 As cool as it gets Raphael Liu
% a)
L = 10;
n = 128;

xspan = linspace(-L, L, n+1);
yspan = linspace(-L, L, n+1);

xx = xspan(1:n);
yy = yspan(1:n);
[X, Y] = meshgrid(xx, yy);


tspan = 0:0.5:25;
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

% 1) initial condition for two oppositely charged Gaussian vortices next to
% each other (space= 1)
w1_a = -3 * exp(-(X - 0.75).^2/6-Y.^2);
w1_b = 3 * exp(-(X + 0.75).^2/6-Y.^2);
% 2) initial condition for two same charged Gaussian vortices next to each
% other (space =2)

w2_a = exp(-(X - 1).^2-Y.^2/20);
w2_b = exp(-(X + 1).^2-Y.^2/20);

% 3) initial condition for two pairs of oppositely charged Gaussian
% vortices, with different amplitudes. ( the center of the four vortices is origin, x space =6; y space = 4; )

w3_a = -2*exp(-(X - 1).^2-(Y+1).^2/3);
w3_b = 1.5*exp(-(X -1).^2-(Y-1).^2/3);
w3_c = 2*exp(-(X+1).^2-(Y + 1).^2/3);
w3_d = -1.5*exp(-(X+1).^2-(Y-1).^2/3);

% 4) initial condition for 15 alternatingly charged Gaussian
% vortices, with different strength and charge
w4_a1 = -2 * exp(-(X - 4).^2-(Y - 2).^2);
w4_b1 = 1 * exp(-(X - 4).^2-(Y - 1).^2);
w4_c1 = -2 * exp(-(X - 4).^2-(Y).^2);
w4_d1 = 1 * exp(-(X - 4).^2-(Y + 1).^2);
w4_e1 = -2 * exp(-(X - 4).^2-(Y + 2).^2);
w4_a2 = 3 * exp(-X.^2-(Y - 2).^2);
w4_b2 = -2 * exp(-X.^2-(Y - 1).^2);
w4_c2 = 5 * exp(-X.^2-(Y).^2);
w4_d2 = -2 * exp(-X.^2-(Y + 1).^2);
w4_e2 = 3 * exp(-X.^2-(Y + 2).^2);
w4_a3 = -2 * exp(-(X + 4).^2-(Y - 2).^2);
w4_b3 = 1 * exp(-(X + 4).^2-(Y - 1).^2);
w4_c3 = -2 * exp(-(X + 4).^2-(Y).^2);
w4_d3 = 1 * exp(-(X + 4).^2-(Y + 1).^2);
w4_e3 = -2 * exp(-(X + 4).^2-(Y + 2).^2);

w1 = w1_a + w1_b;
w2 = w2_a + w2_b;
w3 = w3_a + w3_b + w3_c + w3_d;
w4 = w4_a1 + w4_b1 + w4_c1 + w4_d1 + w4_e1 + w4_a2 + w4_b2 + w4_c2 + w4_d2 + w4_e2 + w4_a3 + w4_b3 + w4_c3 + w4_d3 + w4_e3;

% 5)



w5=w1-w2+2*w3;




% using fft (the fastest)
%%
% 1)
[T, W1] = ode45(@(t, w) rhspde5(t, w, v, A, B, C), tspan, w1);
v1=VideoWriter('as cool as it gets 1.avi');
v1.FrameRate=10;
open(v1);
for i = 1:length(tspan)
    temp = W1(i, :);
    temp = reshape(temp, n, n);

    pcolor( temp);
    title('cool 1')
    shading interp
    frame =getframe(gcf);
    writeVideo(v1,frame);
end
close(v1);
%%
% 2)
[T, W2] = ode45(@(t, w) rhspde5(t, w, v, A, B, C), tspan, w2);
v2=VideoWriter('as cool as it gets 2.avi');
v2.FrameRate=10;
open(v2);
for i = 1:length(tspan)
    temp = W2(i, :);
    temp = reshape(temp, n, n);

    pcolor( temp);
    title('cool 2')
    shading interp
    frame =getframe(gcf);
    writeVideo(v2,frame);
end
close(v2);
%%
% 3)
tspan = 0:0.5:75;

[T, W3] = ode45(@(t, w) rhspde5(t, w, v, A, B, C), tspan, w3);
v3=VideoWriter('as cool as it gets 3.avi');
v3.FrameRate=10;
open(v3);
for i = 1:length(tspan)
    temp = W3(i, :);
    temp = reshape(temp, n, n);

    pcolor( temp);
    title('cool 3')
    shading interp
    frame =getframe(gcf);
    writeVideo(v3,frame);
end
close(v3);
%%
% 4)
[T, W4] = ode45(@(t, w) rhspde5(t, w, v, A, B, C), tspan, w4);
v4=VideoWriter('as cool as it gets 4.avi');
v4.FrameRate=10;
open(v4);
for i = 1:length(tspan)
    temp = W4(i, :);
    temp = reshape(temp, n, n);

    pcolor( temp);
    title('cool 4')
    shading interp
    frame =getframe(gcf);
    writeVideo(v4,frame);
end
close(v4);

%%
% 5)
[T, W5] = ode45(@(t, w) rhspde5(t, w, v, A, B, C), tspan, w5);
v5=VideoWriter('as cool as it gets 5.avi');
v5.FrameRate=10;
open(v5);
for i = 1:length(tspan)
    temp = W5(i, :);
    temp = reshape(temp, n, n);

    pcolor( temp);
    title('cool 5')
    shading interp
    frame =getframe(gcf);
    writeVideo(v5,frame);
end
close(v5);
