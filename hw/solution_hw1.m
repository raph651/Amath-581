
%% Homework 1 Raphael Liu
% Initial value problems for ODEs : Forward-Euler, Heun.
% van der Pol oscillator in ode23 ode45 and ode113
% solve Fitzhugh neuron model 
%% Problem 1

dt = [2^-2, 2^-3, 2^-4, 2^-5, 2^-6, 2^-7, 2^-8];
y0 = pi / sqrt(2);
Error_FE = zeros(1, length(dt));
Error_HE = zeros(1, length(dt));
f = @(a, b) -3 * a * sin(b);

for i = 1:7

    t = [0:dt(i):5];
    y_true{i} = pi .* exp(3*(cos(t) - 1)) / sqrt(2);
    y_FE{i} = zeros(1, length(y_true{1, i}));
    y_HE{i} = zeros(1, length(y_true{1, i}));
    y_FE{1, i}(1) = y0;
    y_HE{1, i}(1) = y0;

    for k = 2:length(t)
        y_FE{1, i}(k) = y_FE{1, i}(k - 1) + dt(i) * f(y_FE{1, i}(k - 1), t(k-1));
        y_HE{1, i}(k) = y_HE{1, i}(k - 1) + 0.5 * dt(i) * (f(y_HE{1, i}(k - 1), t(k-1)) + f(y_HE{1, i}(k - 1)+dt(i)*f(y_HE{1, i}(k - 1), t(k-1)), t(k)));
    end

    Error_FE(i) = mean(abs(y_FE{1, i}-y_true{1, i}));
    Error_HE(i) = mean(abs(y_HE{1, i}-y_true{1, i}));

end
figure(1)
plot(log(dt), log(Error_FE));
p_FE = polyfit(log(dt), log(Error_FE), 1);
hold on;
plot(log(dt), log(Error_HE));
p_HE = polyfit(log(dt), log(Error_HE), 1);
title("log(dt) vs log(E)");
xlabel("log(dt)");
ylabel("log(E)");
legend("FE-forwardEuler", "HE-Heun");
hold off

A1 = y_FE{1, 7}';
A2 = Error_FE;
A3 = p_FE(1)
A4 = y_HE{1, 7}';
A5 = Error_HE;
A6 = p_HE(1)

%% Problem 2
% a)

epsilon = 0.1;
vdp = @(t, y) [y(2); -epsilon * (y(1)^2 - 1) * y(2) - y(1)];
[t1, y1] = ode45(vdp, [0:0.5:32], [sqrt(3); 1]);

epsilon = 1;
vdp = @(t, y) [y(2); -epsilon * (y(1)^2 - 1) * y(2) - y(1)];
[t2, y2] = ode45(vdp, [0:0.5:32], [sqrt(3); 1]);

epsilon = 20;
vdp = @(t, y) [y(2); -epsilon * (y(1)^2 - 1) * y(2) - y(1)];

[t3, y3] = ode45(vdp, [0:0.5:32], [sqrt(3); 1]);


% b)
epsilon = 1;
vdp = @(t, y) [y(2); -epsilon * (y(1)^2 - 1) * y(2) - y(1)];
TOL_start = 1e-4;
tspan = [0 32];
TOL = zeros(1, 7);
stepT = zeros(3, 7);

for i = 1:7

    TOL(i) = TOL_start * 10^(1-i);
    options = odeset('AbsTol', TOL(i), 'RelTol', TOL(i));
    [T, Y] = ode45(vdp, tspan, [2, pi^2], options);
    stepT(1, i) = mean(diff(T));
    [T, Y] = ode23(vdp, tspan, [2, pi^2], options);
    stepT(2, i) = mean(diff(T));
    [T, Y] = ode113(vdp, tspan, [2, pi^2], options);
    stepT(3, i) = mean(diff(T));

end
figure(2)
plot(log(stepT), log(TOL));
title('log(dt) vs log(TOL)');
xlabel('log(dt)');
ylabel('log(TOL)');
legend('ode45', 'ode23', 'ode113');

p_ode45 = polyfit(log(stepT(1, :)), log(TOL), 1);
p_ode23 = polyfit(log(stepT(2, :)), log(TOL), 1);
p_ode113 = polyfit(log(stepT(3, :)), log(TOL), 1);


A7 = [y1(:, 1), y2(:, 1), y3(:, 1)];
A8 = p_ode45(1);
A9 = p_ode23(1);
A10 = p_ode113(1);

%% Problem 3

% --> first part <--
a1 = 0.05;
a2 = 0.25;
b = 0.01;
c = 0.01;
I = 0.1;
d12 = -0.1;
d21 = 0.1;


fitzhugh = @(t, y) [-y(1)^3 + (1 + a1) * y(1)^2 - a1 * y(1) - y(2) + I + d12 * y(3); b * y(1) - c * y(2); -y(3)^3 + (1 + a2) * y(3)^2 - a2 * y(3) - y(4) + I + d21 * y(1); b * y(3) - c * y(4)];

[t, y] = ode15s(fitzhugh, [0: 0.5: 100], [0.1; 0; 0.1; 0]);

%hold on
%plot(t, y);
%plot(t,exp(y))
%plot(t,sin(y))
%legend("v1", "w1", "v2", "w2");

% from the above different graphical representations, we can see that v1 and v2
% are coupled, varying in the same manner, and that w1 and w2 are coupled.

% --> second part <--

d12 = [0, 0, -0.1, -0.3, -0.5];
d21 = [0, 0.2, 0.2, 0.2, 0.2];
y_solution = cell(1, 5);

for i = 1:5
    fitzhugh = @(t, y) [-y(1)^3 + (1 + a1) * y(1)^2 - a1 * y(1) - y(2) + I + d12(i) * y(3); b * y(1) - c * y(2); -y(3)^3 + (1 + a2) * y(3)^2 - a2 * y(3) - y(4) + I + d21(i) * y(1); b * y(3) - c * y(4)];
    [t, y] = ode15s(fitzhugh, [0: 0.5: 100], [0.1; 0; 0.1; 0]);

    y_solution{1, i} = y;
    y_solution{1, i}(:, 2) = y(:, 3);
    y_solution{1, i}(:, 3) = y(:, 2);
end

A11 = y_solution{1, 1};
A12 = y_solution{1, 2};
A13 = y_solution{1, 3};
A14 = y_solution{1, 4};
A15 = y_solution{1, 5}

%%