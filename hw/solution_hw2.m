
%% Homework 2 Raphael Liu
% Quantum Harmonic Oscillator/Boundary Value Problems using shooting scheme
% using Direct Method, solve eigen equations. Bootstrap approach
% double shooting/continuous shooting for different parameters
%% Problem 1

L = 4;
xp = -L:0.1:L;
A = 1;
Tol = 10^-4;

eigen_functions = zeros(5, length(xp));
eigen_values = zeros(1, 5);

e0_start = 0;
de0_start = 1;

for modes = 1:5

    e0 = e0_start;
    de0 = de0_start;

    for i = 1:1000
        y0 = [A; sqrt(L^2-e0) * A];
        [x, y] = ode45(@(x, y) func1_hw2(x, y, e0), xp, y0);
        if (abs(y(end, 2)+sqrt(L^2-e0)*y(end, 1)) < Tol)
            e0;
            break;
        end
        if (-1)^(modes + 1) * (y(end, 2) + sqrt(L^2-e0) * y(end, 1)) > 0
            e0 = e0 + de0;
        else
            e0 = e0 - de0 / 2;
            de0 = de0 / 2;
        end
    end

    e0_start = e0 + 0.1;
    eigen_functions(modes, :) = y(:, 1) ./ sqrt(trapz(x, y(:, 1).^2));
    eigen_values(modes) = e0;
    hold on
    plot(x, abs(eigen_functions(modes, :)));
end

A1 = abs(eigen_functions(1, :))';
A2 = abs(eigen_functions(2, :))';
A3 = abs(eigen_functions(3, :))';
A4 = abs(eigen_functions(4, :))';
A5 = abs(eigen_functions(5, :))';
A6 = eigen_values';

%% Problem 2
L = 4;
xp = -L:0.1:L;
A = zeros(length(xp)-2);


for i = 2:length(A) - 1
    A(i, i) = 2 + xp(i+1)^2 * 0.1^2;
    A(i, i-1) = -1;
    A(i, i+1) = -1;
end


A(1, 1) = 2 / 3 + (-L + 0.1)^2 * 0.1^2;
A(1, 2) = -2 / 3;
A(end, end) = 2 / 3 + (L - 0.1)^2 * 0.1^2;
A(end, end-1) = -2 / 3;


[V, D] = eig(A);
ev = diag(D) ./ 0.1^2;
[ev, idx] = sort(ev, "ascend");
ev = ev(1:5);
V = V(:, idx(1:5));
VW = zeros(81, 5);

for i = 1:5
    VW(1, i) = (4 * V(1, i) - V(2, i)) / (2 * 0.1 * sqrt(L^2-ev(i)) + 3);
    VW(end, i) = (4 * V(end, i) - V(end-1, i)) / (2 * 0.1 * sqrt(L^2-ev(i)) + 3);
    VW(2:end-1, i) = V(:, i);
end


for i = 1:5
    VW(:, i) = abs(VW(:, i)./sqrt(sum(VW(:, i).^2*0.1)));
end

figure();
hold on
plot(xp, VW);

A7 = VW(:, 1);
A8 = VW(:, 2);
A9 = VW(:, 3);
A10 = VW(:, 4);
A11 = VW(:, 5);
A12 = ev;

%% Problem 3


A_start = 1;
dA_start = 1;
e0_start = 0.2;
%e0_start = 0;
de0_start = 1;

L = 2;
xp = [-L :0.1 :L];
Tol = 10^-4;
eigen_functions = zeros(4, length(xp));
eigen_values = zeros(1, 4);


gamma = [0.05, 0.05, -0.05, -0.05];

%%
figure();

for modes = 1:4

    if (modes == 3)
        e0_start = 0.2;
        de0_start = 1;
    end

    e0 = e0_start;
    de0 = de0_start;

    for ite = 1:15
        % if (ite==1 && modes==1)
        if (ite == 1 && (modes == 1 || modes == 3))
            A = A_start;
            dA = dA_start;
        else
            dA = dA * 2^(round(i/10));
            de0 = de0 * 2^(round(j/3));
        end

        for i = 1:1001

            y0 = [A; sqrt(L^2-e0) * A];
            [x, y] = ode45(@(x, y) func3(x, y, e0, gamma(modes)), xp, y0);
            norm_y = sqrt(trapz(x, y(:, 1).^2));

            if (abs(norm_y-1) < Tol)
                A;
                break;

            else if (norm_y < 1)
                    A = A + dA;
            else
                A = abs(A-dA/1.2);
                dA = dA / 1.2;
                %A=A-dA/2;
                %dA=dA/2;
            end
            end
        end

        for j = 1:1001
            y0 = [A; sqrt(L^2-e0) * A];

            [x, y] = ode45(@(x, y) func3(x, y, e0, gamma(modes)), xp, y0);

            if (abs(y(end, 2)+sqrt(L^2-e0)*y(end, 1)) < Tol)
                e0;
                break;
            end
            if (-1)^(modes + 1) * (y(end, 2) + sqrt(L^2-e0) * y(end, 1)) > 0
                e0 = e0 + de0;
            else
                e0 = e0 - de0 / 2;
                de0 = de0 / 2;

            end
        end
        ite = ite + 1;
    end


    e0_start = e0 + 0.2;

    hold on
    plot(x, y(:, 1))


    eigen_functions(modes, :) = y(:, 1);
    eigen_values(modes) = e0;


end

A13 = abs(eigen_functions(1, :))';
A14 = abs(eigen_functions(2, :))';
A15 = eigen_values(1:2)';
A16 = abs(eigen_functions(3, :))';
A17 = abs(eigen_functions(4, :))';
A18 = eigen_values(3:4)';