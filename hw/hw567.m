%% Problem 7 Raphael Liu
% a)
n = 0:200; % pre-define 200 terms, 200 should be enough
z = 0.9;
f = z.^(n);
F = zeros(1, 200);
TOL = 10^-6;
for i = 1:200
    F(i) = sum(f(1:i)); % find the sum of the sequence
    if abs(F(i)-10) < TOL % if the sum is within tolerance, return the number of terms
        i
        break
    end
end

% 153 terms are needed to reach the tolerance, N = 152

%b)
%to calculate G(0.9)
N = i - 1;
G = zeros(1, N+1);
z0 = -0.5;
z = 0.9;
syms F(t)
F(t) = (t^(N + 1) - 1) / (t - 1);
G(1) = F(z0);
for j = 1:N
    F = diff(F) / j;
    G(j+1) = + double(F(z0)) * (z - z0)^j;
end
sum(G(1:153))
%The answer is about 5.7238e+25, the value blows up not agree with the
%solution.


%c)
% to calculate G(1.1)
N = i - 1;
G = zeros(1, N+1);
z0 = -0.5;
z = -1.1;
syms F(t)
F(t) = (t^(N + 1) - 1) / (t - 1);
G(1) = F(z0);

for j = 1:N
    F = diff(F) / j;
    G(j+1) = double(F(z0)) * (z - z0)^j;
end
sum(G(1:153))
% The answer is about 1.0253e+6, which is not close to the correct answer 0.47619.

% d) Based on the data, we see that when m is small, G(z) is sufficiently
% close to the solution. However, due to the precision is not met, the
% value I got didn't converge at N=152. Better try other method to
% implement the data.