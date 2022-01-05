%% AMATH581 Raphael Liu

% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

% your solution code goes here
%% Exercise 1
A = [34 45;17 6];
A1 = A;
%% Exercise 2
A = [1 2;-1 1 ];
B = [2 0;0 2];
C = [2 0 -3;0 0 -1];
D = [1 2;2 3;-1 0];
x = [1 ;0];
y = [0 ;1];
z = [1; 2; -1];


A2 = A+B;

A3 = 3*x -4*y;

A4 = A*x;

A5 = B*(x-y);

A6 = D*x;

A7 = D*y + z;

A8 = A*B;

A9 = B*C;

A10 = C*D;

%% Exercise 3
% part 1 Newton-Raphson method
x=[];
error=[];

value=-0.73908513321516064166;
x(1)=-3;
error(1)=abs(x(1)-value);

i=1;
while( error(i) > 1e-6)

    x(i+1)=x(i)-(-x(i)-cos(x(i)))/(-1+sin(x(i)));
    error(i+1)=abs(x(i+1)-value);
    i=i+1;

end



% part 2 Midpoint method
y=[];
error=[];
value=-0.73908513321516064166;
y(1)=-3;
a=-3;
b=1;
error(1)=abs(y(1)-value);

i=1;
while(error(i)>1e-6)
    c=(a+b)/2;

    if (-a-cos(a))*(-c-cos(c)) <0
        b=c;
    else
        a=c;
    end

    y(i+1)=c;
    error(i+1)=abs(y(i+1)-value);
    i=i+1;
end

A11 = x';
A12 = y(2:end)';
A13 = [5 22];

% your extra functions, if you need them, can be in additional files (don't forget to upload them too!)