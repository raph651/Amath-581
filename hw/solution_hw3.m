
%% Homework 3 Raphael Liu
% sparse matrix for dx, dy and Laplacian using spdiags() method.
%% Problem 1
% a)
xspan = [-10, 10];
yspan = [-10, 10];
n = 8;

x = linspace(xspan(1), xspan(2), n+1);
y = linspace(yspan(1), yspan(2), n+1);

h = x(2) - x(1);
I_array = ones(n, 1);

B = spdiags([I_array, -4 * I_array, I_array, I_array, I_array], [-1, 0, 1, n - 1, 1 - n], n, n);

I = diag(ones(n, 1));

%A_Laplacian=spdiags([full(B) I],[0 -1],n^2,n^2)
Z = zeros(n);
A_Laplacian = zeros(n^2);

A_Laplacian = [B, I, Z, Z, Z, Z, Z, I; ...
    I, B, I, Z, Z, Z, Z, Z; ...
    Z, I, B, I, Z, Z, Z, Z; ...
    Z, Z, I, B, I, Z, Z, Z; ...
    Z, Z, Z, I, B, I, Z, Z; ...
    Z, Z, Z, Z, I, B, I, Z; ...
    Z, Z, Z, Z, Z, I, B, I; ...
    I, Z, Z, Z, Z, Z, I, B] ./ h^2;
dx_matrix =  [Z, I, Z, Z, Z, Z, Z, -I; ...
    -I, Z, I, Z, Z, Z, Z, Z; ...
    Z, -I, Z, I, Z, Z, Z, Z; ...
    Z, Z, -I, Z, I, Z, Z, Z; ...
    Z, Z, Z, -I, Z, I, Z, Z; ...
    Z, Z, Z, Z, -I, Z, I, Z; ...
    Z, Z, Z, Z, Z, -I, Z, I; ...
    I, Z, Z, Z, Z, Z, -I, Z] ./ (2 * h);
Y = spdiags([-I_array, I_array, -I_array, I_array], [-1, 1, n - 1, 1 - n], n, n);
dy_matrix = [Y, Z, Z, Z, Z, Z, Z, Z; ...
    Z, Y, Z, Z, Z, Z, Z, Z; ...
    Z, Z, Y, Z, Z, Z, Z, Z; ...
    Z, Z, Z, Y, Z, Z, Z, Z; ...
    Z, Z, Z, Z, Y, Z, Z, Z; ...
    Z, Z, Z, Z, Z, Y, Z, Z; ...
    Z, Z, Z, Z, Z, Z, Y, Z; ...
    Z, Z, Z, Z, Z, Z, Z, Y] ./ (2 * h);


A1 = full(A_Laplacian);
A2 = full(dx_matrix);
A3 = full(dy_matrix);

%%
% a)
load('permvec.mat');
load('Fmat.mat');
Fmatt = Fmat;

for row = 1:4
    for col = 1:4
        location = permvec((row - 1)*4+col);
        
location
        per_col = mod(location, 4);
        if per_col==0
            per_col=4;
        end
        per_row = (location - per_col) / 4 + 1;
       
        Fmat(141+col*20:160+col*20, 141+row*20:160+row*20) = Fmatt(141+per_col*20:160+per_col*20, 141+per_row*20:160+per_row*20);
 
    end
end


A4 = abs(Fmat);


Fmat = ifftshift(Fmat);


Fmat=ifft2(Fmat);
A5 = abs(Fmat);


b=uint8(A5) ;

figure()
set(gcf,"colormap",gray)
imshow(b)