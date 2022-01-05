%% Homework 5 Raphael Liu

%% Problem 1
%   Ut = λ(A)U - ω(A)V + D1*∇^2U
%   Vt = ω(A)U + λ(A)V + D2*∇^2V
% where A^2 = U^2 + V^2 and ∇^2 = dx^2 + dy^2
% consider λ(A)= 1 - A^2;  ω(A) = -βA^2 ;
clear variables;
L=20; n=64;
x=linspace(-L/2,L/2,n+1); x=x(1:n);
y=x;

beta=1; D1=0.1; D2=0.1;
tspan=0:0.5:4;

kx=2*pi/L*[0:n/2-1 -n/2:-1];kx(1)=1e-6;
ky=kx;

[KX,KY]=meshgrid(kx,ky);
K=KX.^2+KY.^2;
Kvec=reshape(K,n^2,1);
m=1; % number of spirals

%initial conditions
[X,Y]=meshgrid(x,y);
u_init=tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
v_init=tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));

uf_init=fft2(u_init);
vf_init=fft2(v_init);
uf_init=reshape(uf_init,n^2,1);
vf_init=reshape(vf_init,n^2,1);
ff_init_vec=[uf_init;vf_init];

%%
%a)
[t,ff1]=ode45(@(t,f) rea_difu(t,f,Kvec,beta,D1,D2),tspan,ff_init_vec);
A1=real(ff1);
A2=imag(ff1);

%%

%b)
N=31;
[D,x]=cheb(N-1);

D_square=D^2;
D_square(1,:)=zeros(1,N);
D_square(N,:)=zeros(1,N);

I=eye(length(D_square));

Lap=(2/L)^2*(kron(D_square,I)+kron(I,D_square));

x=x.*L./2;
y=x;

[X,Y]=meshgrid(x,y);
u_init=tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
v_init=tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));

f_init_vec=[reshape(u_init,N^2,1);reshape(v_init,N^2,1)];
[t,f2]=ode45(@(t,f) cheb_cal(t,f,Lap,beta,D1,D2,x,y),tspan,f_init_vec);


A3=f2;

%%
function rhs=rea_difu(t,f,Kvec,beta,D1,D2)
n=64;

uf=f(1:n^2);
vf=f(n^2+1:end);


u=real(ifft2(reshape(uf,n,n)));
v=real(ifft2(reshape(vf,n,n)));

%pcolor(u);shading interp;
%drawnow;

lambda=1-(u.^2+v.^2);
omega=-beta*(u.^2+v.^2);

Lap_u=-uf.*Kvec;
Lap_v=-vf.*Kvec;


dt_u=reshape(fft2(lambda.*u-omega.*v),n^2,1)+D1*Lap_u;
dt_v=reshape(fft2(omega.*u+lambda.*v),n^2,1)+D2*Lap_v;

rhs=[dt_u;dt_v];



end

function rhs=cheb_cal(t,f,Lap,beta,D1,D2,x,y)
n=31;
L=20;
u_vec=f(1:n^2);
v_vec=f(n^2+1:end);

u=reshape(u_vec,n,n);
v=reshape(v_vec,n,n);

%pcolor(x,y,u);shading interp;
%drawnow;

A_square=u.^2+v.^2;


dt_u=reshape((1-A_square).*u+(beta*A_square).*v,n^2,1)+D1*Lap*u_vec;
dt_v=reshape((-beta*A_square).*u+(1-A_square).*v,n^2,1)+D2*Lap*v_vec;

rhs=[dt_u;dt_v];

end