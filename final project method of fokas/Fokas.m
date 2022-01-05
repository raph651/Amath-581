%% Final Project 
%Rapahel Liu
%% Suppose given a heat equation u_t = γ*u_xx + β*u_x + α*u 
% 0<x<L, ,L=30, n=128,γ=0.01, β=0.1, α=0;
%given the initial condition u(x,0)= u_0(x)= exp(-(x-16).^2/10)
%and the boundary condition u(0,t)=0
%% Method of Fokas
L=30;
tspan=0:1:10;
gamma=0.01; beta=0.1; alpha=0;

% the path representing partial Omega is k_R^2=k_I^2-beta/gamma*k_I
% use re denote k_R, im denote k_I, then the partial Omega path is 
% described as below: 
num_path_points=61;
re=linspace(-L,L,num_path_points);
im=(re.^2+beta^2/(4*gamma^2)).^0.5+beta/(2*gamma);
path=re+1i.*im;

%W(k)=gamma*k^2-1i*beta*k-alpha
%v(k)=i*beta/gamma-k;
W= @(k) gamma.*k.^2-1i*beta.*k-alpha;
v=@(k) 1i*beta/gamma-k;

I=@(k,x) exp(-1i.*k.*x).*exp(-(x-16).^2/10);
I_v=@(k,x) exp(-1i.*v(k).*x).*exp(-(x-16).^2/10);

u_k=@(k) integral(@(x) I(k,x),0,L);
u_v=@(k) integral(@(x) I_v(k,x),0,L);

n1=128;
x1=linspace(0,L,n1);

u_x=zeros(length(tspan),n1);
for t=0:length(tspan)-1
I1=@(x) integral(@(k) 1/(2*pi).*exp(1i.*k.*x-W(k).*t).*u_k(k),-L,L,'ArrayValued',true,'RelTol',1e-2,'AbsTol',1e-3);
I2=@(x) integral(@(k) -1/(2*pi).*exp(1i.*k.*x-W(k).*t).*u_v(k),path(1),path(num_path_points),"Waypoints",path,"ArrayValued",true,"RelTol",1e-2,"AbsTol",1e-3);
u_x_T=@(x) I1(x)+I2(x);
for i=1:n1
tic
u_x(t+1,i)=u_x_T(x1(i));
toc
end
end





%% Spectral Method ()
% ut=Lu, where L = γ*d^2x+β*dx+α
% ut_hat=(-γk^2+iβk+α)u_hat, where L_hat= -γk^2+iβk+α , in Fourier domain. 

L=30; 
n2=256;
x2=linspace(-L,L,n2+1);x2=x2(1:n2);
gamma=0.01; beta=0.1; alpha=0;

u_init=exp(-(x2-16).^2/10);
kx=2*pi/L*[0:n2/2-1 -n2/2:-1]';kx(1)=1e-6;
uf_init=fft(u_init);
heat_spec= @(t,u) (-gamma*kx.^2+1i*kx.*beta+alpha).*u;

[t,uf]=ode45(heat_spec,tspan,uf_init);



%% Time Stepping method()

L=30; n=128;
x=linspace(0,L,n+1);x=x(1:n);
x3=x;
gamma=0.01; beta=0.1; alpha=0;
dx=x(2)-x(1);
tspan=0:1:10;

u_init=exp(-(x-16).^2/10);

e1=ones(n,1);
Dxx=spdiags([-e1 16*e1 -30*e1 16*e1 -e1],[2 1 0 -1 -2],n,n);
Dxx(1,n)=16;Dxx(1,n-1)=-1;Dxx(2,n)=-1;Dxx(n-1,1)=-1;Dxx(n,1)=16;Dxx(n,2)=-1;
Dxx=Dxx./(12*dx^2);
Dx=spdiags([-e1 8*e1 -8*e1 e1],[2 1 -1 -2],n,n);
Dx(1,n)=-8;Dx(1,n-1)=1;Dx(2,n)=1;Dx(n-1,1)=-1;Dx(n,1)=8;Dx(n,2)=-1;
Dx=Dx./(12*dx);
I=eye(n);

heat_tstep=@(t,u) (gamma.*Dxx+beta.*Dx+alpha.*I)*u;

[t,u]=ode45(heat_tstep,tspan,u_init);

%% Visulization
v=VideoWriter('solutions to heat equation.avi');
v.FrameRate=5;
open(v);

for i=1:length(tspan)
y1=u_x(i,:);
y2=ifft(uf(i,:));
y3=u(i,:);
figure(1)
plot(x1,real(y1),'b--',x2,real(y2),x3,real(y3),'rx');
axis([0 L 0 1]);
xlabel('x')
ylabel('u(x) at t')
legend('fokas','spec','tstep')
drawnow
frame =getframe(gcf);
writeVideo(v,frame);
pause(1)
end