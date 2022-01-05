function rhs = rhspde5(t,w,v,A,B,C)

n=128;
L=10;
w_res=reshape(w,[n n]);
%
wt=fft2(w_res);

kx=pi/L*[0: (n/2-1) (-n/2):-1];
ky=pi/L*[0: (n/2-1) (-n/2):-1];
kx(1) = 1e-6;
ky(1) = 1e-6;
[KX,KY]=meshgrid(kx,kx);
K=KX.^2+KY.^2;

phi=real(ifft2(-wt./K));
phi=reshape(phi,[n^2 1]);



rhs=v.*A*w-((B*phi).*(C*w)-(C*phi).*(B*w));

end
