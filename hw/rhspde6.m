function rhs=rhspde6(t,w,v,A,B,C)
L=20;
n=64;
b2=reshape(w,[n n]);
b2=fft2(b2);

k=2*pi/L*[0: (n/2-1) (-n/2):-1];
k(1)=1e-6;
[K1,K2]=meshgrid(k,k);
K=K1.^2+K2.^2;
phi=-b2./K;

phi=real(ifft2(phi));
phi=reshape(phi,[n^2,1]);

rhs=v.*A*w-((B*phi).*(C*w)-(C*phi).*(B*w));


end