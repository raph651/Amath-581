function rhs = rhspde2(t,w,v,A,B,C)

[L,U,P]=lu(A);
M=L\inv(P)*w;
phi=U\M;

rhs=v*A*w-((B*phi).*(C*w)-(C*phi).*(B*w));

end
