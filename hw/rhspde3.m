function rhs = rhspde3(t,w,v,A,B,C)

phi=bicgstab(A,w,1e-6);

rhs=v*A*w-((B*phi).*(C*w)-(C*phi).*(B*w));

end
