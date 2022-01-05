function rhs = rhspde4(t,w,v,A,B,C)

phi=gmres(A,w,[],[],200);

rhs=v*A*w-((B*phi).*(C*w)-(C*phi).*(B*w));

end
