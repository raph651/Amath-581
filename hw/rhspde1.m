function rhs = rhspde1(t,w,v,A,B,C)

phi=A\w;

rhs=v*A*w-((B*phi).*(C*w)-(C*phi).*(B*w));

end
