function [c1 c2 c3]=Quartic_Surrogate(A,x_t)

c1=trace(A.'*A);

Q=(x_t.'*A*x_t)*A-c1*x_t*x_t.';
Q=Q*2;

[~,n]=eig(Q);

c2=max(diag(n));

c3=-2*x_t.'*(c2*eye(length(x_t))-Q);
c3=c3.';
