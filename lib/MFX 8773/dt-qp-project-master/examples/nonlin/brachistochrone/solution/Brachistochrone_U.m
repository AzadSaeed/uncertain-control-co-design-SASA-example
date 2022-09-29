function out1 = Brachistochrone_U(f,t,t1,t2,theta,w1,w2)
%BRACHISTOCHRONE_U
%    OUT1 = BRACHISTOCHRONE_U(F,T,T1,T2,THETA,W1,W2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    21-Jun-2020 00:46:07

t4 = -t;
out1 = heaviside(t1+t4).*(pi./2.0+t4.*w1)+w2.*heaviside(t-t2).*(f+t4)+theta.*heaviside(t-t1).*heaviside(t2+t4);
