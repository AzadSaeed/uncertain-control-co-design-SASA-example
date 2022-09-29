function out1 = BrysonHo166_Y(t,tf,v0,x0)
%BRYSONHO166_X
%    OUT1 = BRYSONHO166_X(T,TF,V0,X0)

%    This function was generated by the Symbolic Math Toolbox version 6.3.
%    22-Oct-2015 00:01:00

t2 = sin(tf);
t3 = cos(t);
t4 = sin(t);
t7 = tf.*2.0;
t5 = t-t7;
t6 = cos(t5);
t8 = tf.^2;
t9 = t2.^2;
t10 = t8-t9;
t11 = 1.0./t10;
t12 = sin(t5);
out1 = [t11.*(-t3.*x0+t6.*x0-t.*t3.*v0+t.*t6.*v0+t4.*t8.*v0.*2.0-t.*t4.*x0+t3.*t8.*x0.*2.0+t.*t12.*x0+t4.*tf.*x0.*2.0-t.*t4.*tf.*v0.*2.0-t.*t3.*tf.*x0.*2.0).*(1.0./2.0),t11.*(t3.*v0-t6.*v0-t.*t4.*v0-t3.*t8.*v0.*2.0+t.*t12.*v0+t4.*tf.*v0.*2.0+t.*t3.*x0-t.*t6.*x0+t4.*t8.*x0.*2.0+t.*t3.*tf.*v0.*2.0-t.*t4.*tf.*x0.*2.0).*(-1.0./2.0)];
