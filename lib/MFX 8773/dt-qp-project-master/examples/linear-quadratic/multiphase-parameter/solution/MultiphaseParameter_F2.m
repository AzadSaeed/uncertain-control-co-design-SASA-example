function F2 = MultiphaseParameter_F2(P,tf1,tf2,in4)
%MULTIPHASEPARAMETER_F2
%    F2 = MULTIPHASEPARAMETER_F2(P,TF1,TF2,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    10-May-2018 23:44:58

y01 = in4(1,:);
y02 = in4(2,:);
y03 = in4(3,:);
y04 = in4(4,:);
y05 = in4(5,:);
y06 = in4(6,:);
t2 = P.^2;
t3 = tf1.*(5.0e-2-3.0i);
t4 = tf2.*(-5.0e-2+3.0i);
t5 = t3+t4;
t6 = exp(t5);
t7 = tf1.*(5.0e-2+3.0i);
t8 = tf2.*(-5.0e-2-3.0i);
t9 = t7+t8;
t10 = exp(t9);
t11 = tf1.*(5.0e-2-1i);
t12 = tf2.*(-5.0e-2+1i);
t13 = t11+t12;
t14 = exp(t13);
t15 = tf1.*(5.0e-2+1i);
t16 = tf2.*(-5.0e-2-1i);
t17 = t15+t16;
t18 = exp(t17);
t19 = tf1.*(5.0e-2-2.0i);
t20 = tf2.*(-5.0e-2+2.0i);
t21 = t19+t20;
t22 = exp(t21);
t23 = tf1.*(5.0e-2+2.0i);
t24 = tf2.*(-5.0e-2-2.0i);
t25 = t23+t24;
t26 = exp(t25);
t27 = y03.^2;
t28 = y06.^2;
t29 = tf1.*(1.0./1.0e1);
t34 = tf2.*(1.0./1.0e1);
t30 = t29-t34;
t31 = exp(t30);
t32 = y02.^2;
t33 = y05.^2;
t35 = y01.^2;
t36 = y04.^2;
F2 = t2.*2.446245457693481e1+t27.*1.0e1+t28.*1.0e1+t32.*1.0e1+t33.*1.0e1+t35.*1.0e1+t36.*1.0e1-t31.*(t2.*1.4404e7+t27.*6.4836005e7+t28.*6.4836005e7+P.*y03.*4.24918e7-P.*y06.*4.39322e7).*1.542352894815157e-7+t6.*(t2.*(4.87466814773674e5+4.560088864204388e5i)+P.*y03.*(2.36e4+1.416e6i)+P.*y06.*(-1.416e6+2.36e4i)).*7.711764474075786e-8+t10.*(t2.*(4.87466814773674e5-4.560088864204388e5i)+P.*y03.*(2.36e4-1.416e6i)+P.*y06.*(-1.416e6-2.36e4i)).*7.711764474075786e-8-P.*y01.*1.895261845386534e1+P.*y02.*9.993753903810119+P.*y03.*6.775895584559844-P.*y04.*2.094763092269327e1+P.*y05.*2.49843847595253e-1-P.*y06.*6.553735073590669-t6.*(t2.*(4.71466814773674e5-5.039911135795612e5i)+P.*y03.*(1.464e6-2.44e4i)+P.*y06.*(2.44e4+1.464e6i)).*7.711764474075786e-8-t10.*(t2.*(4.71466814773674e5+5.039911135795612e5i)+P.*y03.*(1.464e6+2.44e4i)+P.*y06.*(2.44e4-1.464e6i)).*7.711764474075786e-8-t22.*(t2.*(3.117198474133004e-3-1.559573970798251e-4i)+P.*y02.*(6.242194818119999e-3-1.56054870453e-4i)+P.*y05.*(1.56054870453e-4+6.242194818119999e-3i))-t26.*(t2.*(3.117198474133004e-3+1.559573970798251e-4i)+P.*y02.*(6.242194818119999e-3+1.56054870453e-4i)+P.*y05.*(1.56054870453e-4-6.242194818119999e-3i))-t2.*tf1.*2.467016827392358+t2.*tf2.*2.467016827392358-t14.*(t2.*(-1.664039900249377e5-1.360798004987531e5i)+P.*y01.*(1.52e5-7.6e3i)+P.*y04.*(7.6e3+1.52e5i)).*6.218866798092052e-6-t18.*(t2.*(-1.664039900249377e5+1.360798004987531e5i)+P.*y01.*(1.52e5+7.6e3i)+P.*y04.*(7.6e3-1.52e5i)).*6.218866798092052e-6-t14.*(t2.*(1.504039900249377e5-1.839201995012469e5i)+P.*y01.*(8.4e3+1.68e5i)+P.*y04.*(-1.68e5+8.4e3i)).*6.218866798092052e-6-t18.*(t2.*(1.504039900249377e5+1.839201995012469e5i)+P.*y01.*(8.4e3-1.68e5i)+P.*y04.*(-1.68e5-8.4e3i)).*6.218866798092052e-6+t22.*(t2.*(6.238295883193004e-3+1.246879389653202e-1i)+P.*y02.*(6.242194818119999e-3+2.496877927248e-1i)+P.*y05.*(-2.496877927248e-1+6.242194818119999e-3i))+t26.*(t2.*(6.238295883193004e-3-1.246879389653202e-1i)+P.*y02.*(6.242194818119999e-3-2.496877927248e-1i)+P.*y05.*(-2.496877927248e-1-6.242194818119999e-3i))-t31.*(t2.*1.249219237976265+t32.*5.0+t33.*5.0+P.*y02.*4.996876951905059-P.*y05.*1.249219237976265e-1).*2.0-t31.*(t2.*1.604e6+t35.*8.04005e5+t36.*8.04005e5-P.*y01.*1.6842e6-P.*y04.*1.5238e6).*1.24377335961841e-5;
