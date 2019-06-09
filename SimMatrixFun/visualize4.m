function pFR = visualize4(in1,in2)
%VISUALIZE4
%    PFR = VISUALIZE4(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Jun-2019 23:05:56

L1 = in2(:,1);
L2 = in2(:,2);
L3 = in2(:,3);
L4 = in2(:,4);
LCenter = in2(:,5);
qR1 = in1(:,1);
qR2 = in1(:,2);
qR3 = in1(:,3);
qR4 = in1(:,4);
t2 = cos(qR1);
t3 = sin(qR1);
t4 = qR1+qR2;
t5 = L1.*t2;
t6 = L1.*t3;
t7 = cos(t4);
t8 = qR3+t4;
t9 = sin(t4);
t10 = cos(t8);
t11 = qR4+t8;
t12 = sin(t8);
t13 = L2.*t7;
t14 = L2.*t9;
t15 = L3.*t10;
t16 = L3.*t12;
pFR = reshape([LCenter,0.0,LCenter+t5,t6,LCenter+t5+t13,t6+t14,LCenter+t5+t13+t15,t6+t14+t16,LCenter+t5+t13+t15+L4.*cos(t11),t6+t14+t16+L4.*sin(t11)],[2,5]);