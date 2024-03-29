function D = D5(in1,in2)
%D5
%    D = D5(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Jun-2019 23:28:07

J1 = in1(:,11);
J2 = in1(:,12);
J3 = in1(:,13);
J4 = in1(:,14);
J5 = in1(:,15);
L1 = in1(:,6);
L2 = in1(:,7);
L3 = in1(:,8);
L4 = in1(:,9);
M1 = in1(:,1);
M2 = in1(:,2);
M3 = in1(:,3);
M4 = in1(:,4);
M5 = in1(:,5);
Mr_com_x1 = in1(:,21);
Mr_com_x2 = in1(:,22);
Mr_com_y1 = in1(:,26);
Mr_com_x3 = in1(:,23);
Mr_com_y2 = in1(:,27);
Mr_com_x4 = in1(:,24);
Mr_com_y3 = in1(:,28);
Mr_com_x5 = in1(:,25);
Mr_com_y4 = in1(:,29);
Mr_com_y5 = in1(:,30);
qR2 = in2(:,2);
qR3 = in2(:,3);
qR4 = in2(:,4);
qR5 = in2(:,5);
t2 = cos(qR2);
t3 = cos(qR3);
t4 = cos(qR4);
t5 = cos(qR5);
t6 = sin(qR2);
t7 = sin(qR3);
t8 = sin(qR4);
t9 = sin(qR5);
t10 = qR2+qR3;
t11 = qR3+qR4;
t12 = qR4+qR5;
t13 = L1.^2;
t14 = L2.^2;
t15 = L3.^2;
t16 = L4.^2;
t17 = Mr_com_x2.^2;
t18 = Mr_com_x3.^2;
t19 = Mr_com_y2.^2;
t20 = Mr_com_x4.^2;
t21 = Mr_com_y3.^2;
t22 = Mr_com_x5.^2;
t23 = Mr_com_y4.^2;
t24 = Mr_com_y5.^2;
t25 = cos(t10);
t26 = cos(t11);
t27 = cos(t12);
t28 = qR4+t10;
t29 = qR5+t11;
t30 = sin(t10);
t31 = sin(t11);
t32 = sin(t12);
t35 = t10+t12;
t38 = M3.*t14;
t39 = M4.*t14;
t40 = M5.*t14;
t41 = M4.*t15;
t42 = M5.*t15;
t43 = M5.*t16;
t44 = M2.*t17;
t45 = M2.*t19;
t46 = M3.*t18;
t47 = M3.*t21;
t48 = M4.*t20;
t49 = M4.*t23;
t50 = M5.*t22;
t51 = M5.*t24;
t52 = L1.*L2.*M3.*t2;
t53 = L1.*L2.*M4.*t2;
t54 = L1.*L2.*M5.*t2;
t55 = L2.*L3.*M4.*t3;
t56 = L2.*L3.*M5.*t3;
t57 = L3.*L4.*M5.*t4;
t58 = L1.*M2.*Mr_com_x2.*t2;
t59 = L2.*M3.*Mr_com_x3.*t3;
t60 = L3.*M4.*Mr_com_x4.*t4;
t61 = L4.*M5.*Mr_com_x5.*t5;
t62 = L1.*M2.*Mr_com_y2.*t6;
t63 = L2.*M3.*Mr_com_y3.*t7;
t64 = L3.*M4.*Mr_com_y4.*t8;
t65 = L4.*M5.*Mr_com_y5.*t9;
t33 = cos(t28);
t34 = cos(t29);
t36 = sin(t28);
t37 = sin(t29);
t66 = cos(t35);
t67 = sin(t35);
t68 = t55.*2.0;
t69 = t56.*2.0;
t70 = t57.*2.0;
t71 = t59.*2.0;
t72 = t60.*2.0;
t73 = t61.*2.0;
t74 = t63.*2.0;
t75 = t64.*2.0;
t76 = t65.*2.0;
t77 = -t62;
t78 = -t63;
t80 = -t64;
t82 = -t65;
t84 = L1.*L3.*M4.*t25;
t85 = L1.*L3.*M5.*t25;
t86 = L2.*L4.*M5.*t26;
t87 = L1.*M3.*Mr_com_x3.*t25;
t88 = L2.*M4.*Mr_com_x4.*t26;
t89 = L3.*M5.*Mr_com_x5.*t27;
t90 = L1.*M3.*Mr_com_y3.*t30;
t91 = L2.*M4.*Mr_com_y4.*t31;
t92 = L3.*M5.*Mr_com_y5.*t32;
t79 = -t74;
t81 = -t75;
t83 = -t76;
t93 = L1.*L4.*M5.*t33;
t94 = L1.*M4.*Mr_com_x4.*t33;
t95 = L2.*M5.*Mr_com_x5.*t34;
t96 = L1.*M4.*Mr_com_y4.*t36;
t97 = L2.*M5.*Mr_com_y5.*t37;
t98 = t86.*2.0;
t99 = t88.*2.0;
t100 = t89.*2.0;
t101 = t91.*2.0;
t102 = t92.*2.0;
t105 = L1.*M5.*Mr_com_x5.*t66;
t106 = L1.*M5.*Mr_com_y5.*t67;
t107 = -t90;
t108 = -t91;
t110 = -t92;
t116 = J5+t50+t51+t61+t82;
t103 = t95.*2.0;
t104 = t97.*2.0;
t109 = -t101;
t111 = -t102;
t112 = -t96;
t113 = -t97;
t115 = -t106;
t117 = t89+t110+t116;
t119 = J4+J5+t43+t48+t49+t50+t51+t57+t60+t73+t80+t83+t89+t110;
t114 = -t104;
t118 = t95+t113+t117;
t121 = t86+t88+t95+t108+t113+t119;
t122 = J3+J4+J5+t41+t42+t43+t46+t47+t48+t49+t50+t51+t55+t56+t59+t70+t72+t73+t78+t81+t83+t86+t88+t95+t100+t108+t111+t113;
t120 = t105+t115+t118;
t123 = t93+t94+t105+t112+t115+t121;
t124 = t84+t85+t87+t93+t94+t105+t107+t112+t115+t122;
t125 = J2+J3+J4+J5+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t58+t68+t69+t70+t71+t72+t73+t77+t79+t81+t83+t84+t85+t87+t93+t94+t98+t99+t100+t103+t105+t107+t109+t111+t112+t114+t115;
D = reshape([J1+J2+J3+J4+J5+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52.*2.0+t53.*2.0+t54.*2.0+t58.*2.0-t62.*2.0+t68+t69+t70+t71+t72+t73+t79+t81+t83+t84.*2.0+t85.*2.0+t87.*2.0-t90.*2.0+t93.*2.0+t94.*2.0-t96.*2.0+t98+t99+t100+t103+t105.*2.0-t106.*2.0+t109+t111+t114+M2.*t13+M3.*t13+M4.*t13+M5.*t13+M1.*Mr_com_x1.^2+M1.*Mr_com_y1.^2,t125,t124,t123,t120,t125,J2+J3+J4+J5+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t68+t69+t70+t71+t72+t73+t79+t81+t83+t98+t99+t100+t103+t109+t111+t114,t122,t121,t118,t124,t122,J3+J4+J5+t41+t42+t43+t46+t47+t48+t49+t50+t51+t70+t72+t73+t81+t83+t100+t111,t119,t117,t123,t121,t119,J4+J5+t43+t48+t49+t50+t51+t73+t83,t116,t120,t118,t117,t116,J5+t50+t51],[5,5]);
