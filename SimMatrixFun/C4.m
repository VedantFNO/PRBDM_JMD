function C = C4(in1,in2)
%C4
%    C = C4(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Jun-2019 23:05:55

L1 = in1(:,5);
L2 = in1(:,6);
L3 = in1(:,7);
M2 = in1(:,2);
M3 = in1(:,3);
M4 = in1(:,4);
Mr_com_x2 = in1(:,18);
Mr_com_x3 = in1(:,19);
Mr_com_y2 = in1(:,22);
Mr_com_x4 = in1(:,20);
Mr_com_y3 = in1(:,23);
Mr_com_y4 = in1(:,24);
dqR1 = in2(:,5);
dqR2 = in2(:,6);
dqR3 = in2(:,7);
dqR4 = in2(:,8);
qR2 = in2(:,2);
qR3 = in2(:,3);
qR4 = in2(:,4);
t2 = cos(qR2);
t3 = cos(qR3);
t4 = cos(qR4);
t5 = sin(qR2);
t6 = sin(qR3);
t7 = sin(qR4);
t8 = qR2+qR3;
t9 = qR3+qR4;
t28 = dqR1+dqR2+dqR3+dqR4;
t10 = Mr_com_y4.*t4;
t11 = Mr_com_x4.*t7;
t12 = cos(t8);
t13 = cos(t9);
t14 = qR4+t8;
t15 = sin(t8);
t16 = sin(t9);
t20 = M2.*Mr_com_y2.*t2;
t21 = M3.*Mr_com_y3.*t3;
t22 = L2.*M3.*t5;
t23 = L2.*M4.*t5;
t24 = L3.*M4.*t6;
t26 = M2.*Mr_com_x2.*t5;
t27 = M3.*Mr_com_x3.*t6;
t17 = cos(t14);
t18 = sin(t14);
t19 = L3.*t10;
t25 = L3.*t11;
t31 = L2.*dqR1.*t21;
t32 = L2.*dqR2.*t21;
t33 = L2.*dqR3.*t21;
t35 = L2.*dqR1.*t24;
t36 = L2.*dqR2.*t24;
t37 = L2.*dqR3.*t24;
t38 = L2.*dqR1.*t27;
t39 = L2.*dqR2.*t27;
t40 = L2.*dqR3.*t27;
t42 = L2.*Mr_com_y4.*t13;
t43 = M3.*Mr_com_y3.*t12;
t44 = M4.*Mr_com_y4.*t13;
t45 = L3.*M4.*t15;
t46 = L2.*Mr_com_x4.*t16;
t47 = M3.*Mr_com_x3.*t15;
t48 = M4.*Mr_com_x4.*t16;
t74 = t10+t11;
t29 = M4.*t19;
t30 = M4.*t25;
t49 = L1.*Mr_com_y4.*t17;
t50 = M4.*Mr_com_y4.*t17;
t51 = L1.*Mr_com_x4.*t18;
t52 = M4.*Mr_com_x4.*t18;
t53 = M4.*t42;
t54 = M4.*t46;
t55 = L1.*dqR1.*t43;
t56 = L1.*dqR2.*t43;
t57 = L1.*dqR3.*t43;
t63 = L1.*dqR1.*t45;
t64 = L1.*dqR2.*t45;
t65 = L1.*dqR3.*t45;
t66 = L1.*dqR1.*t47;
t67 = L1.*dqR2.*t47;
t68 = L1.*dqR3.*t47;
t75 = -t33;
t77 = -t37;
t78 = -t40;
t109 = L3.*M4.*dqR4.*t74;
t111 = t19+t25+t42+t46;
t113 = t21+t24+t27+t44+t48;
t34 = dqR4.*t29;
t41 = dqR4.*t30;
t58 = dqR1.*t53;
t59 = dqR2.*t53;
t60 = dqR3.*t53;
t61 = dqR4.*t53;
t62 = M4.*t49;
t69 = dqR1.*t54;
t70 = dqR2.*t54;
t71 = dqR3.*t54;
t72 = dqR4.*t54;
t73 = M4.*t51;
t88 = -t55;
t89 = -t56;
t90 = -t57;
t93 = -t63;
t94 = -t64;
t95 = -t65;
t96 = -t66;
t97 = -t67;
t98 = -t68;
t110 = -t109;
t112 = M4.*dqR4.*t111;
t115 = L2.*dqR1.*t113;
t116 = L2.*dqR2.*t113;
t117 = L2.*dqR3.*t113;
t119 = t49+t51+t111;
t120 = t20+t22+t23+t26+t43+t45+t47+t50+t52;
t76 = -t34;
t79 = -t41;
t80 = dqR1.*t62;
t81 = dqR2.*t62;
t82 = dqR3.*t62;
t83 = dqR4.*t62;
t84 = dqR1.*t73;
t85 = dqR2.*t73;
t86 = dqR3.*t73;
t87 = dqR4.*t73;
t91 = -t60;
t92 = -t61;
t99 = -t71;
t100 = -t72;
t114 = -t112;
t118 = -t117;
t101 = -t80;
t102 = -t81;
t103 = -t82;
t104 = -t83;
t105 = -t84;
t106 = -t85;
t107 = -t86;
t108 = -t87;
C = reshape([-dqR3.*(t53+t54+t62+t73+L2.*t21+L2.*t24+L2.*t27+L1.*t43+L1.*t45+L1.*t47)-L1.*dqR2.*t120-M4.*dqR4.*t119,t114+t118+L1.*dqR1.*t120,t31+t32+t35+t36+t38+t39+t55+t58+t59+t63+t66+t69+t70+t76+t79+t80+t84,dqR3.*(t29+t30)+dqR1.*(t29+t30+t53+t54+t62+t73)+dqR2.*(t29+t30+t53+t54),t75+t76+t77+t78+t79+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98+t99+t100+t101+t102+t103+t104+t105+t106+t107+t108-L1.*dqR1.*t20-L1.*dqR2.*t20-L1.*dqR1.*t22-L1.*dqR1.*t23-L1.*dqR2.*t22-L1.*dqR2.*t23-L1.*dqR1.*t26-L1.*dqR2.*t26,t114+t118,t110+t115+t116,M4.*dqR1.*t111+M4.*dqR2.*t111+L3.*M4.*dqR3.*t74,-t31-t32-t35-t36-t38-t39-t58-t59-t69-t70+t75+t76+t77+t78+t79+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98+t99+t100+t101+t102+t103+t104+t105+t106+t107+t108,t114-t115-t116+t118,t110,L3.*M4.*t74.*(dqR1+dqR2+dqR3),-M4.*t28.*t119,-M4.*t28.*t111,-L3.*M4.*t28.*t74,0.0],[4,4]);
