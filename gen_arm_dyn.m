% Generate 2n-DOF planar dynamic model using Lagrangian method
% Author: Vedant
% Based on work of Prof. Hae-Won Park

%function [B,C,D,G] = gen_arm_dyn(n)

clear
%close all
clc
n=3;
qR = sym('qR', [1,n], 'real'); % Joints of Right array
dqR = sym('dqR', [1,n], 'real'); %Angular rates
ddqR = sym('ddqR', [1,n], 'real'); %Angular rates
M = sym('M', [1,n], 'real'); %mass of each link
L = sym('L', [1,n], 'real'); %length of each link
J = sym('J', [1,n], 'real'); %inertia of each link
K = sym('K', [1,n], 'real'); %inertia of each link
MR_COM_x = sym('Mr_com_x', [1,n], 'real'); % COM position for each link
MR_COM_y = sym('Mr_com_y', [1,n], 'real');

% % For Center Body
syms LCenter JCenter MCenter real % center body mass, intertia and length
% syms MB_COM_x MB_COM_y real
% syms qCenter dqCenter ddqCenter real
% syms xCenter dxCenter ddxCenter real
% syms yCenter dyCenter ddyCenter real

% if(stance_flag_fore)
q = [qR];
dq = [dqR];

Design_Parameters = [M,L,J,K,MR_COM_x,MR_COM_y,LCenter];
Control_States = [q,dq];
% end
disp('Initialization done')

%% forward kinematics
rot_2d = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)]; %counterclockwise rotation matrix for [x;y] vectors

abs_qR = sym('abs_qR',[1 n]);

abs_dqR = sym('abs_dqR',[1 n]);

for(i=1:n)
    %absolute angle
    abs_qR(i) = sum(qR(1:i));
    abs_dqR(i) = sum(dqR(1:i));
end

pFR = sym('pFR',[2 n+1]);
pFR(:,1) = [LCenter;0];
for(i=2:n+1)
    pFR(:,i) = [(pFR(1,(i-1)));(pFR(2,(i-1)))]+rot_2d(abs_qR(i-1))*[L(i-1);0];
end

% calculate COM position for each link
pFR_COM = sym('pFR_COM',[2 n]);

for(i=1:n)
    pFR_COM(:,i) = [LCenter;0]+rot_2d(abs_qR(i))*[ MR_COM_x(i);MR_COM_y(i)];
    if(i>1)
        pFR_COM(:,i) = pFR_COM(:,i)+pFR(:,(i-1));
    end
end
disp('Sat frame done')

%% Velocity kinematics

for(i=1:n)
    VR_COM(:,i)=jacobian(pFR_COM(:,i),q)*dq';
    
end

%% Kinetic Energy, Potential Energy, and Lagrangiandq2abs

for(i=1:n)
    KE_R(i) = 0.5*M(i)*VR_COM(:,i).'*VR_COM(:,i)+0.5*J(i)*abs_dqR(i)^2;
    if i == 1
    PE_R(i) = 0.5*K(i)*(qR(i)^2);
    else
        PE_R(i) = 0.5*K(i)*((qR(i-1)-qR(i))^2);
   end
end

KE = sum(KE_R);
PE = sum(PE_R);


Upsilon = [qR] ; % where control torques go
dimensionparameters = [L,LCenter];
[D, C, G, B] = std_dynamics(KE,PE,q,dq, Upsilon);

matlabFunction(B,'File','B','Vars',{Design_Parameters ,Control_States})
matlabFunction(C,'File','C','Vars',{Design_Parameters ,Control_States})
matlabFunction(D,'File','D','Vars',{Design_Parameters ,Control_States})
matlabFunction(G,'File','G','Vars',{Design_Parameters ,Control_States})
matlabFunction(Jacobian,'File','Jacobian','Vars',{qR,dimensionparameters})
matlabFunction(pFR,'File','visualize','Vars',{qR,dimensionparameters})
%% simulation
clear D C B G

qRn = 0*0.25*ones(1,n); % Joints of Right array
dqRn = -0.0*ones(1,n); %Angular rates
Mn = 2*ones(1,n); %mass of each link
Ln = ones(1,n); %length of each link
Jn = 1*ones(1,n); %inertia of each link
Kn = 0*ones(1,n); %inertia of each link
MR_COM_xn = 0.5*ones(1,n); % COM position for each link
MR_COM_yn = 0*ones(1,n);
K = 1
% % For Center Body
LCentern = 1;
% q = [qR];
% dq = [dqR];

dimensionparametersn = [Ln,LCentern];
Design_Parametersn = [Mn,Ln,Jn,Kn,MR_COM_xn,MR_COM_yn,LCentern];
init = transpose([qRn,dqRn]);

desired_omega = [2*pi*0.5;4.694*2*pi*0.5/1.875];%;;7.885*2*pi*0.5/1.875

Input.omega = desired_omega;
Input.DParams = Design_Parametersn;
Input.x = init;
Input.n = n;


clc

K_list = [] ;
Loss_list = [];
L_list = []
opts = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','FunctionTolerance', 1e-12,'StepTolerance', 1e-12,'FunctionTolerance',1e-12,'MaxFunctionEvaluations',1000);

for i = 0:0.1:pi/2
    qRn = linspace(0,i,n+1);
    qRn = qRn(2:end);
    init = transpose([qRn,dqRn]);
    Input.x = init;
    
%     [K,Loss] = fmincon(@(x) findK(x, Input), ones(1),[],[]);
    [KL,Loss] = fsolve(@(x) findKL(x, Input), [ones(1), Ln, LCentern],opts);
    K
    Loss_list = [Loss_list Loss];
    K_list = [K_list KL(1)];
    L_list = [L_list; KL(2:end)]
%     M1 = D(Input.DParams,transpose(Input.x));
%     G = diag(Kn);
%     A = M1\G
%     eigA = sort(eig(A));
    % F = eigA(2) - desired_omega^2;
%     F = [eigA(1); eigA(2) - desired_omega.^2]
end

Kn = K*ones(1,n); %inertia of each link

%Kn =    1.0 * [   20.3989,17.9315]
qRn = linspace(0,0.25,n+1);
    qRn = qRn(2:end);
    init = transpose([qRn,dqRn]);
    
    %Design_Parametersn = [Mn,L_list(end,1:end-1),Jn,Kn,MR_COM_xn,MR_COM_yn,L_list(end,end)];
% Design_Parametersn = [Mn,Ln,Jn,Kn,MR_COM_xn,MR_COM_yn,LCentern];
[t,xt] = ode113(@(t,x)dynamics(t,x,n,Design_Parametersn),[0 10],init);


h = animatedline('Marker','o');
axis([-0 6 -5 5])
axis square
dt = diff(t);
for (i=1:length(xt))
    p = visualize(xt(i,1:n),dimensionparametersn);
    x = [0 p(1,:)];
    y = [0 p(2,:)];
    clearpoints(h);
    addpoints(h,x,y);
    if(i ~= length(xt))
        pause(dt(i));
    end
end