% Generate 2n-DOF planar dynamic model using Lagrangian method
% Author: Vedant
% Based on work of Prof. Hae-Won Park

function [B,C,D,G,visualize] = gen_arm_dyn_fun(n)


%close all
clc

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

Design_Parameters = [M,L,J,K,MR_COM_x,MR_COM_y,LCenter,JCenter, MCenter];
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
    

    if(i>1)
        pFR_COM(:,i) = rot_2d(abs_qR(i))*[ MR_COM_x(i);MR_COM_y(i)]+pFR(:,(i));
    else
        pFR_COM(:,i) = [LCenter;0]+rot_2d(abs_qR(i))*[ MR_COM_x(i);MR_COM_y(i)];
    end
end
disp('Sat frame done')

%% Velocity kinematics

for(i=1:n)
    VR_COM(:,i)=jacobian(pFR_COM(:,i),q)*dq';
    
end
Jac = jacobian(pFR(:,end),q)
%% Kinetic Energy, Potential Energy, and Lagrangiandq2abs

for(i=1:n)
    KE_R(i) = 0.5*M(i)*VR_COM(:,i).'*VR_COM(:,i)+0.5*J(i)*abs_dqR(i)^2;
    PE_R(i) = 0.5*K(i)*(q(i)^2);
end

KE = sum(KE_R);
PE = sum(PE_R);


Upsilon = [q] ; % where control torques go
dimensionparameters = [L,LCenter];
[D, C, G, B] = std_dynamics(KE,PE,q,dq, Upsilon);
n_name = num2str(n)
cd SimMatrixFun
B = matlabFunction(B,'File',strcat('B',n_name),'Vars',{Design_Parameters ,Control_States})
C = matlabFunction(C,'File',strcat('C',n_name),'Vars',{Design_Parameters ,Control_States})
D = matlabFunction(D,'File',strcat('D',n_name),'Vars',{Design_Parameters ,Control_States})
G = matlabFunction(G,'File',strcat('G',n_name),'Vars',{Design_Parameters ,Control_States})
J = matlabFunction(Jac,'File',strcat('J',n_name),'Vars',{Design_Parameters ,Control_States})

visualize = matlabFunction(pFR,'File',strcat('visualize',n_name),'Vars',{qR,dimensionparameters})
cd ..
end