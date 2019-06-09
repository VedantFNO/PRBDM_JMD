%% initial settings -----------------------------------------------------------------------------
tic
clear all
clc
close all
n = 4;
random = false;
gradient = true;
GA = false;
rectangular = false;
choice = 17572%randi(50000)%3904,17572
% plot variables -----------------------------------------------------------------------------
myxlabel = '$x [m]$'; % x label with latex
myylabel = '$y [m]$'; % y label with latex
myzlabel = '$z [m]$'; % x label with latex
mylegend = ''; % legend with latex


fwidth = 550; % figure width in pixels
fheight = 400; % figure height in pixels
fontlabel = 12; % x,y label font size
fontlegend = 12; % x,y legend font size
fonttick = 12; % x,y rick font size
mycolor1 = [0.8500 0.3250 0.0980]; % custom color 1
mycolor2 = [0.4940 0.1840 0.5560]; % custom color 2
wcolor = [1 1 1]; % white color
bcolor = [0 0 0]; % black color

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
hf = figure; % create a new figure and save handle
%hf.Color = wcolor; % change the figure background color

% -----------------------------------------------------------------------------

%% load/generate matrix -----------------------------------------------------------------------------

f_names = dir (strcat('SimMatrixFun/J',num2str(n),'.m'));
cd SimMatrixFun
if(isempty(f_names))
cd ..
gen_arm_dyn_fun(n)
cd SimMatrixFun
end


B = str2func(strcat('B',num2str(n)));
C = str2func(strcat('C',num2str(n)));
D = str2func(strcat('D',num2str(n)));
G = str2func(strcat('G',num2str(n)));
J = str2func(strcat('J',num2str(n)));
parfiles = {strcat('SimMatrixFun/B',num2str(n),'.m'),strcat('SimMatrixFun/C',num2str(n),'.m'),...
    strcat('SimMatrixFun/D',num2str(n),'.m'),strcat('SimMatrixFun/J',num2str(n),'.m'),strcat('SimMatrixFun/G',num2str(n),'.m'),...
    strcat('SimMatrixFun/visualize',num2str(n),'.m')}
visualize = str2func(strcat('visualize',num2str(n)));
cd ..
% -----------------------------------------------------------------------------
%% generate random panel -----------------------------------------------------------------------------
if(random)
sections = 6
        X = [];
        if sections>2
            X = sort(rand(1,sections-2));
        end
        X = [0,X,1];
        if sections==0
            Y = rand()*ones(1,length(X));
        else
            Y = rand(1,length(X));
        end
        t = 1*rand();
        %         vars = gen_SASA_model_params(X,Y,t,n,samples);
        %vars = double(subs(vars_sym,[Xs Ys ts],[X Y t]));
        
% OR make uniform panel -----------------------------------------------------------------------------
elseif(rectangular)
    
X = [0,1];
% Y = 0.5*rand();
Y = 0.1812;
Y = [Y,Y];
t = rand()/20
t = 0.040532
%  OR read panel from lookup -----------------------------------------------------------------------------
else
data = csvread('lookup.csv',choice,0,[choice 0 choice 27]);
leng = (data(1));
X = data(3:3+leng-1);
Y = data(13:13+leng-1)*0.5;
t = data(23)/20;
[q,w] = sort(Y);
end
% -----------------------------------------------------------------------------

%% define design properties -----------------------------------------------------------------------------
qRn = 0*ones(1,n); % Joints of Right array
dqRn = -0.0*ones(1,n); %Angular rates

% For Center Body
b = 2*Y(1);
factor = sqrt((205e9*(b*t^3/12))/(7850*b*t));
Ln = ones(1,n); %length of each link
LCentern = 1;
L_total = [LCentern,Ln];
L_total = L_total/sum(L_total);
if(~rectangular)
candidate = X(w);
candidate= candidate(X(w)<1 & X(w)>0);
L_total = [0 sort(candidate(1:n)) 1];
L_total = diff(L_total);
desired_omega = [48.342;359.4;1013.5;5711.4;7816];%
end
LCentern = L_total(1);
Ln = L_total(2:end);
figure()
[mass,  COM_x,  COM_y, ~, Iyy]  = gen_SASA_model_params_poly(X,Y,t,L_total);

base_freq = 2;
desired_omega = [(1.875^2);(4.694^2);(7.885^2)]*factor;%;
Mn = mass; %mass of each link
Jn = Iyy; %inertia of each link
Kn = 0*ones(1,n); %inertia of each link
MR_COM_xn = COM_x; % COM position for each link
MR_COM_yn = COM_y;
K = 0;

Design_Parametersn = [Mn,Ln,Jn,Kn,MR_COM_xn,MR_COM_yn,LCentern];
init = transpose([qRn,dqRn]);

Input.omega = desired_omega;
Input.DParams = Design_Parametersn;
Input.x = init;
Input.n = n;
Input.deflection = [1-0.0083; -0.1437];
Input.load = -1e6;
% -----------------------------------------------------------------------------

%% Optimizer inits-----------------------------------------------------------------------------
clc

K_list = [] ;
Loss_list = [];
L_list = [];
Init = [];
for i = 0.1:0.01:(1*pi/(2*n))
    qRn = linspace(0,i,n+1);
    %qRn = qRn(2:end);
    qRn = i*ones(1,n);
    init = transpose([qRn,dqRn]);
%     Input.x = init;
    Init = [Init init];
end
Input.x = Init;
Input.X = X;
Input.Y = Y;
Input.t = t;

design_init = [1e-3*ones(1,size(Input.x,2)), Ln, LCentern];
A = [];
b = [];
A_eq = zeros(length(design_init));
A_eq(1,end - n:end) = 1;
b_eq = zeros(1,length(A_eq));
b_eq(1) = 1;
space = 5;
lb = [1e-4*ones(1,size(Input.x,2)) ,1/(10*space)*ones(1,(length(design_init)-size(Input.x,2)))];
ub = [inf*ones(1,size(Input.x,2)) ,(((10*space)-1)/(10*space))*ones(1,(length(design_init)-size(Input.x,2)))];
% -----------------------------------------------------------------------------

%% Start fmincon

poolobj = parpool('AttachedFiles',parfiles)
pctRunOnAll warning('off')

opts = optimoptions('fmincon','Display','iter','ConstraintTolerance', 1e-6,'StepTolerance', 1e-6,'MaxIterations',10000,...
    'OptimalityTolerance',1e-6,'MaxFunctionEvaluations' , 5000 );%,'Algorithm','sqp'
    problem = createOptimProblem('fmincon','x0',design_init,'objective',@(x) ([0,1]*findKL_eig_def(x, Input,B,C,D,G,J)),...
        'Aeq',[A_eq],'beq',[b_eq],'lb',lb,'ub',ub,'options',opts);
    ms = MultiStart('StartPointsToRun','bounds','UseParallel',true);
    ms.MaxTime = 8*60*60;
    tpoints = CustomStartPointSet(design_init_list);
[K_gard,Loss_gard] = run(ms,problem,{tpoints});
FileName=['results/fmincon01',datestr(now, 'yyyy-mm-dd_HH_MM')];
save(FileName,'K_gard','Loss_gard')


% %% Start Multi -----------------------------------------------------------------------------
% %parpool
% 
% poolobj = parpool('AttachedFiles',parfiles)
% 
% %pctRunOnAll warning('off')
% 
% opts = optimoptions('paretosearch','Display','diagnose','ConstraintTolerance', 1e-8,'MaxTime',24*60*60,'UseParallel',true)%
%                     
% [K_ps,Loss_ps,exitflag,output,residuals] = paretosearch(@(x) findKL_eig_def(x, Input,B,C,D,G,J),length(design_init),A,b,A_eq,b_eq,lb,ub,[],opts);
% 
% % -----------------------------------------------------------------------------
% toc
% FileName=['results/',datestr(now, 'dd-mmm-yyyy')];
% save(FileName,'K_ps','Loss_ps','exitflag','output','residuals')
% 
poolobj.delete