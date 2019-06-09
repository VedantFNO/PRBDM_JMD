function [L] = findKL_eig_def(KL, Input,B,C,D,G,J)

desired_omega = [Input.omega];
Ln = KL(end-Input.n:end-1); %length of each link
Lcenter = KL(end);
L_f = 0;

L_total(1) = Lcenter;
L_total= [L_total Ln];
[mass,  COM_x,  COM_y, Ixx, Iyy]  = gen_SASA_model_params_poly(Input.X,Input.Y,Input.t,L_total);

Mn = mass(2:end); %mass of each link
Jn = Iyy(2:end); %inertia of each link
MR_COM_xn = COM_x(2:end); % COM position for each link
MR_COM_yn = COM_y(2:end);

%cd SimMatrixFun
for i=1:size(Input.x,2)
    
    M1 = D([Mn, Ln,Jn ,KL(i)*ones(1,(Input.n)),MR_COM_xn,MR_COM_yn,Lcenter],transpose(Input.x(:,i)));
    K_diag =sqrt(205e9*Jn./Mn);
    K_diag = KL(i)*K_diag;
    G1 = diag([K_diag]);
    A = pinv(M1)*G1;
    eigA = sort(eig(A));
    
    % Design_Parametersn = [Mn, Ln,Jn ,KL(w)*ones(1,(Input.n)),MR_COM_xn,MR_COM_yn,Lcenter]
    %     x = Input.x(:,w);
    %     eigA = modes(x,Design_Parametersn,KL, Input,B,C,D,G);
    
    idx = min(Input.n,length(desired_omega));
    F = [norm((eigA(1:idx) - desired_omega(1:idx).^2))/norm(desired_omega(1:idx).^2)];
    L_f = L_f+F;
end
L_f = L_f/size(Input.x,2)*100;
L_d=0;
for w = 1:size(Input.load,2)
    for i=1:size(Input.x,2)
        %KL(i)
        K_diag =KL(i)*sqrt(205e9*Jn./Mn);
        Design_Parametersn = [Mn, Ln,Jn ,K_diag,MR_COM_xn,MR_COM_yn,Lcenter];
        target = deflection(Input.load(w),Design_Parametersn, Input,B,C,D,G,J);
        
        error = norm(target-Input.deflection(:,w))/norm(Input.deflection(:,w));
        L_d = L_d +error;
    end
end
L_d = (L_d*100)/(size(Input.load,2)*size(Input.x,2));
L = [L_f ;L_d];
%cd ..
end
