function [deflections] = deflection(load,Design_Parametersn, Input,B,C,D,G,J)

% x = [sym('q',[1,Input.n]),zeros(1,Input.n)];
% [angles, params, conds] = solve(transpose(J(Design_Parametersn,x))*[0;load]-G(Design_Parametersn,x));
x = sym('q',[1,Input.n]);
x_lims = [transpose(-inf*ones(1,Input.n)) transpose(0*ones(1,Input.n))];
[angles] = vpasolve(1*(transpose(J(Design_Parametersn,x))*[0;load]-G(Design_Parametersn,x)),x_lims);

if isempty(angles.q1)
    qRn = -0.0*ones(1,Input.n); % Joints of Right array
    dqRn = 0*ones(1,Input.n); %Angular rates
    init = transpose([qRn,dqRn]);
    %[t,xt] = ode113(@(t,x)dynamicsf(t,x,Input.n,Design_Parametersn,B,C,D,G,J,0,load),[0 50],init);
    [t,xt] = ode45(@(t,x)dynamicsf(t,x,Input.n,Design_Parametersn,B,C,D,G,J,0,load),[0 5],init);
    deflections = visualize(xt(end,:),Design_Parametersn);
    deflections = deflections(:,end)/100;
else
    fn = fieldnames(angles);
    q = zeros(1,Input.n);
    for k=1:numel(fn)
        q_temp = (angles.(fn{k}));
        q(k) = q_temp(1);
    end
    q = double(q);
    x = [q,zeros(size(q))];
    deflections = visualize(x,Design_Parametersn);
    deflections = deflections(:,end)/100;
end
end


