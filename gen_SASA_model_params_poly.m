function [mass,  COM_x,  COM_y, Ixx, Iyy]  = gen_SASA_model_params_poly(X,Y,t,Ln)
% clear all
% clc
colorblind =[0         0    1.0000;...
    1.0000         0         0;...
    1.0000    1.0000         0;...
    0.6602    0.6602    0.6602;...
         0         0         0;...
    1.0000    0.6445         0;...
    1.0000         0    1.0000;...
         0    0.5000    0.5000;...
         0         0    0.5430;...
         0    0.3906         0;...
         0    1.0000    1.0000;...
    0.5977    0.1953    0.7969];
    
close all
%% Analytical (Polygonal only)

Ln = Ln/sum(Ln); % protection against optmiser violating the length constraint
nodes(1) = 0;
mu = 7850;
for k=1:length(Ln)
    nodes(k+1) = nodes(k)+Ln(k);
end
clear xy;
xy = [X' Y'];
xy = [0 0;xy];
xy = [xy;[(fliplr(X))' (fliplr(-Y))']];
xy = [xy;[0 0]];

P = polyshape(xy(:,1),xy(:,2));

for i = 1:(length(nodes)-1)

pgon = polyshape([nodes(i) nodes(i) nodes(i+1) nodes(i+1)],[2 -2 -2 2]);
lengths = nodes(i+1)-nodes(i);
a = nodes(i);


% ------------

% hold on
polyout = intersect(P,pgon);
% x = polyout.Vertices(polyout.Vertices(:,2)>0,1);
% y = polyout.Vertices(polyout.Vertices(:,2)>0,2);
% area(x,y,'LineStyle','none','FaceColor',(0.5*colorblind(i,:)+0.5*ones(1,3)));%,'LineWidth',2
% area(x,-y,'LineStyle','none','FaceColor',(0.5*colorblind(i,:)+0.5*ones(1,3)));%,'LineWidth',2
% set(gca,'FontSize',20)
% xlabel('Length [m]','FontSize',20)%
% ylabel('Width [m]','FontSize',20)%
% axis([-0.05 1.05 -1.0 0.5])

% ----------
polyout.Vertices(:,1) = polyout.Vertices(:,1)-min(polyout.Vertices(:,1));
xy_part = polyout.Vertices;
props = PolygonMoments(xy_part,[],0);
lengths_list(i) = lengths;
mass(i) =  abs(props.Area*t)*1*mu;
COM_y(i) = props.MAy/props.Area;
COM_x(i) = props.MAx/props.Area;
Ixx(i) =   abs(props.Ixx*t*mu);
Iyy(i) =   abs(props.Iyy*t*mu);
% str = {("Section"+" "+dec2rom(i)),strcat("Mass: ",num2str(mass(end))),  strcat("$COM_x$: ",num2str(COM_x(end))),  strcat("$COM_y$: ",num2str(COM_y(end))), strcat("$I_{xx}$: ",num2str(Ixx(end))), strcat("$I_{yy}$: ",num2str(Iyy(end)))};
% text(a+(COM_x(end)),-0.6,str,'FontSize',14,'HorizontalAlignment','center');%,'FontWeight','bold'

end
end
function ans = dec2rom(z)
d = [ 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
c =  {'M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'};
[];
for ii = 1:numel(d)
    if z >= d(ii)    
        ans = [ans,repmat(c{ii},1,fix(z/d(ii)))];
        z = rem(z,d(ii));
    end
end
end