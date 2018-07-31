% % Author:           Diego Gallo
% % Created:          March 2013
% % Last modified:    January 2014
% % Script that takes in input the .dat centerline file, obtains an
% % analytical description of the centerline based on free-knot splines and
% % gives curvature and torsion.
% % Zhao&Shen 2001 J AM Stat Assoc
% % Sangalli et al. 2009 Appl Statist
% % Stein 1981 Ann Statist


%% CENTERLINE FILE SELECTION
%  -------------------------
function [lscurvature] = CenterlineCharacterization(filename, path, init_knots, order)
%filename=sprintf('%s%s',path,filename);
%filename
%path
cl_imp=importdata(path);
% x=flipdim(cl_imp(:,1),1);
% y=flipdim(cl_imp(:,2),1);
% z=flipdim(cl_imp(:,3),1);
x=cl_imp(:,1);  % if you need to flip the centerline, use x=flipdim(cl_imp.data(:,1),1);
y=cl_imp(:,2);  % if you need to flip the centerline, use y=flipdim(cl_imp.data(:,2),1);
z=cl_imp(:,3);  % if you need to flip the centerline, use z=flipdim(cl_imp.data(:,3),1);
d=cat(1,0,cumsum(sqrt(sum(diff([x y z],[],1).^2,2)))); % curvilinear coordinate

%% FREE knot SPLINE DESCRIPTION
%  -----------------------------
fkn=init_knots;
order=order; % continuous derivatives of order 6-2 at the knot
C=4; % penalization constant for Stein unbiased risk estimate, used to add/remove knots CHECK SENSITIVITY

startknts=linspace(min(d),max(d),fkn);
tau=(aptknt(startknts,order));
lsx=spap2(tau,order,d,x);
lsy=spap2(tau,order,d,y);
lsz=spap2(tau,order,d,z);

ls_ka=knot_addition_removal([lsx lsy lsz],[x y z],d,C);
lsfx_ka=fnval(ls_ka(1),d);
lsfy_ka=fnval(ls_ka(2),d);
lsfz_ka=fnval(ls_ka(3),d);


%% FIRST DERIVATIVE
%  ----------------
%  First derivative of x,y,z wrt curvilinear distance d

dlsx=fnder(ls_ka(1),1);
dlsy=fnder(ls_ka(2),1);
dlsz=fnder(ls_ka(3),1);
dlsfx=fnval(dlsx,d);
dlsfy=fnval(dlsy,d);
dlsfz=fnval(dlsz,d);

% Derivatives are also evaluated with central difference method in order to
% validate the analytical description
%% SECOND DERIVATIVE
%  -----------------
%  Second derivative of x,y,z wrt curvilinear distance d

ddlsx=fnder(dlsx,1);
ddlsy=fnder(dlsy,1);
ddlsz=fnder(dlsz,1);
ddlsfx=fnval(ddlsx,d);
ddlsfy=fnval(ddlsy,d);
ddlsfz=fnval(ddlsz,d);
C1xC2_1=ddlsfz.*dlsfy-ddlsfy.*dlsfz;
C1xC2_2=ddlsfx.*dlsfz-ddlsfz.*dlsfx;
C1xC2_3=ddlsfy.*dlsfx-ddlsfx.*dlsfy;

% Curvature and torsion evaluated from derivatives of the free-knot splines derivatives
lscurvature=sqrt(C1xC2_1.^2+C1xC2_2.^2+C1xC2_3.^2)./(dlsfx.^2+dlsfy.^2+dlsfz.^2).^1.5;
end

%% PLOTS
%  -----
% figure, plot(d,cdtorsion), hold on
% plot(d,lstorsion,'r','LineWidth',2);
% xlabel('Curvilinear distance'),ylabel('Torsion'), legend('Discrete','Analitical');
