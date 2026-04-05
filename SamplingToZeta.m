%% Rule: all-semi-discretisation 
dbstop if error
%% Input
clear;clc;

U = 20;
St = 0.1;
[v0,r0,A] = Simlation_TS7(U,St); 
%  Stiffness local
[StiffenssTyre_y,StiffenssTyre_x,Mu_x,SteadyForce] = LocalStiffness_TS2(v0,r0,U,St)
C_ty = StiffenssTyre_y;
C_tx = StiffenssTyre_x;
% Vehicle
Parameters % m-file
%
C_f1 = C_ty(1); C_f2 = C_ty(2); C_r1 = C_ty(3); C_r2 = C_ty(4);
Cx = max(C_tx);
%% Coeeficient matrix;
syms C_1 C_2 C_3 C_4 l_f l_r I m d_T R I_w sigma real
syms beta r omega  real
syms C_x u delta K_r K_v tau real

A = [0 -u 0 C_1/m C_2/m C_3/m C_4/m;
    0 0 C_x*d_T/I C_1*l_f/I C_2*l_f/I -C_3*l_r/I -C_4*l_r/I;
    0 0 -R^2/I_w*C_x/u 0 0 0 0;
    -1/sigma -l_f/sigma 0 -u/sigma 0 0 0;
    -1/sigma -l_f/sigma 0  0 -u/sigma 0 0;
    -1/sigma l_r/sigma  0  0 0 -u/sigma 0;
    -1/sigma l_r/sigma  0  0 0  0 -u/sigma];
%%
TAU = 0.2;
T_s = 0:0.01:2*TAU;
%T_s = 0;
T_c = TAU-T_s/2;
n_s = length(T_s);

    for jjj = 1:length(T_s)
   jjj
    %% Semi-discretisation
        Parameters
        C_f1 = C_ty(1); C_f2 = C_ty(2); C_r1 = C_ty(3); C_r2 = C_ty(4);
        Cx = max(C_tx);
%% Substitue in
A1 = subs(A,{C_x,m,I,l_f,l_r,C_1,C_2,C_3,C_4,d_T,R,u,I_w,sigma},{Cx,M,Iz,a,b,C_f1,C_f2,C_r1,C_r2,1.0,0.4,U,JW,sig});
A1 = double(A1);
%
tau = T_c(jjj)+T_s(jjj);
tau_s = T_s(jjj);
tor = tau;
%tau_s = 0.2;%tau_c = 0.2;
% if T_c(iii)==0.1||T_s(jjj)==0.1
%  h = 0.01;
% elseif T_s(jjj)==0
%     h = T_c(iii)/10;
% else
%     h = min([T_c(iii),T_s(jjj)])/10;
% end
h = 0.05;
if T_s(jjj)/10<h&&T_s(jjj)~=0
    h = T_s(jjj)/10;
elseif T_s(jjj)==0
    h = 0.01;
else
    h = h;
end

stepx = 100; % number of steps of delta
stepy = 100; % number of steps of b0
kv_min = -0.5; % min. value for b0
kv_max = 0.0; % max. value for b0
kr_min = 0.5; % min. value for delta
kr_max = 1; % max. value for delta
options.disp = 0; % options for command eigs
% SEMI-DISCRETIZATION
p = floor(tau_s/h+0.1); % period resolution by the sampling part (not the constant part)
%h=T/p; % discretization step
r = floor(tau/h+1/2); % delay resolution
r1 = r-p;
%% Matrix definition
Ai = A1;dim = length(A1);
%%
D = eye(dim);
Gi = zeros(dim*(r+1)); % constructing matrix Gi
Gi(2*dim+1:dim*(r+1),(dim+1):r*dim) = eye(dim*r-dim);
Gi(1:dim*2,1:dim) = [zeros(dim);D];
for y = 1:stepy+1 % loop for b0
kr = kr_min + (y-1)*(kr_max-kr_min)/stepy;
%B = [0; b0];
for x = 1:stepx+1 % loop for delta
kv = kv_min + (x-1)*(kv_max-kv_min)/stepx;
Phi = eye(dim*(r+1)); % constructing matrix Phi
B = zeros(dim);
%kv = 0.4;
B(3,1:2) = [kv*0.4*0.4*Iz/JW/U/1, -kr*0.4*0.4*Iz/JW/U/1];
for i = 1:p+1
Pi = expm(Ai*h);
if det(Ai) ~= 0 % if inv(Ai) exist
Ri0 = (Pi-eye(dim))*inv(Ai)*B;
else % if inv(Ai) does not exist
Ri0_int = @(s) expm(Ai*(h-s))*B;
Ri0 = quadv(Ri0_int,0,h);
end
 Gi(1:dim,1:dim) = Pi; % constructing matrix Gi
 p_s = i-1;
 Gi(1:dim,dim*r1+dim*p_s+1:dim*r1+dim*(p_s+1)) = Ri0;
Phi0 = Phi; % Phi by effective matrix multiplication
Phi(1:dim,:) = Pi*Phi0(1:dim,:)+[Ri0]*Phi0(dim*r1+dim*p_s+1:dim*r1+dim*(p_s+1),:);
Phi(dim+1:2*dim,:) = D*Phi0(1:dim,:);
Phi(2*dim+1:dim*(r+1),:) = Phi0(dim+1:dim*r,:);
end
kr_m(x,y) = kr; % matrix of delta
kv_m(x,y) = kv; % matrix of b0
if r*dim+dim > 60 % eigenvalue calculation
try % for large matrices (size > ˜60)
eig_m2(x,y) = abs(eigs(Phi,1,'lm',options)); % by eigs
% PP = eigs(Phi);
% eig_m2(x,y) = PP(PP==max(PP,[],'all'));
catch % in case of convergence failure
PP = eig(Phi);
eig_m(x,y) = max(abs(PP)); % by eig
eig_m2(x,y) = PP(PP==max(PP,[],'all'));
end
else % for small matrices (size < ˜60)
PP = eig(Phi);
eig_m(x,y) = max(abs(PP)); % by eig
eig_m2(x,y) = PP(PP==max(PP,[],'all'));
end
end
stepy+1-y; % counter
end

Re_Eig = real(eig_m2);
Im_Eig = imag(eig_m2);
eig_m3 = sqrt(Re_Eig.^2+Im_Eig.^2);
min_eig = min(eig_m3,[],'all');
[x_loc,y_loc]=find(eig_m3 == min_eig);
E_m = eig_m3(x_loc,y_loc);
if p ==0
    p=1;
end
%SD_zeta = real(log(E_m)/(h*(p)));
SD_zeta = log(min_eig)/(h*(p));

min_kr = kr_m(1,y_loc);
min_kv = kv_m(x_loc,1);

ZETA(jjj) = SD_zeta
KV(jjj) = min_kv;
KR(jjj) = min_kr;

   %% Figure plotting     
% hold on;
%         eig_m2 = 1./eig_m;hold on;
%         contour(kv_m,kr_m,eig_m,[1, 1],'-');
%         box on;
%% For calculating the relative change of zeta w.r.t the continuous case
if jjj>1
 kv_c = KV(1);kr_c = KR(1);

Phi = eye(dim*(r+1)); % constructing matrix Phi
B = zeros(dim);
%kv = 0.4;
B(3,1:2) = [kv_c*0.4*0.4*Iz/JW/U/1, -kr_c*0.4*0.4*Iz/JW/U/1];
for i = 1:p+1
Pi = expm(Ai*h);
if det(Ai) ~= 0 % if inv(Ai) exist
Ri0 = (Pi-eye(dim))*inv(Ai)*B;
else % if inv(Ai) does not exist
Ri0_int = @(s) expm(Ai*(h-s))*B;
Ri0 = quadv(Ri0_int,0,h);
end
 Gi(1:dim,1:dim) = Pi; % constructing matrix Gi
 p_s = i-1;
 Gi(1:dim,dim*r1+dim*p_s+1:dim*r1+dim*(p_s+1)) = Ri0;
Phi0 = Phi; % Phi by effective matrix multiplication
Phi(1:dim,:) = Pi*Phi0(1:dim,:)+[Ri0]*Phi0(dim*r1+dim*p_s+1:dim*r1+dim*(p_s+1),:);
Phi(dim+1:2*dim,:) = D*Phi0(1:dim,:);
Phi(2*dim+1:dim*(r+1),:) = Phi0(dim+1:dim*r,:);
end

if r*dim+dim > 60 % eigenvalue calculation
try % for large matrices (size > ˜60)
eig_mc2 = abs(eigs(Phi,1,'lm',options)); % by eigs
catch % in case of convergence failure
PP = eig(Phi);
eig_mc = max(abs(PP)); % by eig
eig_mc2 = PP(PP==max(PP,[],'all'));
end
else % for small matrices (size < ˜60)
PP = eig(Phi);
eig_mc = max(abs(PP)); % by eig
eig_mc2 = PP(PP==max(PP,[],'all'));
end
Re_Eigc = real(eig_mc2);
Im_Eigc = imag(eig_mc2);
eig_mc3 = sqrt(Re_Eigc.^2+Im_Eigc.^2);
min_eigc = min(eig_mc3,[],'all');
if p ==0
    p=1;
end
%SD_zeta = real(log(E_m)/(h*(p)));
SD_zetac(jjj) = log(min_eigc)/(h*(p))
else
    SD_zetac(jjj) = ZETA(jjj)
end
    end

Zeta2 = smooth(T_s,-ZETA,'moving');
Zeta3 = smooth(T_s,-SD_zetac,'moving');
KR2 = smooth(T_s,KR,'moving');
KV2 = smooth(T_s,KV,'moving');

figure;subplot(3,1,1);
%Per = abs(ZETA-SD_zetac)./abs(SD_zetac);
Per = abs(Zeta2-Zeta3)./abs(Zeta3);
plot(T_s,Per,'k','LineWidth',1);
axis([T_s(1) T_s(end) 0 0.5])
%xlabel('$\tau_{s} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);
 ylabel('Loss %','FontSize',16);


subplot(3,1,2);
plot(T_s,Zeta2,'k','LineWidth',1);hold on;plot(T_s,Zeta3,'k--','LineWidth',1);
%xlabel('$\tau_{s} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);
 ylabel('$\zeta \, \mathrm{[1/s]}$','Interpreter','latex','FontSize',16);
%
subplot(3,1,3);
colororder({'r','k'});
yyaxis left;plot(T_s,KR2,'r','LineWidth',1);hold on;%ylabel('$k_r \, \mathrm{[1/s]}$','Interpreter','latex','FontSize',16);
yyaxis right;plot(T_s,KV2,'k','LineWidth',1);%ylabel('$k_v \, \mathrm{[rad/ms]}$','Interpreter','latex','FontSize',16);
legend({'$k_r \, \mathrm{[1/s]}$','$k_v \, \mathrm{[rad/ms]}$'},'Interpreter','latex','FontSize',14)
xlabel('$\tau_{s} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);

% figure;surf(TT_c,TT_s,-ZETA);
% figure;subplot(1,2,1);surf(TT_c,TT_s,KR);xlabel('$\tau_{c} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);
% ylabel('$\tau_{s} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);zlabel('$k_{Pr} \mathrm{[1/s^2]}$','Interpreter','latex','FontSize',16);
% subplot(1,2,2);surf(TT_c,TT_s,KV);
% xlabel('$\tau_{c} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);
% ylabel('$\tau_{s} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);zlabel('$k_{Pv} \mathrm{[rad/ms^2]}$','Interpreter','latex','FontSize',16)
% 
% 
% figure;surf(TT_c,TT_s,KR);xlabel('$\tau_{c} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);
% ylabel('$\tau_{s} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);zlabel('$k_{Pr} \mathrm{[1/s^2]}$','Interpreter','latex','FontSize',16);
% figure;surf(TT_c,TT_s,KV);
% xlabel('$\tau_{c} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);
% ylabel('$\tau_{s} \, \mathrm{[s]}$','Interpreter','latex','FontSize',16);zlabel('$k_{Pv} \mathrm{[rad/ms^2]}$','Interpreter','latex','FontSize',16)


% Zeta2 = smooth(T_s,-ZETA,'moving');
% Zeta3 = smooth(T_s,-SD_zetac,'moving');
% figure;plot(T_s,Zeta2,'r',T_s,Zeta3,'k')

% Opt_zeta = -3.7149
% Opt_Kr2 =0.9084
% Opt_Kv2 =-0.2759