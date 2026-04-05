%% Rule: all-semi- discretisation 
dbstop if error
%% Input
clear;clc;
%figure;
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
syms C_f C_r l_f l_r I m d_T R I_w real
syms r omega real
syms C_x u delta K_r K_v tau real

A=[-(C_f+C_r)/m/u -u-(C_f*l_f-C_r*l_r)/m/u 0;
    -(C_f*l_f-C_r*l_r)/I/u -(C_f*l_f^2+C_r*l_r^2)/u/I C_x*d_T/I;
    0 0 -C_x*R^2/u/I_w];
    %% Semi-discretisation
%% Substitue in
A1 = subs(A,{C_x,m,I,l_f,l_r,C_f,C_r,d_T,R,u,I_w},{Cx,M,Iz,a,b,C_f1+C_f2,C_r1+C_r2,1.0,0.4,U,JW});
A1 = double(A1);
%
TAU = 0.15;
tau_s = 0.3;%tau_c = 0.2;
tau_c = TAU-tau_s/2;
tau = tau_s+tau_c;
h = 0.005;
stepx = 200; % number of steps of delta
stepy = 200; % number of steps of b0
kv_min = -2; % min. value for b0
kv_max = 2; % max. value for b0
kr_min = -4; % min. value for delta
kr_max = 4; % max. value for delta
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
stepy+1-y % counter
end

   %% Figure plotting     
% figure;
% contour(kv_m,kr_m,eig_m,[1, 1],'k')% only the boundary of the eigenvalue equals to 1 
% figure;contourf(kv_m,kr_m,log(eig_m)./tau,[linspace(-10,0,100)],'EdgeColor','none');
% colormap(gray);colorbar;
% hold on;gcf;hold on;
%figure;contourf(kv_m,kr_m,sqrt(real(eig_m2).^2+imag(eig_m2).^2),[linspace(1,0,100)],'EdgeColor','none');
%colormap(gray);colorbar;hold on;
gcf;hold on;contour(kv_m,kr_m,sqrt(real(eig_m2).^2+imag(eig_m2).^2),[1, 1],'y--');

% Find the minmum eigenvalue
min_eig = min(eig_m2,[],'all');
[x_loc,y_loc]=find(eig_m2 == min(eig_m2,[],'all'));
%SD_zeta = log(min_eig^2)/tau/2
SD_zeta = log(min_eig)/((p+1)*h)
min_kr = kr_m(1,y_loc)
min_kv = kv_m(x_loc,1)
hold on;

hold on;scatter(min_kv,min_kr,'*r');
box on;
gcf;
xlabel('$k_{v} \, \, \mathrm{[rad/m \cdot s]}$','Interpreter','latex','FontSize',16);
ylabel('$k_{r} \, \, \mathrm{[1/s]}$','Interpreter','latex','FontSize',16);
hold on;
%% Continuous
%%  Analytical derivation
syms C_f C_r l_f l_r I m d_T R I_w real
syms beta r omega  real
syms C_x u delta K_r K_v tau real
A=[-(C_f+C_r)/m/u -u-(C_f*l_f-C_r*l_r)/m/u 0;
    -(C_f*l_f-C_r*l_r)/I/u -(C_f*l_f^2+C_r*l_r^2)/u/I C_x*d_T/I;
    0 0 -C_x*R^2/u/I_w];
B = [0 0 0;
    0 0 0;
    K_v*R*R*I/I_w/u/d_T -K_r*R*R*I/I_w/u/d_T 0];

syms lambda zeta tau real
Eqn2 = lambda*eye(3)-A-B*exp(-lambda*tau);
CF = det(Eqn2);
CF2 = collect(CF,lambda);
syms omega zeta real
CF3 = subs(CF2,lambda,zeta+1i*omega);
CF4 = simplify(expand(rewrite(CF3,'cos')));
Re_1 = real(CF4);Im_1 = imag(CF4);
SS1_re = subs(Re_1,{omega},{0});
SS1_im = subs(Im_1,{omega},{0});
Eqn1 = SS1_re==0;
SS = solve(Eqn1,K_r);

Eqns = [Re_1==0,Im_1==0];
DS = solve(Eqns,[K_v,K_r]);
DS_V = DS.K_v;DS_R = DS.K_r;
%
Parameters
C_f1 = C_ty(1); C_f2 = C_ty(2); C_r1 = C_ty(3); C_r2 = C_ty(4);
Cx = max(C_tx);

DS_v1 = subs(DS_V,{m,I,l_f,l_r,C_f,C_r,d_T,R,u,I_w},{M,Iz,a,b,C_f1+C_f2,C_r1+C_r2,1.0,0.4,U,JW});
DS_r1 = subs(DS_R,{m,I,l_f,l_r,C_f,C_r,d_T,R,u,I_w},{M,Iz,a,b,C_f1+C_f2,C_r1+C_r2,1.0,0.4,U,JW});
%
tor = TAU;
DS_v2 = subs(DS_v1,{C_x,tau},{Cx,tor});
DS_r2 = subs(DS_r1,{C_x,tau},{Cx,tor});
%- Static boudnary
ww = linspace(1e-3,40,1e4);Ze = 0;
    DS_v3 = subs(DS_v2,{zeta},{Ze});
    DS_r3 = subs(DS_r2,{zeta},{Ze});
 KV = double(subs(DS_v3,omega,ww));
 KR = double(subs(DS_r3,omega,ww));
h = gcf; hold on;
plot(KV,KR,'r');
hold on;
SS1 = subs(SS,{m,I,l_f,l_r,C_f,C_r,d_T,R,u,I_w},{M,Iz,a,b,C_f1+C_f2,C_r1+C_r2,1.0,0.4,U,JW});
SS2 = subs(SS1,{C_x,tau},{Cx,tor});
SS3 = subs(SS2,{zeta},{Ze});
K_vv = linspace(-100,100,100);
KRR_S  = subs(SS3,K_v,K_vv);
hold on;plot(K_vv,KRR_S,'r','Linewidth',1);
axis([-5 5 -12 12])
%%  Optimal point and maximal decay
syms f1(omega) f2(omega)
f1(omega) = DS_v2; f2(omega) = DS_r2;
Df1 = diff(f1,omega); Df2 = diff(f2,omega);

ddf = Df2/Df1;
Lim = limit(ddf,omega,0);
Lim2 = simplify(Lim);

% SS1 = subs(SS,{m,I,l_f,l_r,C_1,C_2,C_3,C_4,d_T,R,u,I_w,sigma},{M,Iz,a,b,C_f1,C_f2,C_r1,C_r2,1.0,0.4,U,JW,sig});
% SS2 = subs(SS1,{tau},{tor});
syms fs(K_v)
fs(K_v) = SS2;
Dsf = diff(fs,K_v);

Eqn3 = Dsf-Lim2==0
SS_zeta = double(solve(Eqn3,zeta))
Opt_zeta = double(SS_zeta(SS_zeta<0&abs(SS_zeta)<=min(abs(SS_zeta)))) 

Opt_Kr = subs(DS_r2,{zeta,omega},{Opt_zeta,1e-8});
Opt_Kv = subs(DS_v2,{zeta,omega},{Opt_zeta,1e-8});

Opt_Kr2 = double(Opt_Kr)
Opt_Kv2 = double(Opt_Kv)
SS_omega = 1e-8;

if Opt_zeta < -23
    Eqn4 = [Df1==0, Df2 ==0];
    SS_2 = vpasolve(Eqn4,[omega,zeta],[5,4]);
    SS_omega = double(SS_2.omega)
    SS_zeta2 = double(SS_2.zeta)
    Opt_zeta = SS_zeta2;
end


Opt_Kr = subs(DS_r2,{zeta,omega},{Opt_zeta,SS_omega});
Opt_Kv = subs(DS_v2,{zeta,omega},{Opt_zeta,SS_omega});

Opt_Kr2 = double(Opt_Kr)
Opt_Kv2 = double(Opt_Kv)
%scatter(1,1,'r','d','filled')
%% Figure plotting 1: zeta = Opt_zeta+2;
%     Set_Zeta = Opt_zeta+0.8;
%     DS_v21 = subs(DS_v2,{zeta},{Set_Zeta});
%     DS_r21 = subs(DS_r2,{zeta},{Set_Zeta});
%     
%  KV2 = double(subs(DS_v21,omega,ww));
%  KR2 = double(subs(DS_r21,omega,ww));
% h = gcf; hold on;
% plot(KV2,KR2,'y');
% hold on;
% 
% SS21 = subs(SS2,{zeta},{Set_Zeta});
% K_vv2 = linspace(-100,100,100);
% KRR_S2  = subs(SS21,K_v,K_vv2);
% hold on;plot(K_vv2,KRR_S2,'y','Linewidth',1);

%% Figure plotting 1: zeta = Opt_zeta+2;
%     Set_Zeta = Opt_zeta+0.4;
%     DS_v21 = subs(DS_v2,{zeta},{Set_Zeta});
%     DS_r21 = subs(DS_r2,{zeta},{Set_Zeta});
%     
%  KV2 = double(subs(DS_v21,omega,ww));
%  KR2 = double(subs(DS_r21,omega,ww));
% h = gcf; hold on;
% plot(KV2,KR2,'g');
% hold on;
% 
% SS21 = subs(SS2,{zeta},{Set_Zeta});
% K_vv2 = linspace(-100,100,100);
% KRR_S2  = subs(SS21,K_v,K_vv2);
% hold on;plot(K_vv2,KRR_S2,'g','Linewidth',1);

%% Figure plotting 1: zeta = Opt_zeta+0.5;
%     Set_Zeta = Opt_zeta+0.1;
%    %Set_Zeta = -3;
%     DS_v21 = subs(DS_v2,{zeta},{Set_Zeta});
%     DS_r21 = subs(DS_r2,{zeta},{Set_Zeta});
%     
%  KV2 = double(subs(DS_v21,omega,ww));
%  KR2 = double(subs(DS_r21,omega,ww));
% h = gcf; hold on;
% plot(KV2,KR2,'b');
% hold on;
% 
% SS21 = subs(SS2,{zeta},{Set_Zeta});
% K_vv2 = linspace(-100,100,100);
% KRR_S2  = subs(SS21,K_v,K_vv2);
% hold on;plot(K_vv2,KRR_S2,'b','Linewidth',1);
%%
hold on; scatter(Opt_Kv2,Opt_Kr2,'rd','filled');
hold on;scatter(min_kv,min_kr,'*r');
box on;
gcf;
xlabel('$k_{v} \, \, \mathrm{[rad/(m \, s)]}$','Interpreter','latex','FontSize',16);
ylabel('$k_{r} \, \, \mathrm{[1/s]}$','Interpreter','latex','FontSize',16);