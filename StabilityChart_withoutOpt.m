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
TAU = 0.2;
tau_s = 0.2;%tau_c = 0.2;
%tau_c = TAU-tau_s/2;
tau_c = TAU-tau_s;
tau = tau_s+tau_c;
h = 0.01;
stepx = 200; % number of steps of delta
stepy = 200; % number of steps of b0
kv_min = -2; % min. value for b0
kv_max = 4; % max. value for b0
kr_min = -8; % min. value for delta
kr_max = 16; % max. value for delta
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
eig_m(x,y) = abs(eigs(Phi,1,'lm',options)); % by eigs
% PP = eigs(Phi);
% eig_m2(x,y) = PP(PP==max(PP,[],'all'));
catch % in case of convergence failure
PP = eig(Phi);
eig_m(x,y) = max(abs(PP)); % by eig
%eig_m2(x,y) = PP(PP==max(PP,[],'all'));
end
else % for small matrices (size < ˜60)
PP = eig(Phi);
eig_m(x,y) = max(abs(PP)); % by eig
%eig_m2(x,y) = PP(PP==max(PP,[],'all'));
end
end
stepy+1-y % counter
end

   %% Figure plotting     
gcf;hold on;
%contour(kv_m,kr_m,sqrt(real(eig_m2).^2+imag(eig_m2).^2),[1, 1],'g--');
contourf(kv_m,kr_m,eig_m,[1, 1],'g');

% Find the minmum eigenvalue
eig_m2 = sqrt(real(eig_m).^2+imag(eig_m).^2);
min_eig = min(eig_m2,[],'all');
[x_loc,y_loc]=find(eig_m2 == min(eig_m2,[],'all'));
%SD_zeta = log(min_eig^2)/tau/2
SD_zeta = log(min_eig)/((p+1)*h);
min_kr = kr_m(1,y_loc);
min_kv = kv_m(x_loc,1);
hold on;

hold on;scatter(min_kv,min_kr,'*r');
box on;
gcf;
xlabel('$k_{v} \, \, \mathrm{[rad/m \cdot s]}$','Interpreter','latex','FontSize',16);
ylabel('$k_{r} \, \, \mathrm{[1/s]}$','Interpreter','latex','FontSize',16);