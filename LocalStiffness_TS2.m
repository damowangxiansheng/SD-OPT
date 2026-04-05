function [StiffenssTyre_y,StiffenssTyre_x,Mu_x,SteadyForce] = LocalStiffness_TS2(v0,r0,u,steering) 
%LocalStiffness_TS2(-0.29,0.5,15,0.1) 
Parameters
delta = steering;
%% Load transfer
phi = M*u*h/(Ks1+Ks2)*r0;dt = 2*dT;
Fz_fl = Fzf-phi*Ks1/dt; 
Fz_fr = Fzf+phi*Ks1/dt; 
Fz_rl = Fzr-phi*Ks2/dt; 
Fz_rr = Fzr+phi*Ks2/dt; 
%% Ackerman steering angle
eta = 1; dt = 2*dT;gamma = 0;% toe-in angle
ds = delta;% steering angle original
delta_1 = gamma+ds+eta*dt*ds^2*sign(ds)/l/2;
delta_2 = -gamma+ds-eta*dt*ds^2*sign(ds)/l/2;
%%
a_1 = delta_1-(v0+a*r0)/(u);
a_2 = delta_2-(v0+a*r0)/(u); 
a_3 =  -(v0-b*r0)/(u);
a_4 =  -(v0-b*r0)/(u);
%% Nonlinear tyre force model
% mu1 = 0.95; cFa1 = 15; C1 = 1.5; E1 = -1;B1 = cFa1/(C1*mu1);
% mu2 = 0.95; cFa2 = 16; C2 = 1.2; E2 = -1;B2 = cFa2/(C2*mu2);%C_f = B*C*D1;% = cFa*Fzf;%B = cFa/(C*mu);%C_r = B*C*D1;% = cFa*Fzr;
%
syms Alpha_fl Alpha_fr Alpha_rl Alpha_rr
F_fl = @(Alpha_fl) mu1*Fz_fl*sin(C1*atan(B1*Alpha_fl-E1*(B1*Alpha_fl-atan(B1*Alpha_fl))));
F_fr = @(Alpha_fr) mu1*Fz_fr*sin(C1*atan(B1*Alpha_fr-E1*(B1*Alpha_fr-atan(B1*Alpha_fr))));
F_rl = @(Alpha_rl) mu2*Fz_rl*sin(C2*atan(B2*Alpha_rl-E2*(B2*Alpha_rl-atan(B2*Alpha_rl))));
F_rr = @(Alpha_rr)  mu2*Fz_rr*sin(C2*atan(B2*Alpha_rr-E2*(B2*Alpha_rr-atan(B2*Alpha_rr))));
%
d_F1 = diff(F_fl,Alpha_fl); C1 = double(subs(d_F1,Alpha_fl,a_1));
d_F2 = diff(F_fr,Alpha_fr); C2 = double(subs(d_F2,Alpha_fr,a_2));
d_F3 = diff(F_rl,Alpha_rl); C3 = double(subs(d_F3,Alpha_rl,a_3));
d_F4 = diff(F_rr,Alpha_rr); C4 = double(subs(d_F4,Alpha_rr,a_4));
StiffenssTyre_y = [C1,C2,C3,C4] % Stiffness output
%
F1_0 = F_fl(a_1);F2_0 = F_fr(a_2);F3_0 = F_rl(a_3);F4_0 = F_rr(a_4);
SteadyForce = [F1_0,F2_0,F3_0,F4_0]   % Lateral force output
%% Calculation of Longitudinal
mu_fl = (sqrt((mu1*Fz_fl)^2 - F1_0^2))/mu1/Fz_fl;
mu_fr = (sqrt((mu1*Fz_fr)^2 - F2_0^2))/mu1/Fz_fr;    
mu_rl = (sqrt((mu2*Fz_rl)^2 - F3_0^2))/mu2/Fz_rl;
mu_rr = (sqrt((mu2*Fz_rr)^2 - F4_0^2))/mu2/Fz_rr;

Mu_x = [mu_fl,mu_fr,mu_rl,mu_rr] % Longitudinal stiffness output

Dx = 14;
StiffenssTyre_x = [mu_fl*Dx*Fz_fl,mu_fr*Dx*Fz_fr,mu_rl*Dx*Fz_rl,mu_rr*Dx*Fz_rr] % Longitudinal force output
% syms a
% F1 = @(a) a^3+sin(a);
% d_F1 = diff(F1,a)
% F1(1)
end