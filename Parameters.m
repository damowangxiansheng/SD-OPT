%% The parameter file for all. Please compile this first. Thanks. 
% (Remember to give the driver's inputs)
%% Vehicle Para
M = 1435; Iz = 2340;
a=1.2;b=1.3;g=10;l=2.5;dT = 1;h=0.6;
%% driver's inputs 
% delta = 0.1;
% U = 15;
%% Tyre f
% C1=1.6;E1=-1.3;Road_u1 = 1;%前轮参数
% B1 = 8.5;
%% Tyre r
% C2=1.6;E2=-1.3;Road_u2 =  0.6;%后轮参数
% B2 = 8.5;

mu1 = 0.95; cFa1 = 15; C1 = 1.5; E1 = -1;B1 = cFa1/(C1*mu1);
mu2 = 0.95; cFa2 = 16; C2 = 1.2; E2 = -1;B2 = cFa2/(C2*mu2);%C_f = B*C*D1;% = cFa*Fzf;%B = cFa/(C*mu);%C_r = B*C*D1;% = cFa*Fzr;
%% Original Load
Fzf = M*g*b/l/2; Fzr = M*g*a/l/2;
%% Load transfer effect

ms = 1200;hs = 0.6;Ks1 = 1.6e4;Ks2 = 1.6e4;

JW = 5;
%% Longitudinal tyre characterititics 
% Muxx = sqrt(1-rho^2);
C_pac = 1.4;B_pac = 8; R0 = 0.4; sig = 0.65;
%Cx = Muxx*B_x*C_x*Fzf*1.15
% %% Part 2: Find the stable points(Steady-state function)
% rho = 0.2;
% [A]=TyreSteadyState(Road_u1,Road_u2,rho);
% Alpha_f = A(1);
% Alpha_r = A(2);
% del_alpha = Alpha_r-Alpha_f;
% %% Part 3: Calculate the system matrix
% syms alpha1 alpha2
% F_yf = Road_u1*Fzf*sin(C1*atan(B1*alpha1-E1*(B1*alpha1-atan(B1*alpha1))));
% F_yr = Road_u2*Fzr*sin(C2*atan(B2*alpha2-E2*(B2*alpha2-atan(B2*alpha2))));
% diff_f1 = diff(F_yf,alpha1);
% diff_f2 = diff(F_yr,alpha2);
% yf_prime = double(subs(diff_f1,alpha1,Alpha_f));
% yr_prime = double(subs(diff_f2,alpha2,Alpha_r));
% % differentiate
% C1 = 2*yf_prime;
% C2 = 2*yr_prime;
% %Cs = C1*a-C2*b;Cq = C1*a^2+C2*b^2;
% % k = sqrt(Iz/m);
% % %  System matrix 
% % a11 = -(C1+C2)/(u*m);a12 = -(1+(Cs)/(m*u^2));
% % a21 = -Cs/Iz; a22 = -Cq/(u*Iz);
% % b0 = a11+a22;c0 = a11*a22-a12*a21;
% %Cyx = double(subs(F_yf,alpha1,Alpha_f))/Road_u1/Fzf;
