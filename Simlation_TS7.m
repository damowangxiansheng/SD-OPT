function [v0,r0,A] = Simlation_TS7(Speed,Steering) 
 dbstop if error
%%
x0 = zeros(6,1);
global u delta
u = Speed;delta = Steering;
%%
tspan = 0:1e-3:100;
%options = odeset('RelTol',1e-3,'AbsTol',1e-4);
[t1,x1]=ode45(@SIM_Vehicle7,tspan,x0);
%%
 v0 = mean(x1(end-20:end,1));
 r0 = mean(x1(end-20:end,2));
 a1 = mean(x1(end-20:end,3));
 a2 = mean(x1(end-20:end,4));
 a3 = mean(x1(end-20:end,5));
 a4 = mean(x1(end-20:end,6));
 A = [a1,a2,a3,a4]
%%
% figure;plot(t1,x1(:,2));  title('Time History');
% xlabel('t/s');ylabel('r ','interpreter','latex','Fontsize',13);
% figure;plot(t1,x1(:,2));  title('Time History');
% xlabel('rad/s');ylabel('r ','interpreter','latex','Fontsize',13);
% hold on 
% plot(tt,Input(:,12)/180*pi,'k');
% figure;plot(tt,Input(:,14)/180*pi,'r');
end

function xdot=SIM_Vehicle7(t,x)
% x 
global delta u
t;
%% 杞﹁締鍙傛暟杈撳叆
Parameters

m = M;
%Nonlinear tyre model:PAC SImpilified
%  Fy = D*sin(C*atan(B*alpha-E*(B*alpha-atan(B*alpha)))); D=mu*Fz; 
%  Cf = B*C*D=CFa*Fz; B = Cfa/(C*mu)
%  Parameter requiring identifications:
%  銆恗u銆戙?併?怌Fa銆戙?併?怌銆戙?併?怑銆?
g = 9.8;
%% Ackerman steering angle
eta = 1; dt = 2*dT;gamma = 0;% toe-in angle
ds = delta;% steering angle original
delta_1 = gamma+ds+eta*dt*ds^2*sign(ds)/l/2;
delta_2 = -gamma+ds-eta*dt*ds^2*sign(ds)/l/2;

%% Load transfer effect
phi = ms*u*h/(Ks1+Ks2-ms*g*h)*x(2);
Fz_fl = Fzf-phi*Ks1/dt; 
Fz_fr = Fzf+phi*Ks1/dt; 
Fz_rl = Fzr-phi*Ks2/dt; 
Fz_rr = Fzr+phi*Ks2/dt; 

%% Transient tyre model (for the delay quantity)
% [Currently not fitted , need to think more]
% [Simple tyre models]
Alpha_fl = delta_1-(x(1)+a*x(2))/(u);
Alpha_fr = delta_2-(x(1)+a*x(2))/(u); 
Alpha_rl =  -(x(1)-b*x(2))/(u);
Alpha_rr =  -(x(1)-b*x(2))/(u);
%% Nonlinear tyre force model
% mu1 = 0.95; cFa1 = 15; C1 = 1.5; E1 = -1;B1 = cFa1/(C1*mu1);
% mu2 = 0.95; cFa2 = 16; C2 = 1.2; E2 = -1;B2 = cFa2/(C2*mu2);%C_f = B*C*D1;% = cFa*Fzf;%B = cFa/(C*mu);%C_r = B*C*D1;% = cFa*Fzr;
%
F_fl = mu1*Fz_fl*sin(C1*atan(B1*x(3)-E1*(B1*x(3)-atan(B1*x(3)))));
F_fr = mu1*Fz_fr*sin(C1*atan(B1*x(4)-E1*(B1*x(4)-atan(B1*x(4)))));
F_rl = mu2*Fz_rl*sin(C2*atan(B2*x(5)-E2*(B2*x(5)-atan(B2*x(5)))));
F_rr = mu2*Fz_rr*sin(C2*atan(B2*x(6)-E2*(B2*x(6)-atan(B2*x(6)))));
%% Dynamic Model of Vehicle 
xdot=[(F_fl*cos(delta_1)+F_fr*cos(delta_2)+(F_rr+F_rl)-x(2)*u*m)/m;
    ((F_fl*cos(delta_1)+F_fr*cos(delta_2))*a-(F_rr+F_rl)*b+(F_fl*sin(delta_1)-F_fr*sin(delta_2))*dt/2)/Iz;
    -u/sig*(x(3)-Alpha_fl);
    -u/sig*(x(4)-Alpha_fr);
    -u/sig*(x(5)-Alpha_rl);
    -u/sig*(x(6)-Alpha_rr);];
end
