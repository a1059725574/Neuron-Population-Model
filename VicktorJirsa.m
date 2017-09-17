clc;
clear;
close all;

dt = 0.05;
time = 3000;

% population 1 parameters
a = 1;
b = 3;
c = 1;
d = 5;
s1 = 8;
I1 = 3.1;

% global parameters
N = 40;

x0 = -3; 
r = 0.0002;


noise2 = 0.4;
noise1 = noise2*20;

Ce = 1;
Gs12 = 0.2;
Gs21 = 0.2;
Gs11 = 0.1;
Gs22 = 0.1;


sigma1 = 1/50;
sigma2 = 50;
Tmax = 1;
alphae = 1.1;
alphai = 5;
betae = 0.19;
betai = 0.18;
Vt = 2;
Kp = 5;
Ee = 0;
Ei = -80;

Cm = 20;
I2 = 0.8;
V1 = -1.2;
V2 = 18;
V3 = 12;
V4 = 17.4;
phi = 0.067;

Eca = 120;
Ek = -84;
El = -60;

gca = 4;
gk = 8;
gl = 2;

x1 = 2.5*rand(40, 1) - 1;
y1 = 5*rand(40, 1) - 5;
z  = 0.1*rand(40, 1)+3.45;
x2 = 2*rand(40, 1) - 1;
V  = zeros(40, 1);
n  = zeros(40, 1);
mv = zeros(40, 1);
tauv = zeros(40, 1);
nv = zeros(40, 1);

u1 = 0;
u2 = 0;
u3 = 0;
u4 = 0;

vth = 0;
firing1 = [];
firing2 = [];

x1s = zeros(N, time/dt);
y1s = zeros(N, time/dt);
mz = zeros(1, time/dt);

for t = 2:time/dt
    
    % Isyn 11
    T1 = Tmax./(1+exp( (1/Kp)*(-mean(x1)*50 + Vt) ));
    u1 = ( alphae*T1*(1-u1) - betae*u1 )*dt + u1;
    Isyn11(1:40, 1) = -Gs11*u1*(x1*50 - Ee*ones(40,1));
    
    % Isyn 21
    T2 = Tmax./(1+exp( (1/Kp)*(-mean(x2)*50 + Vt) ));
    u2 = ( alphai*T2*(1-u2) - betai*u2 )*dt + u2;
    Isyn21(1:40, 1) = -Gs21*u2.*(x1*50 - Ei*ones(40, 1));
    
    % Isyn 12
    T3 = Tmax./(1+exp( (1/Kp)*(-mean(x1)*50 + Vt) ));
    u3 = ( alphae*T3*(1-u3) - betae*u3 )*dt + u3;
    Isyn12(1:40, 1) = -Gs12*u3*(x2*50 - Ee*ones(40, 1));
    
    % Isyn 22
    T4 = Tmax./(1+exp( (1/Kp)*(-mean(x2)*50 + Vt) ));
    u4 = ( alphai*T4*(1-u4) - betai*u4 )*dt + u4;
    Isyn22(1:40, 1) = -Gs22*u4*(x2*50 - Ei*ones(40, 1));
    
    x1_next = ( y1 - a*x1.^3 + b*x1.^2 - z + I1*ones(40, 1) + Ce*(mean(x1)*ones(40, 1) - x1) + sigma1*Isyn11 + sigma1*Isyn21...
    + 2*(rand(40, 1) - 0.5*ones(40, 1))*noise1) * dt + x1; % 2*(rand(40, 1) - 0.5*ones(40, 1))*Wmax
    y1_next = (c - d*x1.^2 - y1)*dt + y1;
    z_next  = r*(s1*(x1 + mean(x2)*ones(40, 1) - x0*ones(40, 1)) - mean(z)*ones(40, 1))* dt + z;

    mv = 1/2*(1+ tanh( (V - V1*ones(40,1))./(V2*ones(40, 1)) ));
    tauv = 1./(cosh( (V - V3*ones(40, 1))./(2*V4*ones(40, 1)) ));
    nv = 1/2*(1 + tanh((V - V3*ones(40, 1))./(V4*ones(40, 1))));
        
    V_next  = (I2*50*ones(40, 1) - gl*(V - El*ones(40, 1)) - gk*n.*(V - Ek*ones(40, 1)) -  gca*mv.*(V - Eca) + sigma2*Ce*(mean(x2)*ones(40, 1) - x2)... 
    + Isyn12 + Isyn22 - sigma2*0.3*(mean(z)*ones(40, 1) - 3*ones(40, 1))  )*dt/Cm + 2*sigma2*(rand(40, 1) - 0.5*ones(40, 1))*noise2*dt + V; %    
    n_next  = (phi*(nv - n)./tauv)*dt + n;
    
    x1 = x1_next; 
    y1 = y1_next;
    z = z_next;
    
    V = V_next;
    n = n_next;
    x2 = V./20;
    
    x1s(:, t) = x1;
    y1s(:, t) = y1;
    mz(1, t) = mean(z);
end

theta = angle(x1s+y1s*i);
kop = mean(exp(theta*i), 1);
r_kop1 = abs(kop);
p1_kop = mean(r_kop1);
p1_kop
t1 = dt:dt:time;  
% N1 = 1:40;
% [t1, N1] = meshgrid(t1, N1);

% pcolor(t1, N1, ReV1);
subplot(3, 1, 1), plot(t1, x1s(1, :), 'b')
subplot(3, 1, 2), plot(t1, y1s(1, :), 'b')
subplot(3, 1, 3), plot(t1, mz, 'b')



