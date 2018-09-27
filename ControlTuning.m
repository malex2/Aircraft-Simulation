clc; clear; close all;

% Constants
m = 2521;
g = 9.81;
W = m*g;
S = 16.72254;
c = 1.3533;
b = 12.8016;
CGmin = 3.4;
CGmax = 3.8;
AR = 9.8;
rho = 0.66;
Ix = 685.9;
Iy = 7501.172;
Iz = 7689;
Ixz = -152.3;
zeta = Ix*Iz - Ixz^2;
ueq = 200*0.5144;
theta_eq = 2*pi/180;

% Nessecary Derivatives
CLalphadot = 0;
Cmalphadot = -3;

Xde    = -1106.1;
Zde    =      -47213;
Mde    =  -1.9965e+5;
Zwdot  = -0.25*rho*S*c*CLalphadot;
Mwdot  = 0.25*rho*S*c^2*Cmalphadot;
XdT    = W*0.1;   % Assume 0.1g of thrust for a unit throttle input
ZdT    = 0;
MdT    = 0;
theta0 = 0;
    
Yda = 1964;
Lda = 60991;
Nda = -36.218;

Ydr = 495.81;
Ldr = 973.79;
Ndr = -6865.8;


%%
ALon = [-0.00958   , 0.0839822 , 0      , -9.81;...
    -0.266257  , -0.858599 , 73.0863,     0;...
    0.000153035, -0.0368076, -0.742218,   0;...
    0          ,          0,         1,   0];

ALon2 = [cos(theta_eq), sin(theta_eq), 0, -ueq*sin(theta_eq);-sin(theta_eq), cos(theta_eq), 0, -ueq*cos(theta_eq);ALon];
ALon2 = [zeros(6,2), ALon2];

BLon = [Xde/m, XdT/m; Zde/(m-Zwdot), ZdT/(m-Zwdot);    
    (Mde + (Mwdot*Zde)/(m-Zwdot))/Iy, (MdT + (Mwdot*ZdT)/(m-Zwdot))/Iy; 0, 0];

BLon2 = [0, 0;0, 0;BLon];

CLon_u = [1, 0, 0, 0];
CLon_w = [0, 1, 0, 0];
CLon = [1, 0, 0, 0;0, 1, 0, 0];
DLon = 0;

CLon2_u = [0, 0, 1, 0, 0, 0];
CLon2_z = [0, 1, 0, 0, 0, 0];
CLon2   = [1, 0, 0, 0;0, 1, 0, 0];
DLon2   = 0;


ALat  =[-0.00695086, -0.25754, -73.9408, 9.81;...
    -0.894067  , -24.262 , 9.62173 , 0   ;...
    0.135007   ,  3.12622, -1.51901, 0   ;...
    0          ,  1      , 0       , 0   ];

BLat = [Yda/m                    ,                   Ydr/m;...
         1/zeta*(Iz*Lda + Ixz*Nda), 1/zeta*(Iz*Ldr+Ixz*Ndr);...
         1/zeta*(Ixz*Lda + Ix*Nda), 1/zeta*(Ixz*Ldr+Ix*Ndr);...
         0                        , 0                     ];
     
CLat = [1, 0, 0, 0];
DLat = 0;

%% Augmented Stability
% Short period desirements
zeta_sp = 0.8;
wn_sp = pi;

LonEig1 = -zeta_sp*wn_sp + wn_sp*sqrt(1-zeta_sp^2)*1i;
LonEig2 = -zeta_sp*wn_sp - wn_sp*sqrt(1-zeta_sp^2)*1i;
% Phugoid period desirements
% S&C longitudinal notes
zeta_p = 0.05;
wn_p = 0.2;
LonEig3 = -zeta_p*wn_sp + wn_p*sqrt(1-zeta_p^2)*1i;
LonEig4 = -zeta_p*wn_sp - wn_p*sqrt(1-zeta_p^2)*1i;

% Roll mode desirements
Tconst_max = 1;
zeta_r = 1;
wn_r  = 15;
LatEig1 = -zeta_r*wn_r + wn_r*sqrt(1-zeta_r^2)*1i;
Tconst_r = 1/(zeta_r*wn_r);

% Spiral mode desirements
Thalf_min = 20;
zw = 0.69/Thalf_min;
zeta_s = 1;
wn_s = 0.01;
Thalf_s = 0.69/(zeta_s*wn_s);
LatEig2 = -zeta_s*wn_s + wn_s*sqrt(1-zeta_s^2)*1i;

% Dutch-roll
zeta_dr = 0.3;
wn_dr = 1.5;
LatEig3 = -zeta_s*wn_dr + wn_dr*sqrt(1-zeta_dr^2)*1i;
LatEig4 = -zeta_s*wn_dr - wn_dr*sqrt(1-zeta_dr^2)*1i;


%% Feedback

pdesLon = [LonEig1, LonEig2, LonEig3, LonEig4];
KLon = place(ALon,BLon,pdesLon);
ACLon = ALon - BLon*KLon;

KLon2 = [zeros(2,2), KLon];
ACLon2 = ALon2 - BLon2*KLon2;

pdesLat = [LatEig1, LatEig2, LatEig3, LatEig4];
KLat = place(ALat,BLat,pdesLat);
ACLat = ALat - BLat*KLat;
%% Response
% LonSys = ss(ALon,BLon(:,2),CLon,0);
% step(LonSys);
% 
% figure
% LonSys = ss(ACLon,BLon(:,2),CLon,0);
% step(LonSys);

% a11 = ALon(1,1);
% b12 = BLon(1,2);
% s = tf('s');
% G = b12/(s-a11);
% Kp = 1;
% Kd = 2;
% Ki = 3;
% K = pid(Kp,Ki,Kd);
% 
% H = feedback(G,K);

% if Ki > 0, then Kp > a11/b12 and Kd > -1/b12

a22 = ALon(2,2);
b21 = BLon(2,1);
s = tf('s');
G = b21/(s-a22);
Kp = 1;
Kd = 1;
Ki = 1;
K = pid(Kp,Ki,Kd);

H = feedback(G*K,1);

% State-Space Tuning
% de
kde = 20;
Rde = 1;
Qde = kde*(CLon_w'*CLon_w);

Kssde = lqr(ACLon,BLon(:,1),Qde,Rde);
Kss2de = [zeros(1,2), Kssde];
pde = -inv(CLon_w*((ACLon-BLon(:,1)*Kssde)\BLon(:,1)));

sys_de = ss(ACLon - BLon(:,1)*Kssde,BLon(:,1)*pde,CLon_w,0);
figure
step(sys_de)

% dT
kdT = 20;
RdT = 1;
QdT = kdT*(CLon_u'*CLon_u);

KssdT = lqr(ACLon,BLon(:,2),QdT,RdT);
pdT = -inv(CLon_u*((ACLon-BLon(:,2)*KssdT)\BLon(:,2)));

sys_dT = ss(ACLon - BLon(:,2)*KssdT,BLon(:,2)*pdT,CLon_u,0);
figure
step(sys_dT)

% try both
kLon = 20;
RLon = eye(2,2);
QLon = kLon*(CLon'*CLon);

KssLon = lqr(ACLon,BLon,QLon,RLon);
pLonde = -inv(CLon_w*((ACLon-BLon*KssLon)\BLon(:,1)));
pLondT = -inv(CLon_u*((ACLon-BLon*KssLon)\BLon(:,2)));

pLon2 = -inv(CLon*((ACLon-BLon*KssLon)\BLon));


sys_de1 = ss(ACLon - BLon*KssLon,BLon(:,1)*pLonde,CLon_w,0);
sys_dT1 = ss(ACLon - BLon*KssLon,BLon(:,2)*pLondT,CLon_u,0);

figure
subplot(2,1,1)
step(sys_de1)

subplot(2,1,2)
step(sys_dT1)

Nmat = [ALon, BLon; CLon, zeros(2,2)]\[zeros;1]