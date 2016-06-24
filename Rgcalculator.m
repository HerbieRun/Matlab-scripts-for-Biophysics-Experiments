%% Note
% This script will provide the value of Rg in the end based on MALS result
% Any questions, contact: tingl1@andrew.cmu.edu


%%
% Load the raw data
load RD.dat;

% Define the angles (stored in xradi vector)
xang = [14.5 25.9 34.8 42.8 51.5 60.0 69.3 79.7 90.0 100.3 110.7 121.2 132.2 142.5 152.5 163.3];
xradi = xang*pi/180;

% (sin (angle/2))^2, stored in x vector
x = (sin(xradi./2)).^2;

% instrument constant K
kcons = 8.1480*10^(-6)/(1.4879^2)*(1.3301^2);



% After checkpoint, truncate those are not in linear dependence of the
% angles

xmod_ang = [69.3 79.7 90.0 100.3 110.7 121.2 132.2];
xmod_radi = xmod_ang*pi/180;
xmod = (sin(xmod_radi./2)).^2;


% K/Rayleigh ratio

cons = kcons/(9.7801*10^(-6));




% Intensity at concentration 1 (get values from 5 repeats) 
% The most important thing is to understand which row corresponds to what
% concentration and what sample!!!

y11 = RD(11,7:13);
y12 = RD(12,7:13);
y13 = RD(13,7:13);
y14 = RD(14,7:13);
y15 = RD(15,7:13);

y1all = [y11;y12;y13;y14;y15];
y1ave = mean(y1all);

c1 = 9.3*10^(-7);
y1 = c1*cons./y1ave;

% Intensity at concentration 2 (get values from 5 repeats)

y21 = RD(21,7:13);
y22 = RD(22,7:13);
y23 = RD(23,7:13);
y24 = RD(24,7:13);
y25 = RD(25,7:13);

y2all = [y21;y22;y23;y24;y25];
y2ave = mean(y2all);

c2= 1.55*10^(-6);
y2 = c2*cons./y2ave;

% Intensity at concentration 3 (get values from 5 repeats)

y31 = RD(26,7:13);
y32 = RD(27,7:13);
y33 = RD(28,7:13);
y34 = RD(29,7:13);
y35 = RD(30,7:13);

y3all = [y31;y32;y33;y34;y35];
y3ave = mean(y3all);

c3 = 1.86*10^(-6);
y3 = c3*cons./y3ave;


% Intensity at concentration 4 (get values from 5 repeats)

y41 = RD(31,7:13);
y42 = RD(32,7:13);
y43 = RD(33,7:13);
y44 = RD(34,7:13);
y45 = RD(35,7:13);

y4all = [y41;y42;y43;y44;y45];
y4ave = mean(y4all);

c4 = 3.1*10^(-6);
y4 = c4*cons./y4ave;


% Intensity at concentration 5 (get values from 5 repeats)

y51 = RD(36,7:13);
y52 = RD(37,7:13);
y53 = RD(38,7:13);
y54 = RD(39,7:13);
y55 = RD(40,7:13);

y5all = [y51;y52;y53;y54;y55];

y5ave = mean(y5all);

c5 = 3.72*10^(-6);
y5 = c5*cons./y5ave;

% Check Linear Dependence of concentrations 

%plot (x,y1,x,y2,x,y3,x,y4,x,y5)

%Check new linear relationship
plot (xmod,y1,'d',xmod,y2,'d',xmod,y3,'d',xmod,y4,'d',xmod,y5,'d');

%Find polynormial (linear) coefficient of these linear fit
p1c = polyfit(xmod,y1,1);
p2c = polyfit(xmod,y2,1);
p3c = polyfit(xmod,y3,1);
p4c = polyfit(xmod,y4,1);
p5c = polyfit(xmod,y5,1);

% Rgs

A2=4.30595*5*10^(-1);


Rg1 = sqrt(p1c(1)/(1/(5.217*10^7)+2*c1*A2)*3*(690*10^(-9))^2/(16*pi^2))*10^9
Rg2 = sqrt(p2c(1)/(1/(5.217*10^7)+2*c2*A2)*3*(690*10^(-9))^2/(16*pi^2))*10^9
Rg3 = sqrt(p3c(1)/(1/(5.217*10^7)+2*c3*A2)*3*(690*10^(-9))^2/(16*pi^2))*10^9
Rg4 = sqrt(p4c(1)/(1/(5.217*10^7)+2*c4*A2)*3*(690*10^(-9))^2/(16*pi^2))*10^9
Rg5 = sqrt(p5c(1)/(1/(5.217*10^7)+2*c5*A2)*3*(690*10^(-9))^2/(16*pi^2))*10^9



