%% Note:
% This Script is for anyone using MALS, knowing their angles, you need to
% define it, a Zim-plot will be generated in the end.
% Any questions contact: tingl1@andrew.cmu.edu

%% Scripts start
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
plot (xmod,y1,xmod,y2,xmod,y3,xmod,y4,xmod,y5);

%Find polynormial (linear) coefficient of these linear fit
p1c = polyfit(xmod,y1,1);
p2c = polyfit(xmod,y2,1);
p3c = polyfit(xmod,y3,1);
p4c = polyfit(xmod,y4,1);
p5c = polyfit(xmod,y5,1);

%Intercepts of these linear fits are the Kc/R at theta=0

theta_inter = [p1c(2) p2c(2) p3c(2) p4c(2) p5c(2)];

%concentrations

c = [c1 c2 c3 c4 c5];

% Linear relationship of c relative to theta_inter

plot(c,theta_inter);


% Mw of the particle = 1/intercept

fitconstant1 = polyfit(c,theta_inter,1);


mw = 1/fitconstant1(2);
A2 = fitconstant1(1);

% The rest part of this calculation looks at Rg:

% Angles collected in this experiment

xinc0_ang = [0 69.3 79.7 90.0 100.3 110.7 121.2 132.2];
xinc0_radi = xinc0_ang*pi/180;
xinc0 = (sin(xinc0_radi./2)).^2;

% valid angles

ang0 = [p1c(2) p2c(2) p3c(2) p4c(2) p5c(2)];
ang1 = [y1(1) y2(1) y3(1) y4(1) y5(1)];
ang2 = [y1(2) y2(2) y3(2) y4(2) y5(2)];
ang3 = [y1(3) y2(3) y3(3) y4(3) y5(3)];
ang4 = [y1(4) y2(4) y3(4) y4(4) y5(4)];
ang5 = [y1(5) y2(5) y3(5) y4(5) y5(5)];
ang6 = [y1(6) y2(6) y3(6) y4(6) y5(6)];
ang7 = [y1(7) y2(7) y3(7) y4(7) y5(7)];


% Check linear of centration

plot(c,ang0,c,ang1,c,ang2,c,ang3,c,ang4);


% find fit (p1a, polynomial fit by angle)

p0a = polyfit(c,ang0,1);
p1a = polyfit(c,ang1,1);
p2a = polyfit(c,ang2,1);
p3a = polyfit(c,ang3,1);
p4a = polyfit(c,ang4,1);
p5a = polyfit(c,ang5,1);
p6a = polyfit(c,ang6,1);
p7a = polyfit(c,ang7,1);



% collect intercepts at concentration=0

c_inter = [p0a(2),p1a(2) p2a(2) p3a(2) p4a(2) p5a(2) p6a(2) p7a(2)];

% Linear fit of sin(theta/2)^2 versus concentration0 points

plot(xinc0,c_inter,'d');

fitconstant2 = polyfit(xinc0,c_inter,1);

mwrepeat = 1/fitconstant2(2);
rgsq = fitconstant2(1)/fitconstant2(2)/(16*pi^2*42/(3*(690*10^(-9))^2));
rg=sqrt(rgsq)*10^9;




% Zimm Plot

zimang0 = [p0a(2) p1c(2) p2c(2) p3c(2) p4c(2) p5c(2)];
zimang1 = [p1a(2) y1(1) y2(1) y3(1) y4(1) y5(1)];
zimang2 = [p2a(2) y1(2) y2(2) y3(2) y4(2) y5(2)];
zimang3 = [p3a(2) y1(3) y2(3) y3(3) y4(3) y5(3)];
zimang4 = [p4a(2) y1(4) y2(4) y3(4) y4(4) y5(4)];
zimang5 = [p5a(2) y1(5) y2(5) y3(5) y4(5) y5(5)];
zimang6 = [p6a(2) y1(6) y2(6) y3(6) y4(6) y5(6)];
zimang7 = [p7a(2) y1(7) y2(7) y3(7) y4(7) y5(7)];


cmax = 1/10000000;
c0 = [0 c1 c2 c3 c4 c5];
zimx0 = xinc0(1) + c0./cmax;
zimx1 = xinc0(2) + c0./cmax;
zimx2 = xinc0(3) + c0./cmax;
zimx3 = xinc0(4) + c0./cmax;
zimx4 = xinc0(5) + c0./cmax;
zimx5 = xinc0(6) + c0./cmax;
zimx6 = xinc0(7) + c0./cmax;
zimx7 = xinc0(8) + c0./cmax;


fit0=polyfit(zimx0,zimang0,1);
fit1=polyfit(zimx1,zimang1,1);
fit2=polyfit(zimx2,zimang2,1);
fit3=polyfit(zimx3,zimang3,1);
fit4=polyfit(zimx4,zimang4,1);
fit5=polyfit(zimx5,zimang5,1);
fit6=polyfit(zimx6,zimang6,1);
fit7=polyfit(zimx7,zimang7,1);



line0=polyval(fit0,zimx0);
line1=polyval(fit1,zimx1);
line2=polyval(fit2,zimx2);
line3=polyval(fit3,zimx3);
line4=polyval(fit4,zimx4);
line5=polyval(fit5,zimx5);
line6=polyval(fit6,zimx6);
line7=polyval(fit7,zimx7);






plot(zimx0,zimang0,'*r',zimx0,line0,'k',zimx1,zimang1,'*r',zimx1,line1,'k',zimx2,zimang2,'*r',zimx2,line2,'k',zimx3,zimang3,'*r',zimx3,line3,'k',zimx4,zimang4,'*r',zimx4,line4,'k',zimx5,zimang5,'*r',zimx6,line6,'k',zimx7,zimang7,'*r',zimx7,line7,'k')



