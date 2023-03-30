% Santorac 1.36
clc
clear all
u1=67.5032;    %linear velocity [m/s]
a=0.001;
b=0.001;
Qmax=18800; %[kw] motor power
j=0;
r1=3.4995e-3;%[m] Sun Shaft size
r2=17.523e-3; % [m] Fixed roller 1 size
r3=18.175e-3; % [m] Fixed roller 2 size
r4=17.523e-3; % [m] Moveable roller size
gamma=10.96; % speed ratio
w1=u1/r1 %[rad/s] Impeller Speed
wsp=w1/gamma; %spin rotation velocity of motor
pmax=1.36e9;%maximum pressure
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus

%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************

%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1

Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact


nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;

% %coefficients to be found from book Engineering Tribology by Batchelor,2013
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
    fc=2*pi*a*b/3*pmax;
    a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
    b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
    tmpa=a;
    tmpb=b;
end
%**********************************************************
%preprocessing for contact point area geometry

%**********************************************************
fiid=fopen('ellipse-connectivity-oil-Santotrac50.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fiid,'%d%d%d%d');
condata=cell2mat(c);
fclose(fiid);

fiid=fopen('ellipse-nodes.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fiid,'%f%f');
nodecoords=cell2mat(c);
fclose(fiid);

[qq,p]=size(condata);
x=zeros(qq,p); %memory preallocation
y=zeros(qq,p); %memory preallocation

for i=1:qq
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);

    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;
xx=centroidx;
yy=centroidy;
%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas
kk=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
% %*******************************************************************
tepry=0;
tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);
%******************************************************************************
delx =2*b/60;      %heat generation length [m], approximate element length
%in the direction of rolling
km=45;     %We have steel;thermal conductivity of rolling element material [W/(m.K)]
rhom=7860;   %steel; density of rolling element material [kg/m3];
cm= 445; %steel; specific heat constant of rolling element material [j/(kg.K)]
alpha=19.3e-9;  %[Pa^-1], value at 80C, reference paper: The determination of the pressure-viscosity coefficient
%of two traction oils using film thickness measurements', by Leeuwen
k=a/b;
%***********
s=0.886;% specific gravity of Santotrac-50 @ 40C
Toil= 413;% [degree K]supply oil temperature
%******************************************************************************
% I used 'The determination of the pressure-viscosity coefficient
%of two traction oils using film thickness measurements' by Harry Van Leeuen as reference
%for atmospheric viscosity below, Santotrac 50
eta0=26.81e-3;%[Pa.s]
%**************************************************************
kf0=0.13; %[W/(m.K)]thermal conductivity of oil under atm pressure
%see page 35 of book Engineering tribology by Batchelor, 2013
%****************************************************************
%Reference: please see paper 'Lubricant thermal conductivity and heat capacity under high
%pressure' by Larsson, I approximated the graph by linear formulation
kf=0.15+0.111*(p*(1e-9)-0.3); %[W/(m.k)]thermal conductivity of oil under contact pressure
for cr=0.001:0.001:0.02
    Umean=0.5*cr*u1;
    tassume=1e5*ones(qq,1);
    q=tassume*Umean;          % heating amount per unit area (q=Tau*u) [W/m2]

    for iteration=1:10 %convergence loop
        %central thickness of oil film (Harock-Dowson equation) [m], see book,
        %Engineering Tribology, by Batchelor, 2013
        hc=Ry*2.69*(Umean*eta0/(Eprime*Ry)).^(0.67)*(alpha*Eprime)^(0.53)*(fc/(Eprime*Ry^2))^(0.067)*(1-0.61*exp(-0.73*k));
        delTinlet=Umean.^2*eta0/(5*kf0);
        delTsurf=sqrt(delx./(2*pi*km*rhom*cm*Umean))*q;
        Tfilm=zeros(qq,1);
        for ii=1:qq
            if (alpha*p(ii)<=25)
                delTfilm=hc/8*(q(ii)/kf(ii));
            elseif ((alpha*p(ii))>25)
                delTfilm=hc/4*(q(ii)/kf(ii));
            end
            Tfilm(ii)=Toil+delTinlet+delTsurf(ii)+delTfilm;
        end
        %The relationships for finding tl (limiting shear) have been taken from
        %paper 'Geometrical optimization of half toroidal continuously variable
        %transmission using particle swarm optimization', by M. Delkhosh et. al
        tau0=2.7e7+3.5e5*(Tfilm-113-273);%it appears in literature that Tfilm has to be in degC so 273 subtracted
        %here and in next relation, as Tfilm is in Kelvin (K)
        a=0.0765-4e-4*(Tfilm-113-273)-1e-3*Umean;
        tl=(tau0+a.*p);
        %*************||||||||||
        teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
        te=teprx;
        for i=1:qq

            if te(i)>=tl(i)
                te(i)=tl(i);
            end
        end
        q=te*Umean;
        % disp(max(Tfilm));
        % disp(max(te));
    end
    ftraction=sum(te.*area);
    mu(kk)=ftraction/fc;
    kk=kk+1;
    %****************|||||||||||
    % disp(max(tl));
end
disp(size(mu));
%**************************************
figure (5)
cr=0.001:0.001:0.02;
length(cr);
length(mu);
p=plot(cr,mu,'r-*','LineWidth',1.5);
% set(gca,'ytick',0:0.02:0.12);
% axis([0 0.02 0 0.106])
title1='Pure Rolling-Oil Santorac,U1=67.9 m/s'
title2='Pmax=1.36 GPa,Toil=140C'
title3='temperature effects included'
title({title1 title2 title3})
xlabel('Creeping Rate in %age (Cr) ');
ylabel('Traction coefficient \mu');

saveas(p,'graph1.png')
%%  Santorac  2.36
clc
clear all
u1=67.5032;    %linear velocity [m/s]
a=0.001;
b=0.001;
Qmax=18800; %[kw] motor power
j=0;
r1=3.4995e-3;%[m] Sun Shaft size
r2=17.523e-3; % [m] Fixed roller 1 size
r3=18.175e-3; % [m] Fixed roller 2 size
r4=17.523e-3; % [m] Moveable roller size
gamma=10.96; % speed ratio
w1=u1/r1 %[rad/s] Impeller Speed
wsp=w1/gamma; %spin rotation velocity of motor
pmax=2.36e9;%maximum pressure
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus

%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************

%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1

Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact


nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;

% %coefficients to be found from book Engineering Tribology by Batchelor,2013
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
    fc=2*pi*a*b/3*pmax;
    a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
    b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
    tmpa=a;
    tmpb=b;
    % disp(2*a);
    % disp(2*b);
end
%**********************************************************
%preprocessing for contact point area geometry

%**********************************************************
fiid=fopen('ellipse-connectivity-oil-Santotrac50.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fiid,'%d%d%d%d');
condata=cell2mat(c);
fclose(fiid);

fiid=fopen('ellipse-nodes.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fiid,'%f%f');
nodecoords=cell2mat(c);
fclose(fiid);

[qq,p]=size(condata);
x=zeros(qq,p); %memory preallocation
y=zeros(qq,p); %memory preallocation

for i=1:qq
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);

    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;
xx=centroidx;
yy=centroidy;
%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas
kk=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
% %*******************************************************************
tepry=0;
tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);
%******************************************************************************
delx =2*b/60;      %heat generation length [m], approximate element length
%in the direction of rolling
km=45;     %We have steel;thermal conductivity of rolling element material [W/(m.K)]
rhom=7860;   %steel; density of rolling element material [kg/m3];
cm= 445; %steel; specific heat constant of rolling element material [j/(kg.K)]
alpha=19.3e-9;  %[Pa^-1], value at 80C, reference paper: The determination of the pressure-viscosity coefficient
%of two traction oils using film thickness measurements', by Leeuwen
k=a/b;
%***********
s=0.886;% specific gravity of Santotrac-50 @ 40C
Toil= 353;% [degree K]supply oil temperature
%******************************************************************************
% I used 'The determination of the pressure-viscosity coefficient
%of two traction oils using film thickness measurements' by Harry Van Leeuen as reference
%for atmospheric viscosity below, Santotrac 50
eta0=26.81e-3;%[Pa.s]
%**************************************************************
kf0=0.13; %[W/(m.K)]thermal conductivity of oil under atm pressure
%see page 35 of book Engineering tribology by Batchelor, 2013
%****************************************************************
%Reference: please see paper 'Lubricant thermal conductivity and heat capacity under high
%pressure' by Larsson, I approximated the graph by linear formulation
kf=0.15+0.111*(p*(1e-9)-0.3); %[W/(m.k)]thermal conductivity of oil under contact pressure
for cr=0.001:0.001:0.02
    Umean=0.5*cr*u1;
    tassume=1e5*ones(qq,1);
    q=tassume*Umean;          % heating amount per unit area (q=Tau*u) [W/m2]

    for iteration=1:10 %convergence loop
        %central thickness of oil film (Harock-Dowson equation) [m], see book,
        %Engineering Tribology, by Batchelor, 2013
        hc=Ry*2.69*(Umean*eta0/(Eprime*Ry)).^(0.67)*(alpha*Eprime)^(0.53)*(fc/(Eprime*Ry^2))^(0.067)*(1-0.61*exp(-0.73*k));
        delTinlet=Umean.^2*eta0/(5*kf0);
        delTsurf=sqrt(delx./(2*pi*km*rhom*cm*Umean))*q;
        Tfilm=zeros(qq,1);
        for ii=1:qq
            if (alpha*p(ii)<=25)
                delTfilm=hc/8*(q(ii)/kf(ii));
            elseif ((alpha*p(ii))>25)
                delTfilm=hc/4*(q(ii)/kf(ii));
            end
            Tfilm(ii)=Toil+delTinlet+delTsurf(ii)+delTfilm;
        end
        %The relationships for finding tl (limiting shear) have been taken from
        %paper 'Geometrical optimization of half toroidal continuously variable
        %transmission using particle swarm optimization', by M. Delkhosh et. al
        tau0=2.7e7+3.5e5*(Tfilm-113-273);%it appears in literature that Tfilm has to be in degC so 273 subtracted
        %here and in next relation, as Tfilm is in Kelvin (K)
        a=0.0765-4e-4*(Tfilm-113-273)-1e-3*Umean;
        tl=(tau0+a.*p);
        %*************||||||||||
        teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
        te=teprx;
        for i=1:qq

            if te(i)>=tl(i)
                te(i)=tl(i);
            end
        end
        q=te*Umean;
        % disp(max(Tfilm));
        % disp(max(te));
    end
    ftraction=sum(te.*area);
    mu(kk)=ftraction/fc;
    kk=kk+1;
    %****************|||||||||||
    % disp(max(tl));
end
disp(size(mu));
%**************************************
figure (6)
cr=0.001:0.001:0.02;
length(cr);
length(mu);
p=plot(cr,mu,'r-*','LineWidth',1.5);
% set(gca,'ytick',0:0.02:0.12);
% axis([0 0.02 0 0.106])
title1='Pure Rolling-Oil Santorac,U1=67 m/s'
title2='Pmax=2.36 GPa,Toil=80C'
title3='temperature effects included'
title({title1 title2 title3})
xlabel('Creeping Rate in %age (Cr) ');
ylabel('Traction coefficient \mu');

saveas(p,'graph2.png')
%% OIl II 2.36 %*************************************************************
%Given geometry and material data
%*************************************************************
clc
clear all
u1=67.5032;    %linear velocity [m/s]
a=0.001;
b=0.001;
Qmax=18800; %[kw] motor power
j=0;
r1=3.4995e-3;%[m] Sun Shaft size
r2=17.523e-3; % [m] Fixed roller 1 size
r3=18.175e-3; % [m] Fixed roller 2 size
r4=17.523e-3; % [m] Moveable roller size
gamma=10.96; % speed ratio
w1=u1/r1 %[rad/s] Impeller Speed
wsp=w1/gamma; %spin rotation velocity of motor
pmax=2.36e9;%maximum pressure
mumax=0.087;
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus

%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************

%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1

Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact

%Poisson ratio, assumed similar materials; steel
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;

% %coefficients to be found from book tribology by Batchelor
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
    fc=2*pi*a*b/3*pmax;
    a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
    b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
    tmpa=a;
    tmpb=b;
    disp(2*a);
    disp(2*b);
    disp(fc);
end

%**********************************************************
%preprocessing for contact point area geometry

%**********************************************************
fid=fopen('ellipse-connectivity-Oil-A-II.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fid,'%d %d %d %d');
condata=cell2mat(c);
fclose(fid);

fid=fopen('ellipse_nodes-oil-A-II.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fid,'%f%f');
nodecoords=cell2mat(c);
fclose(fid);

[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation

for i=1:q
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);

    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;

xx=centroidx;
yy=centroidy;

%**********************************************************************

%calculate shear stress distribution

%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas

k=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
%find limiting shear
tl=mumax*p;

for cr=0:0.001:0.02

    teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
    tepry=0;

    tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
    tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);

    te=teprx;

    for i=1:q

        if te(i)>=tl(i)
            te(i)=tl(i);
        end
    end

    ftraction=sum(te.*area);

    mu(k)=ftraction/fc;
    k=k+1;
end
cr=0:0.001:0.02;
figure (7)
p=plot(cr,mu,'r--','LineWidth',1.5);
set(gca,'ytick',0:0.02:0.12);
axis([0 0.02 0 0.15]);
title1='Pure Rolling, Oil A-II'
title2='U1=67 m/s,Pmax=3.5 Gpa,T=80'
title({title1 title2})
xlabel('Creeping Rate in %age (Cr) ');
ylabel('Traction coefficient \mu');
ylim([0 0.08])
saveas(p,'graph3.png')
%% %% %*************************************************************
%Given geometry and material data
%*************************************************************
clc
clear all
u1=67.5032;    %linear velocity [m/s]
a=0.001;
b=0.001;
Qmax=18800; %[kw] motor power
j=0;
r1=3.4995e-3;%[m] Sun Shaft size
r2=17.523e-3; % [m] Fixed roller 1 size
r3=18.175e-3; % [m] Fixed roller 2 size
r4=17.523e-3; % [m] Moveable roller size
gamma=10.96; % speed ratio
w1=u1/r1 %[rad/s] Impeller Speed
wsp=w1/gamma; %spin rotation velocity of motor
pmax=3.18e9;%maximum pressure
mumax=0.087;
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus

%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************

%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1

Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact

%Poisson ratio, assumed similar materials; steel
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;

% %coefficients to be found from book tribology by Batchelor
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
    fc=2*pi*a*b/3*pmax;
    a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
    b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
    tmpa=a;
    tmpb=b;
    disp(2*a);
    disp(2*b);
    disp(fc);
end

%**********************************************************
%preprocessing for contact point area geometry

%**********************************************************
fid=fopen('ellipse-connectivity-Oil-A-II.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fid,'%d %d %d %d');
condata=cell2mat(c);
fclose(fid);

fid=fopen('ellipse_nodes-oil-A-II.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fid,'%f%f');
nodecoords=cell2mat(c);
fclose(fid);

[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation

for i=1:q
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);

    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;

xx=centroidx;
yy=centroidy;

%**********************************************************************

%calculate shear stress distribution

%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas

k=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
%find limiting shear
tl=mumax*p;

for cr=0:0.001:0.02

    teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
    tepry=0;

    tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
    tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);

    te=teprx;

    for i=1:q

        if te(i)>=tl(i)
            te(i)=tl(i);
        end
    end

    ftraction=sum(te.*area);

    mu(k)=ftraction/fc;
    k=k+1;
end
cr=0:0.001:0.02;
figure (8)
p=plot(cr,mu,'r--','LineWidth',1.5);
set(gca,'ytick',0:0.02:0.12);
axis([0 0.02 0 0.15]);
title1='Pure Rolling, Oil A-II';
title2='U1=67 m/s,Pmax=3.18 Gpa';
title({title1 title2});
xlabel('Creeping Rate in %age (Cr) ');
ylabel('Traction coefficient \mu');
xlim([0 0.014])
saveas(p,'graph4.png')
%% %% %% %*************************************************************
%Given geometry and material data
%*************************************************************
clc
clear all
u1=67.5032;    %linear velocity [m/s]
a=0.001;
b=0.001;
Qmax=18800; %[kw] motor power
j=0;
r1=3.4995e-3;%[m] Sun Shaft size
r2=17.523e-3; % [m] Fixed roller 1 size
r3=18.175e-3; % [m] Fixed roller 2 size
r4=17.523e-3; % [m] Moveable roller size
gamma=10.96; % speed ratio
w1=u1/r1 %[rad/s] Impeller Speed
wsp=w1/gamma; %spin rotation velocity of motor
pmax=2.8e9;%maximum pressure
mumax=0.087;
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus

%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************

%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1

Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact

%Poisson ratio, assumed similar materials; steel
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;

% %coefficients to be found from book tribology by Batchelor
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
    fc=2*pi*a*b/3*pmax;
    a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
    b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
    tmpa=a;
    tmpb=b;
    disp(2*a);
    disp(2*b);
    disp(fc);
end

%**********************************************************
%preprocessing for contact point area geometry

%**********************************************************
fid=fopen('ellipse-connectivity-Oil-A-II.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fid,'%d %d %d %d');
condata=cell2mat(c);
fclose(fid);

fid=fopen('ellipse_nodes-oil-A-II.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fid,'%f%f');
nodecoords=cell2mat(c);
fclose(fid);

[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation

for i=1:q
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);

    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;

xx=centroidx;
yy=centroidy;

%**********************************************************************

%calculate shear stress distribution

%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas

k=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
%find limiting shear
tl=mumax*p;

for cr=0:0.001:0.02

    teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
    tepry=0;

    tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
    tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);

    te=teprx;

    for i=1:q

        if te(i)>=tl(i)
            te(i)=tl(i);
        end
    end

    ftraction=sum(te.*area);

    mu(k)=ftraction/fc;
    k=k+1;
end
cr=0:0.001:0.02;
figure (9)
p=plot(cr,mu,'r--','LineWidth',1.5);
set(gca,'ytick',0:0.02:0.12);
axis([0 0.02 0 0.15]);
title1='Pure Rolling, Oil A-II'
title2='U1=67.5032 m/s,Pmax=2.3 Gpa,T=140'
title({title1 title2})
xlabel('Creeping Rate in %age (Cr) ');
ylabel('Traction coefficient \mu');
xlim([0 0.014])
saveas(p,'graph5.png')
%%
clc
clear all
u1=20;    %linear velocity [m/s]
a=0.001; % assumed contact elliptical area dimensions,
%later corrected thru iteration. see line 49 of this code
b=0.001;
Qmax=21000; %[kw] power
j=0;
r1=3.4995e-3;%[m] Sun Shaft size
r2=17.523e-3; % [m] Fixed roller 1 size
r3=18.175e-3; % [m] Fixed roller 2 size
r4=17.523e-3; % [m] Moveable roller size
gamma=0;
w1=u1/r1;%[rad/s]
wsp=w1*gamma;%spin rotation velocity
pmax=2.5e9;%maximum pressure
mumax=0.087;
nua=0.3; % poisson ratio, both materials are assumed similar
nub=0.3;
Ea=2.308e11; %[Pa] Young's modulus of contacting parts
Eb=2.308e11;
Eprime=(0.5*((1-nua^2)/Ea+(1-nub^2)/Eb))^(-1); %Equivalent Young's modulus
Gprime=0.5*Eprime/(1+nua); %Equivalent shear modulus

%*********************************************************************
%calculate area of contact based on geometry and material properties of
%the contacting bodies
%*********************************************************************

%Radii of the contacting bodies in x and y directions
Rax=10e-3; %[m] Roller 2
Rbx=inf; % [m] roller 1
Ray=70e-3;%[m] roller 2
Rby=70e-3;% [m] roller 1

Rx=Rax; % Rbx is inf
Ry=Ray*Rby/(Ray+Rby);
Rprime=Rx*Ry/(Rx+Ry);
phi=pi/2; %angle between planes of medium radii of the bodies in contact

%Poisson ratio, assumed similar materials; steel
nua=0.3;
nub=0.3;
Ea=2.308e11; %Young's modulus of contacting parts
Eb=2.308e11;

% %coefficients to be found from book tribology by Batchelor
kdash=1.0339*(Ry/Rx)^0.636;
epsdash=1.0003+(0.5968*Rx)/Ry;
zetadash=1.5277+0.623*log(Ry/Rx);
for ii=1:100
    fc=2*pi*a*b/3*pmax;
    a=(6*kdash^2*epsdash*fc*Rprime/(pi*Eprime))^(1/3);
    b=(6*epsdash*fc*Rprime/(pi*kdash*Eprime))^(1/3);
    tmpa=a;
    tmpb=b;
    disp(2*a);
    disp(2*b);
    disp(fc);
end

%**********************************************************
%preprocessing for contact point area geometry

%**********************************************************
fid=fopen('ellipse-connectivity-Oil-A-II.txt','r','l','UTF-8'); %get connectivity data of areas
c=textscan(fid,'%d %d %d %d');
condata=cell2mat(c);
fclose(fid);

fid=fopen('ellipse_nodes-oil-A-II.txt','r','l','UTF-8'); %get coordinate data of areas
c=textscan(fid,'%f%f');
nodecoords=cell2mat(c);
fclose(fid);

[q,p]=size(condata);
x=zeros(q,p); %memory preallocation
y=zeros(q,p); %memory preallocation

for i=1:q
    for jj=1:p
        x(i,jj)=nodecoords(condata(i,jj),1);
        y(i,jj)=nodecoords(condata(i,jj),2);

    end
end
%find areas and centroids of areas (shear is found out at centroids of areas)
area=0.5*((x(:,1).*y(:,2)+x(:,2).*y(:,3)+x(:,3).*y(:,4)+x(:,4).*y(:,1))-(x(:,2).*y(:,1)+x(:,3).*y(:,2)+x(:,4).*y(:,3)+x(:,1).*y(:,4)));
centroidx=(x(:,1)+x(:,2)+x(:,3)+x(:,4))/4;
centroidy=(y(:,1)+y(:,2)+y(:,3)+y(:,4))/4;

xx=centroidx;
yy=centroidy;

%**********************************************************************

%calculate shear stress distribution

%*********************************************************************
%find shear in x and y directions due to pure rolling at all small areas

k=1;
%find load distribution
p=pmax*sqrt(1-(xx/b).^2-(yy/a).^2);
%find limiting shear
tl=mumax*p;

for cr=0:0.001:0.02

    teprx=32*Gprime/(pi*(4-3*nua))*2*cr/(2-cr)*(xx+b)/(2*b); %pure rolling
    tepry=0;

    tespx=-64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*yy/sqrt(a*b); %spin component
    tespy=64*(2-nua)*Gprime*j/(9*pi*(3-2*nua))*xx/sqrt(a*b);

    te=teprx;

    for i=1:q

        if te(i)>=tl(i)
            te(i)=tl(i);
        end
    end

    ftraction=sum(te.*area);

    mu(k)=ftraction/fc;
    k=k+1;
end
cr=0:0.001:0.02;
figure (10)
p=plot(cr,mu,'r--','LineWidth',1.5);
set(gca,'ytick',0:0.02:0.12);
axis([0 0.02 0 0.15]);
title1='Pure Rolling, Oil A-II'
title2='U1=20 m/s,Pmax=2.5 Gpa,Toil=120C'
title3=' Temperature Effects Included'
title({title1 title2 title3})
xlabel('Creeping Rate in %age (Cr) ');
ylabel('Traction coefficient \mu');
saveas(p,'graph6.png')


