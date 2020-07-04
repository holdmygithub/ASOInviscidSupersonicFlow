%program written for Weak Shock
%inputs
close all
clc
clear

tic;
global R
R = 0.0635; %radius of cone
global L
%L = 0.240; %length of cone initial
constant = 4.5;
L = R*constant;


global m1
m1 = 2;
global gamma
gamma = 1.4;
pstatic1 = 97000;
tstatic1 = 300;

global subd
subd = 100; %no. of points, subdivide

cp = 1006.43;
Rconstant = 287.5514; %gasconstant for air 287 J/kg
global h
%h = -0.0005; %spacing of theeta in radians
h = -0.005;
global shockanglespacing
%shockanglespacing = 0.0005;
shockanglespacing = 0.005;

tstagnation = tstatic1*(1+(((gamma-1)/2)*(m1.^2)));
global pstagnation1
pstagnation1 = pstatic1/((1+(((gamma-1)/2)*(m1.^2))).^(gamma/(1-gamma)));
densitystagnation1 = pstagnation1/(Rconstant*tstagnation);
global density1
density1 = densitystagnation1*((1+(((gamma-1)/2)*(m1.^2))).^(1/(1-gamma)));
global v1
v1 = m1*sqrt(gamma*Rconstant*tstatic1);

% h0 = cp*tstagnation;
% vmax = sqrt(2*h0);


%optimization
global lambda
lambda = 0.001;
global de
de = 0.001;

xb = L/3; %initial values of parameters
yb = (R/3);
xc = (2/3)*L;
yc = (2/3)*R;

X = [xb yb xc yc];

tolerance = 0.00000001;
iter = 500; %number of iterations
% hold on;
% plot()
% xlabel('Cycle Number');
% ylabel('X');
% xticks(1:1:iter);
pstatics = zeros(1,subd);
waste = 0
%for optimisation
% run =1;
cdn = 1;
cdnew = 1;
for run = 1:iter

figure(1)    
ax1 = subplot(2,2,1);
ylabel(ax1,'X_{2}')
xlabel(ax1,'Optimization Cycle')
xticks(0:50:iter);
hold(ax1,'on')
plot(ax1,run,X(1),'--*')
drawnow;

figure(1)    
ax2 = subplot(2,2,2);
ylabel(ax2,'Y_{2}')
xlabel(ax2,'Optimization Cycle')
xticks(0:50:iter);
hold(ax2,'on')
plot(ax2,run,X(2),'--*')
drawnow;

figure(1)    
ax3 = subplot(2,2,3);
ylabel(ax3,'X_{3}')
xlabel(ax3,'Optimization Cycle')
xticks(0:50:iter);
hold(ax3,'on')
plot(ax3,run,X(3),'--*')
drawnow;

figure(1)    
ax4 = subplot(2,2,4);
ylabel(ax4,'Y_{3}')
xlabel(ax4,'Optimization Cycle')
xticks(0:50:iter);
hold(ax4,'on')
plot(ax4,run,X(4),'--*')
drawnow;

Xoldcheck = X;
cdn = cdnew;



for k = 1:4    
X(k) = optimize(run,X,k);
disp(X)
end
%cdn for n    


cntrl = [0 X(1) X(3) L;0 X(2) X(4) R];
p = curvy(cntrl,subd);
x = p(1,:);
y = p(2,:);

figure(2)
ax2 = axes;
hold(ax2,'off')
plot(ax2,x,y)
drawnow;

drag=zeros(1,subd);
for i = 2:subd
   [drag(i),pstatics(i)]=elementdrag(run,X,x(i-1),x(i),y(i-1),y(i)); 
end
cdnew = (2*sum(drag))/(density1*(pi*(R.^2))*(v1.^2));

if(abs((cdn-cdnew)/cdn)<=tolerance)
     break;
end

if(cdnew>cdn)
     lambda = 0.5*lambda;
else
     lambda = 1.4*lambda;
end

% plot(x,y)
% drawnow;
end
fprintf('\nThe optimised value of X is: %f %f %f %f\n and cdnew is: %f\n ',X,cdnew);
plot(x,y)
toc;


function pass = optimize(run,X,k)
global lambda
global subd
global density1
global L
global R
global v1
global de


cntrl = [0 X(1) X(3) L;0 X(2) X(4) R];
p = curvy(cntrl,subd);
x = p(1,:);
y = p(2,:);

drag=zeros(1,subd);
for i = 2:subd
   [drag(i),waste]=elementdrag(run,X,x(i-1),x(i),y(i-1),y(i)); 
end
cdn = (2*sum(drag))/(density1*(pi*(R.^2))*(v1.^2));

Xold = X(k);
X(k) = X(k)+de;
cntrl = [0 X(1) X(3) L;0 X(2) X(4) R];
p = curvy(cntrl,subd);
x = p(1,:);
y = p(2,:);


drag=zeros(1,subd);
for i = 2:subd
   [drag(i),waste] =elementdrag(run,X,x(i-1),x(i),y(i-1),y(i)); 
end
cdperturb = (2*sum(drag))/(density1*(pi*(R.^2))*(v1.^2));

fprintf('cdperturb =%f , cdn =%f ',cdperturb,cdn)
% pause(5);
pass = Xold - ((lambda)*((cdperturb-cdn)/de));

end

%for drag calculation of an element
function [drag,pstatics] = elementdrag(run,X,x1,x2,y1,y2)
global m1
global gamma
global h
global shockanglespacing
global pstagnation1
global lambda

coneangle = atan((y2-y1)/(x2-x1));
disp('PLEASE WAIT..... PROGRAM RUNNING...')
for shockangle = asin(1/m1):shockanglespacing:(pi/2)
    devangle = atan(abs((2*cot(shockangle)*(((m1^2)*(sin(shockangle)^2))-1))/(((m1^2)*(gamma+cos(2*shockangle)))+2)));

%calculations for final mach number
    mn1 = m1*sin(shockangle);
    mn2 = sqrt(((mn1^2)+(2/(gamma-1)))/(((mn1^2)*((2*gamma)/(gamma-1)))-1));
    m2 = mn2/sin(shockangle-devangle);

    %calculating vdash and components (velocity normalized to vmax = sqrt(2h0))
    vdash = 1/sqrt((2/((gamma-1)*(m2^2)))+1);
    vrdashinitial = vdash*cos(shockangle-devangle);
    vthetadashinitial = -1*vdash*sin(shockangle-devangle);

    %getting ready for iteration
    theta = shockangle;
    yi = vrdashinitial;
    ydashi = vthetadashinitial;

    while(ydashi<0)

        %calculating ydashiplus1
        f1 = f(gamma,theta,yi,ydashi);
        f2 = f(gamma,(theta+(h/2)),(yi+(ydashi*h/2)),(ydashi+(f1*h/2)));
        f3 = f(gamma,(theta+(h/2)),(yi+(ydashi*h/2)+(f1*h*h/4)),(ydashi+(f2*h/2)));
        f4 = f(gamma,(theta+h),(yi+(ydashi*h)+(f2*h/2)),(ydashi+(f3*h)));

        ydashi = ydashi + ((f1 + (2*f2)+(2*f3)+f4)*h/6);
        yi = yi + (ydashi*h) + ((f1+f2+f3)*h*h/6);
        theta = theta + h;

    end
    if(theta>coneangle)
        break;
    end
end
clc; %clears the pleasewait..

ms = sqrt(2/((gamma-1)*((yi.^(-2))-1))); %mach number at the surface

% disp('theta'); 
% theta*180/pi
% disp("coneangle "); 
% coneangle*180/pi
% disp("shockangle "); 
% shockangle*180/pi

%calculating pressure static
% pstatic2 = pstatic1*(((2*gamma*(power(m1,2))*(power(sin(shockangle),2))) - (gamma-1))/(gamma+1));
pstagnation2 = pstagnation1*((((gamma+1)*(m1.^2)*(sin(shockangle).^2))/((gamma-1)*((m1.^2)*(sin(shockangle).^2))+2)).^(gamma/(gamma-1)))*(((gamma+1)/(2*gamma*((m1*sin(shockangle)).^2)-gamma+1)).^(1/(gamma-1)));
% density2 = density1*(((gamma+1)*(m1.^2)*(sin(shockangle).^2))/(((gamma-1)*((m1.^2)*(sin(shockangle).^2))) + 2  ));
% tstatic2 = tstagnation/(1+(((gamma-1)/2)*(m2.^2)));

%calculations on surface
pstatics = pstagnation2/((1+(((gamma-1)/2)*(ms.^2))).^(gamma/(gamma-1)));
% densitys = density1*(((gamma+1)*(m1.^2)*(sin(shockangle).^2))/(((gamma-1)*((m1.^2)*(sin(shockangle).^2))) + 2  ));
% tstatics = tstagnation/(1+(((gamma-1)/2)*(m2.^2)));

basearea = pi*((y2.^2)-(y1.^2));

drag = (pstatics*sin(coneangle)*basearea);
% cdn = drag/(0.5*density1*basearea*(v1.^2));

%ansys support
% k = sin(coneangle)*basearea/(0.5*density1*basearea*(v1.^2));

fprintf("Run Number = %i\nconeangle = %f deg\nshockangle = %f deg\nM Surface = %f\nPstatic on Surface = %f Pa\nElement Drag = %f\nlambda = %f\nX = %f %f %f %f ",run,(coneangle*180/pi),(shockangle*180/pi),ms,pstatics,drag,lambda,X)
fprintf("\n")
end 

%for taylor maccoll equation solution
function x = f(gamma,theta,yi,ydashi)
    x = ((((gamma-1)/2)*(1-power(yi,2)-power(ydashi,2))*((2*yi)+(ydashi*cot(theta))))-(ydashi*ydashi*yi))/((((1-gamma)/2)*(1-power(yi,2)-power(ydashi,2)))+(power(ydashi,2)));
end


