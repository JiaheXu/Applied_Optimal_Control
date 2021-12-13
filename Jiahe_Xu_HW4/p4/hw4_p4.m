clc;
clear;

% x y theta v phi
x0 = [0,0,0,0,0]';
% simulate dynamics for 5 seconds
[ts, xs] = ode45(@point_ode, [0 5], x0,[]);

% legend for trajectory
%uncomment the following line if you want to plot trajectory 
%legend('real path','desired path')

% legend for u
%uncomment the following line if you want to plot u
legend('u1','u2')
xlabel('time s');

% x y theta v phi
function dx = point_ode(t, x)

theta = x(3);
v = x(4);
phi = x(5);

xd = [ t , 2*t , atan2(2,1),sqrt(5),0]';

%uncomment the following line if you want to plot trajectory 
% plot(x(1),x(2),'*b');
% hold on;
% plot(xd(1),xd(2),'*r');
% hold on;

ud = [0,0]';

thetad = xd(3);
vd = xd(4);
phid = xd(5);
A = [0, 0, -vd* sin(thetad), cos(thetad), 0;
 0, 0, vd* cos(thetad), sin(thetad), 0;
 0, 0, 0, tan(phid), vd/(cos(phid)^2);
 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0];

B = [ 0 0; 0 0;0 0; 1 0; 0 1];

Q = diag([5, 5, 0.01, 0.1, 0.1]);
R = diag([0.5, 0.1]);
[K, P] = lqr(A, B, Q, R);
u = -K*(x-xd) + ud;

%uncomment the following line if you want to plot u
plot(t,u(1),'*b');
hold on;
plot(t,u(2),'*r');
hold on;

dx = [ v*cos(theta); v*sin(theta); v*tan(phi); u(1); u(2);];
end