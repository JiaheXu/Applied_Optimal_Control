clear;
clc;

% dynamic x_dot = Ax - Bu
A = [ 0 , 1 ; 2 , -1];
B = [ 0 ; 1];

% 1/2 * x'*Q*x
Q = [ 2 , 0 ; 0 , 1 ];
R = 0.005;
Pf = zeros(4,1);
tf = 20;
x0 = [-5 , 5]';
dt = 0.001;

[t,P] = ode45( @(t,P)get_Pdot(t,P,A,B,Q,R) , tf : -dt : 0 , Pf );

len = length(t);
x = zeros(2,len);
u = zeros(len,1);
x(:,len) = [-5,5]';
for i = len:-1:2
    [xdot,u(i)] = get_xdot( t(i) , x(:,i) , A , B , P(i,:) , R );
    x(:,i-1) = x(:,i) + dt*xdot;
end

%plotting P
for i=1:4
    figure;
    plot(t,P(:,i));
    ylabel(strcat('P ',int2str(i),'(t)'));
    xlabel('time s');
end 
% plotting x
for i = 1:2
    figure;
    plot(t,x(i,:));
    ylabel(strcat('x ',int2str(i),'(t)'));
    xlabel('time s');
end
% plotting u
figure;
plot(t,u);
ylabel('u(t)');
xlabel('time s');


function Pdot = get_Pdot(t,P,A,B,Q,R)
    P = reshape(P,2,2);
    Pdot = -(A')*P - P*A + P*B*(B')*P*(1/R) - Q; 
    Pdot = Pdot(:);
end

function [xdot, u] = get_xdot(t,x,A,B,P,R)
    Pt = reshape(P',2,2);
    u = -(1/R)*(B')*Pt*x;
    xdot = A*x + B*u;
end