%% Juego para el tutorial 
graphics_toolkit('gnuplot') 
clear all
close all
clc

% Definimos las variables del problema
T = 5;
dt = 0.01;
tspan = 0:dt:T;
x0 = [1;0];

a1 = -1;
a2 = -1;

b1 = 1;
b2 = 1;

A = [a1 0; 0 a2];
B = [b1 0;0 b2];
B1 = B(:,1);
B2 = B(:,2);

Q1 = 5*[1 -1; -1 1];
Q2 = 10*[1 -1; -1 1];

R1 = 1;
R2 = 1;

S1 = B(:,1)*inv(R1)*B(:,1)';
S2 = B(:,2)*inv(R2)*B(:,2)';

D = [A' zeros(size(A)); zeros(size(A)) A']; 
S = [S1 S2];

Q = [Q1;Q2];

M = [A -S1 -S2; 
        -Q1 -A' zeros(size(A)); 
        -Q2 zeros(size(A)) -A'];
        
P = zeros(3*size(A,1),3*size(A,2));
P(1:size(A,1),1:size(A,2)) = eye(size(A,1));

Q_ = [zeros(size(A,1),3*size(A,2));
          -Q1 eye(size(A,1)) zeros(size(A));
          -Q2 zeros(size(A)) eye(size(A,1))];

x = zeros(size(A,1),length(tspan));
x(:,1) = x0;
d = zeros(size(A,1),length(tspan));
d(2,200) = 0.5;

for i = 1:length(tspan)-1
	i
	y_0 = (P + Q_*expm(T*M))\[x(:,i);zeros(4,1)];
	for k = 1:length(tspan)
	y(:,k) = expm(M*tspan(k))*y_0;
	end
	u1(i) = -inv(R1)*B1'*y(3:4,1);
	u2(i) = -inv(R2)*B2'*y(5:6,1);
	x(:,i+1) = x(:,i) + dt*(A*x(:,i) + B1*u1(i) + B2*u2(i)) + d(:,i);
end

figure
subplot(2,1,1)	
plot(x')
%legend("x1","x2","leader")
subplot(2,1,2)
plot(u1)
hold on
plot(u2,'r')
%legend('u1','u2')