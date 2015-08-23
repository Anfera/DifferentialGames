%% Prueba de Differential Game
graphics_toolkit('gnuplot') 
clear all
close all
clc

% Definimos las variables del problema
T = 1;
dt = 0.01;
tspan = 0:dt:T;
x0 = [1;0;2];

a1 = -0;
a2 = -0;
a0 = -0;

b1 = 1;
b2 = 1;

A = [a1 0 0; 0 a2 0; 0 0 a0];
B = [b1 0;0 b2;0 0];
B1 = B(:,1);
B2 = B(:,2);

Q1 = 50*[1 -1 0; -1 1 0; 0 0 0];
Q2 = 5*[0 0 0; 0 1 -1; 0 -1 1];

R1 = 1;
R2 = 1;
S1 = B(:,1)*inv(R1)*B(:,1)';
S2 = B(:,2)*inv(R2)*B(:,2)';

D = [A' zeros(size(A)); zeros(size(A)) A']; 
S = [S1 S2];

Q = [Q1;Q2];

P0 = zeros(size(Q));
%P0 = P0(:);
%[T P] = ode45(@(t,P)mRiccati(t, P, A, S, Q, D), flipud(tspan), P0);
%
%[m n] = size(P);
%PP = mat2cell(P, ones(m,1), n);
%fh_reshape = @(p)reshape(p,size(Q));
%PP = cellfun(fh_reshape,PP,'UniformOutput',false);
%PP = flipud(PP);
%
%for i = 1:length(tspan)
	%aux = PP{i};
	%P1{i} = aux(1:3,:);
	%P2{i} = aux(4:6,:);
	%phi{i} = expm(A-S1*P1{i}-S2*P2{i});
	%u1(i) = -inv(R1)*B1'*P1{i}*phi{i}*[1;0;2];
	%u2(i) = -inv(R2)*B2'*P2{i}*phi{i}*[1;0;2];
%end
%
%x(:,1) = [1;0;2];	
%for i = 1:length(tspan)-1
	%x(:,i+1) = x(:,i) + dt*(A*x(:,i) + B*[u1(i);u2(i)]);
%end

M = [A -S1 -S2; 
        -Q1 -A' zeros(size(A)); 
        -Q2 zeros(size(A)) -A'];
        
P = zeros(3*size(A,1),3*size(A,2));
P(1:3,1:3) = eye(size(A,1));

Q_ = [zeros(size(A,1),3*size(A,2));
          -Q1 eye(3) zeros(3,3);
          -Q2 zeros(3,3) eye(3)];
          
%y_0 = inv(P + Q_*expm(T*M))*[x0;zeros(6,1)];
%
%for i = 1:length(tspan)
	%y(:,i) = expm(M*tspan(i))*y_0;
%end
%
%x = y(1:3,:);

x = zeros(3,5*length(tspan));
x(:,1) = x0;
d = zeros(3,5*length(tspan));
d(2,300) = 0.5;

for i = 1:5*length(tspan)-1
	i
	y_0 = (P + Q_*expm(T*M))\[x(:,i);zeros(6,1)];
	for k = 1:length(tspan)
	y(:,k) = expm(M*tspan(k))*y_0;
	end
	u1(i) = -inv(R1)*B1'*y(4:6,1);
	u2(i) = -inv(R2)*B2'*y(7:9,1);
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