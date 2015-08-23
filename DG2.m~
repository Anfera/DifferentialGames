%% Prueba de Differential Game con 5 agentes y 1 lider
graphics_toolkit('gnuplot') 
clear all
close all
clc

% Definimos las variables del problema
T = 1;
dt = 0.1;
tspan = 0:dt:T;
x0 = [rand(1,5) 1]';
veces = 1;

a1 = -0;
a2 = -0;
a3 = -0;
a4 = -0;
a5 = -0;
a0 = -1

b1 = 1;
b2 = 1;
b3 = 1;
b4 = 1;
b5 = 1;

A = diag([a1 a2 a3 a4 a5 a0]);
B = [-b1 0 0 0 0; b1 -b2 0 0 0; 0 b2 -b3 0 0; 0 0 b3 -b4 0; 0 0 0 b4 -b5; 0 0 0 0 b5];
B1 = B(:,1);
B2 = B(:,2);
B3 = B(:,3);
B4 = B(:,4);
B5 = B(:,5);

Q1 = 20*[1 0 -1 0 0 0; 0 0 0 0 0 0; -1 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
Q2 = 20*[0 0 0 0 0 0; 0 1 -1 0 0 0; 0 -1 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
Q3 = 10*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 -1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 -1 0 0 1];
Q4 = 20*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 -1; 0 0 0 0 0 0; 0 0 0 -1 0 1];
Q5 = 20*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 -1; 0 0 0 0 -1 1];

R1 = 1;
R2 = 1;
R3 = 1;
R4 = 1;
R5 = 1;

S1 = B1*inv(R1)*B1';
S2 = B2*inv(R2)*B2';
S3 = B3*inv(R3)*B3';
S4 = B4*inv(R4)*B4';
S5 = B5*inv(R5)*B5';

D = [A' zeros(size(A)) zeros(size(A)) zeros(size(A)) zeros(size(A)); 
zeros(size(A)) A' zeros(size(A)) zeros(size(A)) zeros(size(A));
zeros(size(A)) zeros(size(A)) A' zeros(size(A)) zeros(size(A));
zeros(size(A)) zeros(size(A)) zeros(size(A)) A' zeros(size(A));
zeros(size(A)) zeros(size(A)) zeros(size(A)) zeros(size(A)) A']; 

S = [S1 S2 S3 S4 S5];

Q = [Q1; Q2; Q3; Q4; Q5];

M = [A -S1 -S2 -S3 -S4 -S5; 
        -Q1 -A' zeros(size(A)) zeros(size(A)) zeros(size(A)) zeros(size(A)); 
        -Q2 zeros(size(A)) -A' zeros(size(A)) zeros(size(A)) zeros(size(A));
        -Q3 zeros(size(A)) zeros(size(A)) -A' zeros(size(A)) zeros(size(A));
        -Q4 zeros(size(A)) zeros(size(A)) zeros(size(A)) -A' zeros(size(A));
        -Q5 zeros(size(A)) zeros(size(A)) zeros(size(A)) zeros(size(A)) -A'];
        
P = zeros(size(A,1)*size(A,1),size(A,1)*size(A,2));
P(1:size(A,1),1:size(A,2)) = eye(size(A,1));

Q_ = [zeros(size(A,1),size(A,1)*size(A,2));
          -Q1 eye(size(A,1)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2));
          -Q2 zeros(size(A,1),size(A,2)) eye(size(A,1)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2));
          -Q3 zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) eye(size(A,1)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2));
          -Q4 zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) eye(size(A,1)) zeros(size(A,1),size(A,2));
          -Q5 zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) zeros(size(A,1),size(A,2)) eye(size(A,1))];

x = zeros(size(A,1),veces*length(tspan));
x(:,1) = x0;
d = zeros(size(A,1),veces*length(tspan));
%d(size(A,1),25) = 0.5;

for i = 1:veces*length(tspan)-1
	i
	y_0 = (P + Q_*expm(T*M))\[x(:,i);zeros(5*size(A,1),1)];
	for k = 1:length(tspan)
	y(:,k) = expm(M*tspan(k))*y_0;
	end
	u1(i) = -inv(R1)*B1'*y((1*size(A,1)+1:2*size(A,1)),1);
	u2(i) = -inv(R2)*B2'*y((2*size(A,1)+1:3*size(A,1)),1);
	u3(i) = -inv(R3)*B3'*y((3*size(A,1)+1:4*size(A,1)),1);
	u4(i) = -inv(R4)*B4'*y((4*size(A,1)+1:5*size(A,1)),1);
	u5(i) = -inv(R5)*B5'*y((5*size(A,1)+1:6*size(A,1)),1);
	x(:,i+1) = x(:,i) + dt*(A*x(:,i) + B1*u1(i) + B2*u2(i) + B3*u3(i) + B4*u4(i) + B5*u5(i)) + d(:,i);
end

figure
subplot(2,1,1)	
plot(x')
try
  legend(["v 1"; "v 2"; "v 3"; "v 4"; "v 5";"v 6"], 
          "location", "northeast",
          "orientation", "vertical");
catch
end_try_catch
subplot(2,1,2)
U = [u1' u2' u3' u4' u5'];
stairs(U)
try
  legend(["q 1"; "q 2"; "q 3"; "q 4"; "q 5"], 
          "location", "northeast",
          "orientation", "vertical");
catch
end_try_catch