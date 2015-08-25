%% Prueba de Differential Game con 5 agentes y 1 lider
clear all
close all
clc

load('Inflows.mat'); %Cargamos las variables de lluvia del problema
% Definimos las variables del problema
T = 0.05; 	% Esta variable es el horizonte de prediccion en horas
tspan = t;	% Tiempo de simulacion definido en SWMM
x0 = zeros(6,1); % Los tanques empiezan descargados

a1 = -0;
a2 = -0;
a3 = -0;
a4 = -0;
a5 = -0;
a0 = -K(6);  % Aqui es necesario añadir el K del último tanque

A = 3600*diag([a1 a2 a3 a4 a5 a0]); % K fue calculado en segundos, por eso el 3600
B = [-1 0 0 0 0; 
    0 -1 0 0 0; 
    1 1 -1 0 0;
    0 0 0 -1 0; 
    0 0 0 0 -1; 
    0 0 1 1 1]; % Esta matriz contiene la interaccion entre tanques, recordar el grafo.
B1 = B(:,1);
B2 = B(:,2);
B3 = B(:,3);
B4 = B(:,4);
B5 = B(:,5);

% Empezamos el differential Game. Todo es sacado de Dynamic Optimization and LQ Differential Games
Q1 = 50*[1 0 -1 0 0 0; 0 0 0 0 0 0; -1 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
Q2 = 50*[0 0 0 0 0 0; 0 1 -1 0 0 0; 0 -1 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
Q3 = 50*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 -1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 -1 0 0 50];
Q4 = 50*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 -1; 0 0 0 0 0 0; 0 0 0 -1 0 50];
Q5 = 50*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 -1; 0 0 0 0 -1 50];

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

% Ahora empezamos a resolver el sistema      
x = zeros(size(A,1),length(tspan));
x(:,1) = x0;

for i = 1:length(tspan)-1
	[i length(tspan)-1]
	
	y_0 = (P + Q_*expm(T*M))\[x(:,i);zeros(5*size(A,1),1)];

	u1(i) = max(0,-inv(R1)*B1'*y_0((1*size(A,1)+1:2*size(A,1)),1));
	u2(i) = max(0,-inv(R2)*B2'*y_0((2*size(A,1)+1:3*size(A,1)),1));
	u3(i) = max(0,-inv(R3)*B3'*y_0((3*size(A,1)+1:4*size(A,1)),1));
	u4(i) = max(0,-inv(R4)*B4'*y_0((4*size(A,1)+1:5*size(A,1)),1));
	u5(i) = max(0,-inv(R5)*B5'*y_0((5*size(A,1)+1:6*size(A,1)),1));
	x(:,i+1) = x(:,i) + dt*(A*x(:,i) + B1*u1(i) + B2*u2(i) + B3*u3(i) + B4*u4(i) + B5*u5(i)) + dt*[inflows(i,1) inflows(i,2) 0 inflows(i,3) inflows(i,4) 0]';
end

% Plots
figure
subplot(2,1,1)	
plot(x')
legend(['v_1'; 'v_2'; 'v_3'; 'v_4'; 'v_5';'v_6']);
subplot(2,1,2)
U = [u1' u2' u3' u4' u5'];
stairs(U)
legend(['q_1'; 'q_2'; 'q_3'; 'q_4'; 'q_5']);
