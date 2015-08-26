%% Prueba de Differential Game con 5 agentes y 1 lider
clear all
close all
clc

load('Inflows.mat'); %Cargamos las variables de lluvia del problema
% Definimos las variables del problema
T = 1;               % Esta variable es el horizonte de prediccion en horas
tspan = t;           % Tiempo de simulacion definido en SWMM
x0 = zeros(6,1);     % Los tanques empiezan descargados
vmax = vmax/3600;    % El 3600 es por la integración utilizada
%vmax(6) = 1;

a1 = -0;
a2 = -0;
a3 = -0;
a4 = -0;
a5 = -0;
a0 = -K(6);          % Aqui es necesario añadir el K del último tanque

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
Q3 = 50*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 -1; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 -1 0 0 1];
Q4 = 50*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 -1; 0 0 0 0 0 0; 0 0 0 -1 0 1];
Q5 = 50*[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 -1; 0 0 0 0 -1 1];

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

	u1 = max(0,-inv(R1)*B1'*y_0((1*size(A,1)+1:2*size(A,1)),1));
	u2 = max(0,-inv(R2)*B2'*y_0((2*size(A,1)+1:3*size(A,1)),1));
	u3 = max(0,-inv(R3)*B3'*y_0((3*size(A,1)+1:4*size(A,1)),1));
	u4 = max(0,-inv(R4)*B4'*y_0((4*size(A,1)+1:5*size(A,1)),1));
	u5 = max(0,-inv(R5)*B5'*y_0((5*size(A,1)+1:6*size(A,1)),1));
    
	x(:,i+1) = x(:,i) + dt*(A*x(:,i) + B1*u1 + B2*u2 + B3*u3 + B4*u4 + B5*u5 + [inflows(i,1) inflows(i,2) 0 inflows(i,3) inflows(i,4) 0]')./vmax';
    U(i,:) = [u1 u2 u3 u4 u5].*vmax(1,end-1);
    v(i,:) = x(:,i).*vmax';
end

%% Plots
figure
subplot(2,2,1)	
plot(t',x')
legend(['v_1'; 'v_2'; 'v_3'; 'v_4'; 'v_5';'v_6']);
ylabel('% of occupancy')
subplot(2,2,3)
stairs(t(1:end-1),U)
legend(['q_1'; 'q_2'; 'q_3'; 'q_4'; 'q_5']);
xlabel('time [h]')
ylabel('flow [m^3/s]')

%figure
subplot(2,2,[2,4])
plot(t(1:end-1),v)
legend(['v_1'; 'v_2'; 'v_3'; 'v_4'; 'v_5';'v_6']);
xlabel('time [h]')
ylabel('volume [m^3]')


% figure
% subplot(2,3,1)
% plot(t(1:end-1),v(:,1))
% hold on
% plot(t(1:end-1),ones(size(v,1),1)*vmax(1),'--')
% legend('v_1','v_{1max}');
% ylabel('volume [m^3]')
% 
% subplot(2,3,2)
% plot(t(1:end-1),v(:,2))
% hold on
% plot(t(1:end-1),ones(size(v,1),1)*vmax(2),'--')
% legend('v_2', 'v_{2max}');
% 
% subplot(2,3,3)
% plot(t(1:end-1),v(:,3))
% hold on
% plot(t(1:end-1),ones(size(v,1),1)*vmax(3),'--')
% legend('v_3', 'v_{3max}');
% 
% subplot(2,3,4)
% plot(t(1:end-1),v(:,4))
% hold on
% plot(t(1:end-1),ones(size(v,1),1)*vmax(4),'--')
% legend('v_4', 'v_{4max}');
% ylabel('volume [m^3]')
% xlabel('time [h]')
% 
% subplot(2,3,5)
% plot(t(1:end-1),v(:,5))
% hold on
% plot(t(1:end-1),ones(size(v,1),1)*vmax(5),'--')
% legend('v_5', 'v_{5max}');
% xlabel('time [h]')
% 
% subplot(2,3,6)
% plot(t(1:end-1),v(:,6))
% hold on
% plot(t(1:end-1),ones(size(v,1),1)*vmax(6),'--')
% legend('v_6', 'v_{6max}');
% xlabel('time [h]')
