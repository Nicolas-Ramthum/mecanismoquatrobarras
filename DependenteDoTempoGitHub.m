%ANÁLISE DE MECANISMOS QUATRO BARRAS
%VERSÃO: ANÁLISE DEPENDENTE DO TEMPO
%Matlab v1.0
 
%TRABALHO DE CONCLUSÃO DE CURSO DE ENGENHARIA MECÂNICA
%UNOCHAPECÓ - UNIVERSIDADE COMUNITÁRIA DA REGIÃO DE CHAPECÓ
%TÍTULO: DESENVOLVIMENTO DE PROGRAMA PARA SIMULAÇÃO E ANÁLISE DE MECANISMOS QUATRO BARRAS
%TCC 1
 
%ACADÊMICO: Nícolas Alexandre Gregol Ramthum
%ORIENTADOR: Prof. André Luiz Grando Santos, Dr. Eng.

clear all; % Limpar memoria
clc;
disp('Comeco do programa')

R_1 = 0.4573 %Barra "Terra"
R_2 = 0.15242
R_3 = 0.40444
R_4 = 0.30479
R_5 = 0.40444/2
k = 1; %Incremento do tempo
Teta_2 = 0:k:360
Alpha = 0
Rot = 120
Omega_2 =  Rot*0.104719755120


for i=1:1:length(Teta_2)
S(i) = ((R_1.^2) + (R_2.^2) - 2*R_1*R_2*cosd(Teta_2(i))).^0.5
Beta(i) = acosd(((R_1.^2)+(S(i).^2)-(R_2.^2))/(2*(R_1)*S(i)))
phi(i) = acosd((R_3.^2 - R_4.^2 + S(i).^2)/(2*R_3*S(i)))
lambda(i) = acosd((R_4.^2 - R_3.^2 + S(i).^2)/(2*R_4*S(i)))
gama(i) = acosd((R_3.^2 + R_4.^2 - S(i).^2)/(2*R_3*R_4))
end

for i=1:1:length(Teta_2)
  
if Teta_2(i)<=180 & Teta_2(i)>=0 
Teta_3(i) = phi(i) - Beta(i)
Teta_4(i) = 180 - lambda(i) - Beta(i)
else 
Teta_3(i) = phi(i) + Beta(i)
Teta_4(i) = 180 - lambda(i) + Beta(i) 
end
    
P_x(i) = R_2*cosd(Teta_2(i)) + R_5*cosd(Teta_3(i)+Alpha)
P_y(i) = R_2*sind(Teta_2(i)) + R_5*sind(Teta_3(i)+Alpha)
A_x(i) = R_2*cosd(Teta_2(i)) 
A_y(i) = R_2*sind(Teta_2(i))
B_x(i) = R_2*cosd(Teta_2(i)) + R_3*cosd(Teta_3(i))
B_y(i) = R_2*sind(Teta_2(i)) + R_3*sind(Teta_3(i))

end

N = Omega_2/(2*pi())
DeltaT = (Teta_2(2)- Teta_2(1))/(360*N)
Tempo1V = 2*pi()/Omega_2
NumDiv = Tempo1V/DeltaT
Tempo = 0:DeltaT:Tempo1V

for i=2:1:(length(Teta_2))
   V_Px(i) = (P_x(i)-P_x(i-1))/DeltaT
   V_Py(i) = (P_y(i)-P_y(i-1))/DeltaT
   V_P(i) = (V_Px(i).^2 + V_Py(i).^2).^0.5
   
   V_Ax(i) = (A_x(i)-A_x(i-1))/DeltaT
   V_Ay(i) = (A_y(i)-A_y(i-1))/DeltaT
   V_A(i) = (V_Ax(i).^2 + V_Ay(i).^2).^0.5
   
   V_Bx(i) = (B_x(i)-B_x(i-1))/DeltaT
   V_By(i) = (B_y(i)-B_y(i-1))/DeltaT
   V_B(i) = (V_Bx(i).^2 + V_By(i).^2).^0.5
      
end

V_Px(1) = (P_x(1)-P_x(length(Teta_2)-1))/DeltaT
V_Py(1) = (P_y(1)-P_y(length(Teta_2)-1))/DeltaT
V_P(1) = (V_Px(1).^2 + V_Py(1).^2).^0.5

V_Ax(1) = (A_x(1)-A_x(length(Teta_2)-1))/DeltaT
V_Ay(1) = (A_y(1)-A_y(length(Teta_2)-1))/DeltaT
V_A(1) = (V_Ax(1).^2 + V_Ay(1).^2).^0.5

V_Bx(1) = (B_x(1)-B_x(length(Teta_2)-1))/DeltaT
V_By(1) = (B_y(1)-B_y(length(Teta_2)-1))/DeltaT
V_B(1) = (V_Bx(1).^2 + V_By(1).^2).^0.5


for i=2:1:length(Teta_2)
 A_Px(i) = (V_Px(i)-V_Px(i-1))/DeltaT
 A_Py(i) = (V_Py(i)-V_Py(i-1))/DeltaT
 A_P(i) = (A_Px(i).^2 + A_Py(i).^2).^0.5
 
 A_Ax(i) = (V_Ax(i)-V_Ax(i-1))/DeltaT
 A_Ay(i) = (V_Ay(i)-V_Ay(i-1))/DeltaT
 A_A(i) = (A_Ax(i).^2 + A_Ay(i).^2).^0.5
 
 A_Bx(i) = (V_Bx(i)-V_Bx(i-1))/DeltaT
 A_By(i) = (V_By(i)-V_By(i-1))/DeltaT
 A_B(i) = (A_Bx(i).^2 + A_By(i).^2).^0.5
end

figure(1)

plot(Teta_2,A_x,Teta_2,A_y)
title('Posição do Ponto A em Função de Teta_2')
xlabel('Teta_2 (graus)');
ylabel('Posição (m)');
legend('Posição em X','Posição em Y')

figure(2)

plot(Teta_2,B_x,Teta_2,B_y)
title('Posição do Ponto B em Função de Teta_2')
xlabel('Teta_2 (graus)');
ylabel('Posição (m)');
legend('Posição em X','Posição em Y')

figure(3)

plot(Teta_2,V_A,Teta_2,V_B)
title('Velocidade Linear Resultante no Ponto A e Ponto B em Função de Teta_2')
xlabel('Teta_2 (graus)');
ylabel('Velocidade (m/s)');
legend('Velocidade Ponto A','Velocidade Ponto B')

figure(4)

plot(Teta_2,A_A,Teta_2,A_B)
title('Aceleração Linear Resultante no Ponto A e Ponto B em Função de Teta_2')
xlabel('Teta_2 (graus)');
ylabel('Aceleração (m/s^2)');
legend('Aceleração Ponto A','Aceleração Ponto B')

figure(4)

plot(Teta_2,A_A,Teta_2,A_B)
title('Aceleração Linear Resultante no Ponto A e Ponto B em Função de Teta_2')
xlabel('Teta_2 (graus)');
ylabel('Aceleração (m/s^2)');
legend('Aceleração Ponto A','Aceleração Ponto B')

figure(5)

plot(Teta_2,V_P)
title('Velocidade Linear Resultante no Ponto P em Função de Teta_2')
xlabel('Teta_2 (graus)');
ylabel('Velocidade (m/s)');
legend('Velocidade Ponto P')

figure(6)

plot(Teta_2,A_P)
title('Aceleração Linear Resultante no Ponto P em Função de Teta_2')
xlabel('Teta_2 (graus)');
ylabel('Aceleração (m/s^2)');
legend('Aceleração Ponto P')