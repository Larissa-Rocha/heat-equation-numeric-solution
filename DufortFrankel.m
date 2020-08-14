clc
clear all

%% Equação do calor - Dufort-Frankel

%Parâmetros (constantes)
alfa = 0.2;         %
L = 1;              %comprimento do domínio
t_f = 0.5;          %tempo final

%Variáveis
NX = 11;            %número de pontos da malha (11 ou 16)
dx = L / (NX-1);    %passo de espaço 
r = 2.0;           %parâmetro de estabilidade (0.25, 0.50, 1, 2)
dt = r*(dx^2)/alfa; %passo de tempo
t = 0;              %inicialização da variável "tempo" em zero

%Discretização do domínio
i = 1 : NX;
x = (i - 1)*dx;     %domínio espacial
NT = t_f / dt;      %número do loop do tempo

%Condição inicial
u_0 = 100 * sin(pi * x / L);
p_0 = u_0; q_0 = u_0;
UUU(1,:) = u_0;     %armazenando os valores iniciais

   
%Calculando o primeiro passo de tempo usando Barakat e Clark
%Condição de contorno
p_n(1) = 0; q_n(1)=0;
p_n(NX) = 0; q_n(NX)=0;


for j = 2:NX-1  %loop no espaço para calcular p
        p_n(j) = (r*p_n(j-1) + (1-r)*p_0(j) + r*p_0(j+1))/(1+r);
end
    
for j = NX-1:-1:2  %loop no espaço para calcular q
        q_n(j) = (r*q_0(j-1) + (1-r)*q_0(j) + r*q_n(j+1))/(1+r);
end   
   
u_1 = (p_n + q_n)/2; %solução nós internos
UUU(2,:) = u_1;     %armazenando o primeiro passo

%Calculando os passos seguintes
for k = 3 : NT 
    for j = 2 : NX-1
        u_n(j) = (u_0(j) + 2*r*(u_1(j+1) - u_0(j) +u_1(j-1)))/(1+2*r);
    end
    
    %contornos 
    u_n(1)=0; u_n(NX)=0;
    
    %atualizações
    u_0 = u_1;
    u_1 = u_n;
    t=t+dt;
    
    trem=-alfa*pi^2*0.5;
    u_ana = 100*exp(trem)*sin(pi*x);
    
end

%plotando os resultados
plot(u_n,'-b');
    hold on
    plot(u_ana,'-r');
    hold off
    %pause (dt);
    title('Dufort-Frankel: r = 2.0');
    grid on;
    xlabel('Número de nós');
    ylabel('u(x,t)');
    legend('Sol. numérica','Sol analítica');
