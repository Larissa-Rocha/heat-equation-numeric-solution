clc
clear all

%% Equa��o do calor - Barakat e clark

%Par�metros (constantes)
alfa = 0.2;         %
L = 1;              %comprimento do dom�nio
t_f = 0.5;          %tempo final
t=0;

%Vari�veis
NX = 11;            %n�mero de pontos da malha (11 ou 16)
dx = L / (NX-1);    %passo de espa�o 
r = 2.0;           %par�metro de estabilidade (0.25, 0.50, 1, 2)
dt = r*(dx^2)/alfa; %passo de tempo

%Discretiza��o do dom�nio
i = 1 : NX;
x = (i - 1)*dx;     %dom�nio espacial
NT = t_f / dt;      %n�mero do loop do tempo

%Condi��o inicial
u_0 = 100 * sin(pi * x / L);
p_0 = u_0;
q_0 = u_0;
UUU(1,:) = u_0;     %armazenando a condi��o inicial

%Condi��o de contorno
% u_n(1) = 0;
% u_n(NX) = 0;
p_n(1) = 0; q_n(1) = 0;
p_n(NX) = 0; q_n(NX) = 0;

%Implementa��o do m�todo
for k = 2 : NT      %loop no tempo
    
    for j = 2:NX-1  %loop no espa�o para calcular p
        p_n(j) = (r*p_n(j-1) + (1-r)*p_0(j) + r*p_0(j+1))/(1+r);
    end
    
    for j = NX-1:-1:2  %loop no espa�o para calcular q
        q_n(j) = (r*q_0(j-1) + (1-r)*q_0(j) + r*q_n(j+1))/(1+r);
    end   
   
    u_n = (p_n + q_n)/2; %solu��o
    
    %Armazenando os resultados
    %UUU(k,:) = u_n;
    
    %atualizando p e q
    p_0 = p_n;
    q_0 = q_n;
    t=t+dt;
    
    trem=-alfa*pi^2*0.5;
    u_ana = 100*exp(trem)*sin(pi*x);
    
    %plotando os resultados
    plot(u_n,'-b');
    hold on
    plot(u_ana,'-r');
    hold off
    %pause (dt);
    title('Barakat e Clark: r = 2.0');
    grid on;
    xlabel('N�mero de n�s');
    ylabel('u(x,t)');
    legend('Sol. num�rica','Sol anal�tica');
    
end


