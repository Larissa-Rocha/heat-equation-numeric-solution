clc
clear all

%% Equa��o do calor - Euler expl�cito

%Par�metros (constantes)
alfa = 0.2;         %condutividade
L = 1;              %comprimento do dom�nio
t_f = 0.5;          %tempo final

%Vari�veis
NX = 11;            %n�mero de pontos da malha (11 ou 16)
dx = L / (NX-1);    %passo de espa�o 
t = 0;              %inicializa��o da vari�vel "tempo"
r = 2.0;           %par�metro de estabilidade (0.25, 0.50, 1, 2)
dt = r*(dx^2)/alfa; %passo de tempo

%Discretiza��o do dom�nio
i = 1 : NX;
x = (i - 1)*dx;     %dom�nio espacial

%Condi��o inicial
u_0 = 100 * sin(pi * x / L);
u_num = u_0;        %inicializa��o de u_num na condi��o inicial

NT = t_f / dt;      %n�mero do loop do tempo
%Implementa��o do m�todo
for k = 2 : NT %loop no tempo
    
    for j = 2 : NX-1
        u_num(j) = u_0(j) + r*(u_0(j+1) - 2*u_0(j) + u_0(j-1)); 
    end
    
    %Condi��o de contorno
    u_num(1) = 0;
    u_num(NX) = 0;
    
    %Atualizar u e t
    t = t + dt;
    u_0 = u_num;
    
    %Solu��o anal�tica  
    trem=-alfa*pi^2*0.5;
    u_ana = 100*exp(trem)*sin(pi*x);
 
end

    %plotando os resultados
    plot(u_num,'-b');
     hold on
     plot(u_ana,'-r');
     hold off
    title('Expl�cito Simples: r = 2.0');
    grid on;
    xlabel('N�mero de n�s');
    ylabel('u(x,t)');
    legend('Sol. num�rica','Sol anal�tica');

