clc
clear all

%% Equação do calor - Euler explícito

%Parâmetros (constantes)
alfa = 0.2;         %condutividade
L = 1;              %comprimento do domínio
t_f = 0.5;          %tempo final

%Variáveis
NX = 11;            %número de pontos da malha (11 ou 16)
dx = L / (NX-1);    %passo de espaço 
t = 0;              %inicialização da variável "tempo"
r = 2.0;           %parâmetro de estabilidade (0.25, 0.50, 1, 2)
dt = r*(dx^2)/alfa; %passo de tempo

%Discretização do domínio
i = 1 : NX;
x = (i - 1)*dx;     %domínio espacial

%Condição inicial
u_0 = 100 * sin(pi * x / L);
u_num = u_0;        %inicialização de u_num na condição inicial

NT = t_f / dt;      %número do loop do tempo
%Implementação do método
for k = 2 : NT %loop no tempo
    
    for j = 2 : NX-1
        u_num(j) = u_0(j) + r*(u_0(j+1) - 2*u_0(j) + u_0(j-1)); 
    end
    
    %Condição de contorno
    u_num(1) = 0;
    u_num(NX) = 0;
    
    %Atualizar u e t
    t = t + dt;
    u_0 = u_num;
    
    %Solução analítica  
    trem=-alfa*pi^2*0.5;
    u_ana = 100*exp(trem)*sin(pi*x);
 
end

    %plotando os resultados
    plot(u_num,'-b');
     hold on
     plot(u_ana,'-r');
     hold off
    title('Explícito Simples: r = 2.0');
    grid on;
    xlabel('Número de nós');
    ylabel('u(x,t)');
    legend('Sol. numérica','Sol analítica');

