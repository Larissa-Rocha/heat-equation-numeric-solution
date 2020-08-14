clc
clear all

%% Equação do calor - Crank-Nicolson

%Parâmetros (constantes)
alfa = 0.2;         %
L = 1;              %comprimento do domínio
t_f = 0.5;          %tempo final

%Variáveis
NX = 11;            %número de pontos da malha (11 ou 16)
dx = L / (NX-1);    %passo de espaço 
r = 0.25;           %parâmetro de estabilidade (0.25, 0.50, 1, 2)
dt = r*(dx^2)/alfa; %passo de tempo
t = 0;              %inicialização da variável "tempo" em zero

%Discretização do domínio
i = 1 : NX;
x = (i - 1)*dx;     %domínio espacial
NT = t_f / dt;      %número do loop do tempo

%Condição inicial
u_0 = 100 * sin(pi * x / L);

%Condição de contorno
u_n(1) = 0;
u_n(NX) = 0;

%Armazenando os resltados da C.I. para uso futuro
UUU(1,:) = u_0;

%Coeficientes da matriz tridiagonal
b = r;              %diagonal superior
c = r;              %diagonal inferior
a = 2 + 2*r;        %diagonal principal

%Loop no tempo
for k = 2 : NT
    for j = 1 : NX-2
        if j==1
            d(j) = r*u_0(j) + (2-2*r)*u_0(j+1) + r*u_0(j+2) + r*u_n(1);
        elseif j==NX-2
            d(j) = r*u_0(j) + (2-2*r)*u_0(j+1) + r*u_0(j+2) + r*u_n(NX);
        else
            d(j) = r*u_0(j) + (2-2*r)*u_0(j+1) + r*u_0(j+2);
        end
        %d é um vetor linha
    end
    
    %transformar os escalares a, b e c em vetores coluna
    bb = b*ones(NX-3,1);
    cc = bb;
    aa = a*ones(NX-2,1);

    %Construindo a matriz tridiagonal
    AA = diag(aa)+diag(-bb,1)+diag(-cc,-1);

    %Encontrar a solução para os nós internos
    UU = AA\d';

    %Solução completa que inclui as C.C.
    u_n = [u_n(1),UU',u_n(NX)];

    %Armazenar resultados para uso futuro
    UUU(k,:)=u_n;
    
    %Atualizando u e t
    u_0 = u_n;
    t = t + dt;
  
    plot(u_0,'-r');
    hold on
    plot(u_n,'bo');
    hold off
    pause (dt);

end




