clc
clear all

%% Equa��o do calor - Crank-Nicolson

%Par�metros (constantes)
alfa = 0.2;         %
L = 1;              %comprimento do dom�nio
t_f = 0.5;          %tempo final

%Vari�veis
NX = 11;            %n�mero de pontos da malha (11 ou 16)
dx = L / (NX-1);    %passo de espa�o 
r = 0.25;           %par�metro de estabilidade (0.25, 0.50, 1, 2)
dt = r*(dx^2)/alfa; %passo de tempo
t = 0;              %inicializa��o da vari�vel "tempo" em zero

%Discretiza��o do dom�nio
i = 1 : NX;
x = (i - 1)*dx;     %dom�nio espacial
NT = t_f / dt;      %n�mero do loop do tempo

%Condi��o inicial
u_0 = 100 * sin(pi * x / L);

%Condi��o de contorno
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
        %d � um vetor linha
    end
    
    %transformar os escalares a, b e c em vetores coluna
    bb = b*ones(NX-3,1);
    cc = bb;
    aa = a*ones(NX-2,1);

    %Construindo a matriz tridiagonal
    AA = diag(aa)+diag(-bb,1)+diag(-cc,-1);

    %Encontrar a solu��o para os n�s internos
    UU = AA\d';

    %Solu��o completa que inclui as C.C.
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




