%Função para resolver a matriz tridiagonal
%Input:
%P = vetor subdiagonal
%Q = vetor diagonal
%R = vetor superdiagonal
%S = vetor da direita
%Output:
%u = vetor solução

function u = Tridiag(P,Q,R,S)

n = numel(Q);

%eliminação para frente
for k = 2 : n
    factor = P(k)/Q(k-1);
    Q(k) = Q(k)-factor*R(k-1);
    S(k) = S(k) - factor*S(k-1);
end

%back substitution
u(n) = S(n)/Q(n);
for k = n-1 : -1 : 1
    u(k) = (S(k) - R(k)*u(k+1)) / Q(k);
end

end