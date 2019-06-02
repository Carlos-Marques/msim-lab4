function [z,w,Par] = qrrls(U,d,Par)
%
% Filtro QR-RLS transversal.
%
% Simon Haykin, "Adaptive Filter Theory", 3rd ed., Prentice-Hall, 1996
% (Tabela 14.2)
%
% U - Matriz com vectores de entrada (um por linha)
% d - Vector de referências (uma por linha)
%
% Par - Parâmetros iniciais/finais do algoritmo RLS
%     lam - Factor de esquecimento \in (0,1]
%     p - Vector de coeficientes modificado
%     F - Matriz de covariância (factor de Cholesky triangular inferior)
%
% z - Erro 'a posteriori'/'a priori'
% w - Vector de coeficientes do filtro transversal convencional
%


%---------------------------------------------------------------------------
%
% Inicialização das matrizes e vectores
%
%---------------------------------------------------------------------------
					
[M,N] = size(U);

if ~isfield(Par,'lam'), Par.lam = 1; end
slam = sqrt(Par.lam);

if ~isfield(Par,'p'), Par.p = zeros(N,1); end
if ~isfield(Par,'F'), Par.F = diag(repmat(1e-4,[N 1])); end
p = Par.p; F = Par.F;

z = zeros(size(U,1),1);


%---------------------------------------------------------------------------
%
% Algoritmo QR-RLS
%
%---------------------------------------------------------------------------

for i = 1:M
  
  u = U(i,:).';
  %-------------------------------------------------------------------------
  % Inicializar pré-array
  X = [slam*F u; slam*p' d(i); zeros(size(p')) 1];

  %-------------------------------------------------------------------------
  % Filtragem
  for k = 1:N
    L = rotg([X(k:end,k) X(k:end,end)]);
    X(k:end,k) = L(:,1); X(k:end,end) = L(:,2);
  end

  %-------------------------------------------------------------------------
  % Extrair variáveis do pós-array
  p = X(N+1,1:N)';
  F = X(1:N,1:N);
  e = X(N+1,N+1);			% Angle-normalized error
  g = X(end,end);			% Conversion factor

  z(i) = e/g;				% A priori error
%   z(i) = e*g;				% A posteriori error
end


%---------------------------------------------------------------------------
%
% Finalização (housekeeping)
%
%---------------------------------------------------------------------------

w = (p'/F)';
Par.p = p; Par.F = F;


function [L,U] = rotg(X,U)
%
% Aniquilação de um elemento de uma matriz (complexa) usando rotações
%

if nargin == 1,  U = givens(X(1,1),X(1,2)).'; end
L = X*U;
