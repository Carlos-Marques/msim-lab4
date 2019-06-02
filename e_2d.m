%% Exercicio 2d
%Reset do ambiente de trabalho
close all;

MarkovChainDraw();

% 2d i)
close all;

%Alteracao dos pesos das ligacoes de maneira a melhorar a situacao
%- diminuicao dos pesos das ligacoes que formam os subgrupos onde o token passa mais tempo
%- aumento dos pesos das ligacoes que formam os subgrupos onde o token passa menos tempo
Pi = [ 
    [1, 6, 0.3]; 
    [1, 7, 0.35]; 
    [1, 20, 0.35]; 
    [3, 12, 0.5]; 
    [3, 19, 0.5]; 
    [6, 1, 0.3]; 
    [6, 11, 0.35]; 
    [6, 15, 0.35]; 
    [12, 3, 0.3]; 
    [12, 8, 0.35]; 
    [12, 10, 0.35] 
];

Pc = P;

Pi_length = size(Pi, 1);

%Alteracao da matriz de pesos dada
for n=1:Pi_length
    Pc(Pi(n, 1), Pi(n, 2)) = Pi(n, 3);
end

%Calcula os vectores e valores proprios da matriz Pc transposta
[v_c, u_c] = eig(Pc');

%Encontra o indice do valor proprio 1
[~, i_c] = min(abs(u_c(:) - 1));

i_c = mod(i_c, size(u_c, 1));

%Normaliza o vector
v_c_norm = v_c(:, i_c) / sum(v_c(:, i_c));

figure;
bar([v_norm, v_c_norm]);
grid on;
xlabel('Estados');
ylabel('Probabilidade')

%Potencia da fonte
Pw0 = 100;
%Desvio padrao
sig = 10^-1;

%Numero de medidas
M = 1000;

No = round(v_norm.*M);
a = zeros(sum(No), 2);
k1 = 1;
k2 = 0;

%Cria observacoes para cada ancora
for i=1:size(No)
  k2 = k2+No(i);
  a(k1:k2,:) = repmat([nodePos(i,2) nodePos(i, 3)], No(i), 1);
  k1 = k1 + No(i);
end

%Obtem a posicao da fonte
x = sourcePos';

D = squareform(pdist([x zeros(size(x)) a']'));

%Calcula distancia entre fonte e a ancora 
d = D(1, 3:end);
%Calcula normas das acoras
an = D(2, 3:end);

%Calcula potencia nas ancoras sem ruido
Pw = Pw0 ./ (d.^2);
%Aplica ruido
Pw = Pw.*exp(sig*rand(size(Pw)));
QPw = 1e-2;
%Quantitiza as potencias
Pw = QPw*round(Pw/QPw);

%Aplica metodo dos minimos quadrados
A = [-2*repmat(Pw, [2 1]).*a'; -ones(size(Pw)); Pw]';
b = (-Pw.*(an.^2))';

z = A\b;
xe = z(1:2);
fprintf('Distancia entre posicao real e calculada da fonte: %f\n', norm(x-xe));

figure;

plot(a*[1; 1i],'o', 'DisplayName', 'Ancoras'); 
hold all;
grid on;
plot(x'*[1; 1i],'x', 'DisplayName', 'Real'); 
plot(xe'*[1; 1i],'s', 'DisplayName', 'Calculada');
axis('square')
title('Estimativa da posicao da fonte');
legend('Location', 'northeastoutside');

% RLS formulation (one-shot)
RlsPar = struct('lam',1);
[e,w,RlsPar] = qrrls(A,b,RlsPar);
fprintf('Erro da Recursive Least Squares (one-shot): %f\n', norm(z-w));

% RLS formulation (incremental)
RlsPar = struct('lam',1);
for i = 1:size(A,1)
  [e,w,RlsPar] = qrrls(A(i,:),b(i),RlsPar);
end
fprintf('Erro da Recursive Least Squares (incremental): %f\n', norm(z-w));

state = repmat([1:20], ttotal, 1);
t = repmat(linspace(0, ttotal, ttotal), 20, 1);
prob = zeros(20, ttotal);

for n = 1:i_size+2
    prob(:, 1) = prob0(:, n);

    for i = 2:ttotal
        prob(:, i) = prob(:, i-1)'*P;
    end

    figure;
    plot3(t', state, prob);
    grid on;
end

% 2d ii)
close all;

Pi = [ 
    [3, 12, 0]; 
    [3, 19, 1]; 
];

Pc = P;

Pi_length = size(Pi, 1);

for n=1:Pi_length
    Pc(Pi(n, 1), Pi(n, 2)) = Pi(n, 3);
end

[v_c, u_c] = eig(Pc');

[~, i_c] = min(abs(u_c(:) - 1));

i_c = mod(i_c, size(u_c, 1));

v_c_norm = v_c(:, i_c) / sum(v_c(:, i_c));

figure;
bar([v_norm, v_c_norm]);

Pw0 = 100;
sig = 10^-2;

a = nodePos(:, [2 3])';
x = sourcePos';

D = squareform(pdist([x zeros(size(x)) a]'));
d = D(1, 3:end);
an = D(2, 3:end);

Pw = Pw0 ./ (d.^2);
Pw = Pw.*exp(sig*rand(size(Pw)));
QPw = 1e-2;
Pw = QPw*round(Pw/QPw);

A = [-2*repmat(Pw, [2 1]).*a; -ones(size(Pw)); Pw]';
b = (-Pw.*(an.^2))';

v_p = diag(v_c_norm);

A = v_p*A;
b = v_p*b;

z = A\b;
xe = z(1:2);
err = norm(x-xe);

figure();
hold on
plot(a'*[1; 1i], 'o');
plot(x'*[1; 1i],'x'); 
plot(xe'*[1; 1i],'s');

prob = zeros(20, ttotal);

for n = 1:i_size+2
    prob(:, 1) = prob0(:, n);

    for i = 2:ttotal
        prob(:, i) = prob(:, i-1)'*P;
    end

    figure;
    plot3(t', state, prob);
    grid on;
end