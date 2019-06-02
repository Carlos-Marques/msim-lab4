%% Exercicio 2a
%Reset do ambiente de trabalho
clear;
close all;

load MarkovChain;

%Calcula os vectores e valores proprios da matriz P transposta
[v, u] = eig(P');

%Encontra o indice do valor proprio 1
[~, i] = min(abs(u(:)-1));

i = mod(i, size(u, 1));

%Normaliza o vector
v_norm = v(:, i) / sum(v(:, i));

figure();
bar(v(:, i), 'DisplayName', 'nao normalizado');
hold on;
grid on;
bar(v_norm , 'DisplayName', 'normalizado');
title('Probabilidades Limide da Cadeia de Markov');
xlabel('Estados');
ylabel('Probabilidade');
legend('Location', 'northeastoutside');

%%
% *Comentarios:*
% Observa-se que os estados com maior probabilidade sao o estado 7 e 19 
% (9.649%) e que os estados com menor probabilidade são o 8 e o 17 
% (1.072%). Tal não é surpreendente visto que 7 e 19 
% pertencem a um subgrupo onde o token customa ficar preso e 
% 8 e 17 pertencem a um subgrupo onde o token não tende a ficar.