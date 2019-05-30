%% Exercicio 2a
%Reset do ambiente de trabalho
clear;
close all;

load MarkovChain;

[v, u] = eig('P');

[~, i] = min(abs(U(:)-1));

i = mod(i, size(U, 1));

v_norm = v(:, i) / sum(v(:, i));

figure();
bar(v_norm);
