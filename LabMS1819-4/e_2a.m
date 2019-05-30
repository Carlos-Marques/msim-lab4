%% Exercicio 2a
%Reset do ambiente de trabalho
clear;
close all;

load MarkovChain;

[v, u] = eig('P');

[~, i] = min(abs(u(:)-1));

i = mod(i, size(u, 1));

v_norm = v(:, i) / sum(v(:, i));

figure();
bar(v_norm);
