%% Exercicio 2c
%Reset do ambiente de trabalho
close all;

%Iteracoes
ttotal = 250;

%Vectores para o plot3
state = repmat([1:20], ttotal, 1);
t = repmat(linspace(0, ttotal, ttotal), 20, 1);

%Matriz para guardar probabilidades
prob = zeros(20, ttotal);

%Indices para testar como estado inicial
i_test = [1 2 4 6 8 10];

%Calcula numero de indices individuais a ser testados
[~, i_size] = size(i_test); 

%Matriz de varios conjuntos de probabilidades a serem testadas
prob0 = zeros(20, i_size+2);

%Define indices a serem alterados para estados iniciais
for n = 1:i_size
    prob0(i_test(n), n) = 1;
end

%Caso de probabilidade inicial com distribuicao uniforme
prob0(:, i_size+1) = 1/20;

%Caso de probabilidade inicial com distribuicao de equilibrio
prob0(:, i_size+2) = v_norm;

for n = 1:i_size+2
    prob(:, 1) = prob0(:, n);

    %Simulacao para condicoes iniciais escolhidas
    for i = 2:ttotal
        prob(:, i) = prob(:, i-1)'*P;
    end

    figure;
    plot3(t', state, prob);
    xlabel('Tempo [s]');
    ylabel('Estados');
    zlabel('Probabilidades');
    if(n <= i_size)
        title(sprintf('Evolucao das probabilidades para o estado inicial: %d', i_test(n)));
    elseif(n == i_size+1)
        title('Evolucao das probabilidades para distribuicao uniforme');
    elseif(n == i_size+2)
        title('Evolucao das probabilidades para distribuicao de equilibrio');
    end
    grid on;
end

%%
% *Comentarios:*
% Atraves da analise dos graficos podemos concluir que dado um tempo suficientemente grande as probabilidades de cada estado tendem para um certo valor limite,
% que são aproximadamente iguais aos valores obtidos na pergunta 2a