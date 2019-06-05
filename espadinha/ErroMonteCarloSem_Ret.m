function step_erro = ErroMonteCarloSem_Ret(runs, steps, Lambda, P_Fonte, desvio_pad, QP, fonte, nodePos, P)


% Matriz que vai conter a estimativa do erro para cada run e instante de
% tempo
erro_matrix = zeros(steps, runs);
step_erro = zeros(1, steps);

% Inicializa barra de loading
progresso = 0;
f = waitbar(progresso, 'Starting');



%Possiveis Estados iniciais (Remover Subconjuntos 2 e 3 [2d])
possivel_estado_inicial = [5, 6, 8, 9, 10, 11, 12, 15, 17];

for i=1:runs % realiza as M runs
    
    pos=zeros(steps, 2); % Inicializa vetor de posições a 0 para cada run
    
    ind = randi(length(possivel_estado_inicial));
    estado_inicial = possivel_estado_inicial(ind); % Gera um estado inicial aleatório
    
    estado_atual = estado_inicial;
    
    % Inicializa RlsPar em cada run
    RlsPar = struct('lam',Lambda);
    
    
    %---------
    pos(1,:) = [nodePos(estado_atual,2) nodePos(estado_atual,3)];
        
    % estimativa do erro
    D = squareform(pdist([fonte zeros(size(fonte)) pos']'));
    d = D(1,3:end);				% distâncias fonte-âncora
    an = D(2,3:end);			% norma das âncoras
        
    Pot = P_Fonte./(d.^2);				% potências sem ruído
    Pot = Pot.*exp(desvio_pad*randn(size(Pot)));	% potência com ruído
    Pot = QP*round(Pot/QP);			% quantização das potências
        
    A = [-2*repmat(Pot,[2 1]).*pos'; -ones(size(Pot)); Pot]';
    b = (-Pot.*(an.^2))';
    
    % Estimação RLS incremental
    for k = 1:size(A,1)
        [~,w,RlsPar] = qrrls(A(k,:),b(k),RlsPar);
    end

    
    %----------------------
    erro_matrix(1, runs) = norm(fonte-w(1:2));
    
    for j = 2:steps % realiza N passos para um dado estado inicial
        
        estado_atual = find (cumsum(P(estado_atual, :))>rand, 1 , 'first');
        

        
        % matriz de posições das âncoras vai aumentar de tamanho com cada
        % iteração
        pos(j,:) = [nodePos(estado_atual,2) nodePos(estado_atual,3)];
            
        D = squareform(pdist([fonte zeros(size(fonte)) pos']'));
        d = D(1,3:end);				
        an = D(2,3:end);			
        
        Pot = P_Fonte./(d.^2);			
        Pot = Pot.*exp(desvio_pad*randn(size(Pot)));
        Pot = QP*round(Pot/QP);			
        
        A = [-2*repmat(Pot,[2 1]).*pos'; -ones(size(Pot)); Pot]';
        b = (-Pot.*(an.^2))';
        
        for k = 1:size(A,1)
            [~,w,RlsPar] = qrrls(A(k,:),b(k),RlsPar);
        end
        
        % Matriz de erro vai conter a norma entre a posição da fonte e da
        % estimação atual
        erro_matrix(j, runs) = norm(fonte-w(1:2));
        
        % Incrementa o progresso da loading bar
        progresso = progresso+1;
        waitbar (progresso/(runs*steps), f, sprintf('Processing'));      
    end
end

% Fecha loading bar
delete(f);

for i = 1:steps
    step_erro(i) = mean(erro_matrix(i,:));
end

end