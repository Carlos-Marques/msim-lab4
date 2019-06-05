function [ocorrencias_estados, caminho_estados, step_erro] = MonteCarlo(runs, steps, Lambda, P_Fonte, desvio_pad, QP, fonte, nodePos, P)

% Realizar M runs com N passos
% Gera matriz que vai conter os estados resultantes das simula��es de Monte Carlo
caminho_estados = zeros(runs,steps);  
ocorrencias_estados = zeros (1, 20);

% Matriz que vai conter a estimativa do erro para cada run e instante de
% tempo

erro_matrix = zeros(steps, runs);
step_erro = zeros(1, steps);

% Inicializa barra de loading
progresso = 0;
f = waitbar(progresso, 'Starting');

for i=1:runs % realiza as M runs
    
    pos=zeros(steps, 2); % Inicializa vetor de posi��es a 0 para cada run
    
    estado_inicial = randi([1,20]); % Gera um estado inicial aleat�rio
    
    caminho_estados (i, 1) = estado_inicial;
    
    estado_atual = estado_inicial;
    
    % Inicializa RlsPar em cada run
    RlsPar = struct('lam',Lambda);
    
    
    %---------
    pos(1,:) = [nodePos(estado_atual,2) nodePos(estado_atual,3)];
        
    % estimativa do erro
    D = squareform(pdist([fonte zeros(size(fonte)) pos']'));
    d = D(1,3:end);				% dist�ncias fonte-�ncora
    an = D(2,3:end);			% norma das �ncoras
        
    Pot = P_Fonte./(d.^2);				% pot�ncias sem ru�do
    Pot = Pot.*exp(desvio_pad*randn(size(Pot)));	% pot�ncia com ru�do
    Pot = QP*round(Pot/QP);			% quantiza��o das pot�ncias
        
    A = [-2*repmat(Pot,[2 1]).*pos'; -ones(size(Pot)); Pot]';
    b = (-Pot.*(an.^2))';
    
    % Estima��o RLS incremental
    for k = 1:size(A,1)
        [~,w,RlsPar] = qrrls(A(k,:),b(k),RlsPar);
    end

    % Vetor que conta as ocorr�ncias de cada estado
    ocorrencias_estados (estado_atual) = ocorrencias_estados (estado_atual) + 1;
    
    %----------------------
    erro_matrix(1, runs) = norm(fonte-w(1:2));
    
    for j = 2:steps % realiza N passos para um dado estado inicial
        
        estado_atual = find (cumsum(P(estado_atual, :))>rand, 1 , 'first');
        
        caminho_estados (i, j) = estado_atual;
        
        % matriz de posi��es das �ncoras vai aumentar de tamanho com cada
        % itera��o
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
        
        % Vetor que conta as ocorr�ncias de cada estado
        ocorrencias_estados (estado_atual) = ocorrencias_estados (estado_atual) + 1;
        
        % Matriz de erro vai conter a norma entre a posi��o da fonte e da
        % estima��o atual
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

