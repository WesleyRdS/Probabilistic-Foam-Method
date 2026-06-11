clc; clear; close all;


%% ==================== PFM_3D_COMPLETE.M ====================
% Implementação dos 4 algoritmos do artigo em 3D:
% - PFM (Probabilistic Foam Method)
% - GBPF (Goal-Biased Probabilistic Foam)
% - RBPF (Radius-Biased Probabilistic Foam)
% - HPF (Heuristic-Guided Probabilistic Foam)
%
% Baseado em: Sensors 2021, 21, 4156
% Map 4 - Complex 3D Environment

%% ==================== ALGORITMO 1: PFM 3D ====================
function [rosary, success, stats] = ...
    PFM_STL(q_init, q_goal, expandir_bolha, ...
            r_min, K, z_minimo, z_maximo)
    % Algoritmo 1 do artigo: Probabilistic Foam Method (versão 3D)
    
    tic;
    
    % Inicialização
    F = {};           % Lista de bolhas
    Q = {};           % Fila de índices para expansão
    
    % Bolha inicial
    r_init = max(expandir_bolha(q_init),r_min);
    r_init = min(r_init, q_init(3) - z_minimo);
    F{1} = struct('center', q_init, 'radius', r_init, 'parent', 0);
    Q{1} = 1;
    
    max_iter = 3000;
    success = false;
    goal_bubble_idx = [];
    
    for iter = 1:max_iter
        if isempty(Q); break; end
        
        parent_idx = Q{1};
        Q(1) = [];
        parent = F{parent_idx};
        
        % Número máximo de bolhas filhas (Equação 2)
        N = K * floor(parent.radius / r_min);
        N = max(1, min(N, 30));
        
        for n = 1:N
            % Direção aleatória uniforme na esfera
            direction = randn(1, 3);
            if norm(direction) < 1e-6; continue; end
            direction = direction / norm(direction);
            
            q_surface = parent.center + direction * parent.radius;
            
            if q_surface(3) < z_minimo || ...
               q_surface(3) > z_maximo
                continue;
            end
            
            % Verifica se está dentro de outra bolha
            inside = false;
            for i = 1:length(F)
                if norm(q_surface - F{i}.center) < F{i}.radius - 0.1
                    inside = true;
                    break;
                end
            end
            if inside; continue; end
            
            % Expande nova bolha
            r_new = expandir_bolha(q_surface);
            r_max_permitido = q_surface(3) - z_minimo;
            r_new = min(r_new, r_max_permitido);
            if r_new >= r_min
                new_idx = length(F) + 1;
                F{new_idx} = struct('center', q_surface, 'radius', r_new, ...
                                    'parent', parent_idx);
                Q{end+1} = new_idx;
                
                % Verifica goal
                if norm(q_surface - q_goal) <= r_new
                    success = true;
                    goal_bubble_idx = new_idx;
                    break;
                end
            end
        end
        
        if success; break; end
    end
    
    if success
        rosary = extract_rosary_3d(F, goal_bubble_idx);
        stats.path_length = compute_path_length_3d(rosary);
    else
        rosary = {};
        stats.path_length = inf;
    end
    
    stats.time = toc;
    stats.num_bubbles = length(F);
    stats.success = success;
    stats.max_iter = iter;
end

%% ==================== ALGORITMO 2: GBPF 3D ====================
function [rosary,success,stats] = ...
GBPF_STL(...
    q_init,...
    q_goal,...
    expandir_bolha,...
    r_min,...
    bias,...
    z_minimo,...
    z_maximo,...
    vertices)
    % Algoritmo 2 do artigo: Goal-Biased Probabilistic Foam (3D)
    
    tic;
    
    F = {};
    
    % Bolha inicial
    r_init = max(expandir_bolha(q_init),r_min);
    r_init = min(r_init, q_init(3) - z_minimo);
    F{1} = struct('center', q_init, 'radius', r_init, 'parent', 0);
    
    max_iter = 3000;
    success = false;
    goal_bubble_idx = [];
    
    for iter = 1:max_iter
        % Amostra q_aux com viés (lines 5-9 do Algorithm 2)
        if rand() < bias
            q_aux = q_goal;
        else
            xmin = min(vertices(:,1));
            xmax = max(vertices(:,1));
            
            ymin = min(vertices(:,2));
            ymax = max(vertices(:,2));
            
            q_aux = [ ...
                xmin + rand*(xmax-xmin), ...
                ymin + rand*(ymax-ymin), ...
                z_minimo + rand*(z_maximo-z_minimo)];
        end
        
        % Encontra bolha mais próxima (nearest_bubble)
        nearest_idx = 1;
        nearest_dist = norm(F{1}.center - q_aux);
        for i = 2:length(F)
            dist = norm(F{i}.center - q_aux);
            if dist < nearest_dist
                nearest_dist = dist;
                nearest_idx = i;
            end
        end
        
        parent = F{nearest_idx};
        
        % Encontra ponto na superfície mais próximo (nearest_config)
        direction = q_aux - parent.center;
        if norm(direction) > 0
            direction = direction / norm(direction);
        else
            direction = randn(1, 3);
            direction = direction / norm(direction);
        end
    
        q_near = parent.center + direction * parent.radius;
        
        if q_near(3) < z_minimo || ...
           q_near(3) > z_maximo
            continue;
        end
        
        % Verifica interior
        inside = false;
        for i = 1:length(F)
            if norm(q_near - F{i}.center) < F{i}.radius - 0.1
                inside = true;
                break;
            end
        end
        if inside; continue; end
        
        % Expande nova bolha
        r_new = expandir_bolha(q_near);
        r_max_permitido = q_near(3) - z_minimo;
        r_new = min(r_new, r_max_permitido);
        if r_new >= r_min
            new_idx = length(F) + 1;
            F{new_idx} = struct('center', q_near, 'radius', r_new, ...
                                'parent', nearest_idx);
            
            if norm(q_near - q_goal) <= r_new
                success = true;
                goal_bubble_idx = new_idx;
                break;
            end
        end
    end
    
    if success
        rosary = extract_rosary_3d(F, goal_bubble_idx);
        stats.path_length = compute_path_length_3d(rosary);
    else
        rosary = {};
        stats.path_length = inf;
    end
    
    stats.time = toc;
    stats.num_bubbles = length(F);
    stats.success = success;
    stats.max_iter = iter;
end

%% ==================== ALGORITMO 3: RBPF 3D ====================
function [rosary, success, stats] =     RBPF_STL(...
    q_init,...
    q_goal,...
    expandir_bolha,...
    r_min,...
    K,...
    z_minimo,...
    z_maximo)
    % Algoritmo 3 do artigo: Radius-Biased Probabilistic Foam (3D)
    % Roulette Wheel baseado no raio (Equação 6)
    
    tic;
    
    F = {};
    OpenList = [];
    
    % Bolha inicial
    r_init = max(expandir_bolha(q_init),r_min);
    r_init = min(r_init, q_init(3) - z_minimo);
    F{1} = struct('center', q_init, 'radius', r_init, 'parent', 0);
    OpenList = 1;
    
    max_iter = 3000;
    success = false;
    goal_bubble_idx = [];
    
    for iter = 1:max_iter
        if isempty(OpenList); break; end
        
        % Seleção por roleta (Equação 6)
        radii = zeros(length(OpenList), 1);
        for i = 1:length(OpenList)
            radii(i) = F{OpenList(i)}.radius;
        end
        
        if sum(radii) == 0; break; end
        prob = radii / sum(radii);
        
        % Escolhe baseado na probabilidade
        r = rand();
        cumsum_prob = cumsum(prob);
        chosen_idx = find(r <= cumsum_prob, 1, 'first');
        parent_idx = OpenList(chosen_idx);
        parent = F{parent_idx};
        
        % Remove da OpenList
        OpenList(chosen_idx) = [];
        
        % Número máximo de bolhas filhas
        N = K * floor(parent.radius / r_min);
        N = max(1, min(N, 30));
        
        for n = 1:N
            direction = randn(1, 3);
            if norm(direction) < 1e-6; continue; end
            direction = direction / norm(direction);

            q_surface = parent.center + direction * parent.radius;
            
            if q_surface(3) < z_minimo || ...
               q_surface(3) > z_maximo
                continue;
            end
            
            % Verifica cobertura
            covered = false;
            for i = 1:length(F)
                if norm(q_surface - F{i}.center) < F{i}.radius - 0.1
                    covered = true;
                    break;
                end
            end
            if covered; continue; end
            
            r_new = expandir_bolha(q_surface);
            r_max_permitido = q_surface(3) - z_minimo;
            r_new = min(r_new, r_max_permitido);
            if r_new >= r_min
                new_idx = length(F) + 1;
                F{new_idx} = struct('center', q_surface, 'radius', r_new, ...
                                    'parent', parent_idx);
                OpenList(end+1) = new_idx;
                
                if norm(q_surface - q_goal) <= r_new
                    success = true;
                    goal_bubble_idx = new_idx;
                    break;
                end
            end
        end
        
        if success; break; end
    end
    
    if success
        rosary = extract_rosary_3d(F, goal_bubble_idx);
        stats.path_length = compute_path_length_3d(rosary);
    else
        rosary = {};
        stats.path_length = inf;
    end
    
    stats.time = toc;
    stats.num_bubbles = length(F);
    stats.success = success;
    stats.max_iter = iter;
end

%% ==================== ALGORITMO 4: HPF 3D ====================
function [rosary, success, stats] = HPF_STL(...
    q_init,...
    q_goal,...
    expandir_bolha,...
    r_min,...
    K,...
    z_minimo,...
    z_maximo)
    % Algoritmo 4 do artigo: Heuristic-Guided Probabilistic Foam (3D)
    % f(q) = g(q) + h(q) (Equação 7)
    
    tic;
    
    F = {};
    OpenList = struct('idx', {}, 'g', {}, 'h', {}, 'f', {});
    
    % Bolha inicial
    r_init = max(expandir_bolha(q_init),r_min);
    r_init = min(r_init, q_init(3) - z_minimo);
    F{1} = struct('center', q_init, 'radius', r_init, 'parent', 0);
    
    g_init = r_init;
    h_init = norm(q_init - q_goal);
    f_init = g_init + h_init;
    
    OpenList(1).idx = 1;
    OpenList(1).g = g_init;
    OpenList(1).h = h_init;
    OpenList(1).f = f_init;
    
    max_iter = 3000;
    success = false;
    goal_bubble_idx = [];
    
    for iter = 1:max_iter
        if isempty(OpenList); break; end
        
        % Seleciona bolha com menor f
        f_values = [OpenList.f];
        [~, min_idx] = min(f_values);
        current = OpenList(min_idx);
        OpenList(min_idx) = [];
        
        parent = F{current.idx};
        
        % Número máximo de bolhas filhas
        N = K * floor(parent.radius / r_min);
        N = max(1, min(N, 30));
        
        for n = 1:N
            direction = randn(1, 3);
            if norm(direction) < 1e-6; continue; end
            direction = direction / norm(direction);

            
            q_surface = parent.center + direction * parent.radius;
            
            if q_surface(3) < z_minimo || ...
               q_surface(3) > z_maximo
                continue;
            end
            
            % Verifica cobertura
            covered = false;
            for i = 1:length(F)
                if norm(q_surface - F{i}.center) < F{i}.radius - 0.1
                    covered = true;
                    break;
                end
            end
            if covered; continue; end
            
            r_new = expandir_bolha(q_surface);
            r_max_permitido = q_surface(3) - z_minimo;
            r_new = min(r_new, r_max_permitido);
            if r_new >= r_min
                new_idx = length(F) + 1;
                F{new_idx} = struct('center', q_surface, 'radius', r_new, ...
                                    'parent', current.idx);
                
                % Calcula custos (Equação 7)
                g_new = current.g + r_new;
                h_new = norm(q_surface - q_goal);
                f_new = g_new + h_new;
                
                OpenList(end+1).idx = new_idx;
                OpenList(end).g = g_new;
                OpenList(end).h = h_new;
                OpenList(end).f = f_new;
                
                if norm(q_surface - q_goal) <= r_new
                    success = true;
                    goal_bubble_idx = new_idx;
                    break;
                end
            end
        end
        
        if success; break; end
    end
    
    if success
        rosary = extract_rosary_3d(F, goal_bubble_idx);
        stats.path_length = compute_path_length_3d(rosary);
    else
        rosary = {};
        stats.path_length = inf;
    end
    
    stats.time = toc;
    stats.num_bubbles = length(F);
    stats.success = success;
    stats.max_iter = iter;
end

%% ==================== FUNÇÕES AUXILIARES 3D ====================

function rosary = extract_rosary_3d(F, goal_idx)
    % Extrai o rosário (sequência de bolhas)
    rosary = {};
    current_idx = goal_idx;
    while current_idx > 0
        rosary{end+1} = F{current_idx};
        current_idx = F{current_idx}.parent;
    end
    rosary = fliplr(rosary);
end

function path_length = compute_path_length_3d(rosary)
    % Comprimento do caminho (distância entre centros)
    if length(rosary) < 2
        path_length = 0;
        return;
    end
    path_length = 0;
    for i = 1:length(rosary)-1
        path_length = path_length + norm(rosary{i+1}.center - rosary{i}.center);
    end
end

function sm = compute_safety_metric(rosary, r_min)
    % Safety Metric SM (Equação 8 do artigo)
    % SM = (1/k) * sum((r_i - r_min)^2)
    if isempty(rosary); sm = 0; return; end
    
    k = length(rosary);
    sum_sq = 0;
    for i = 1:k
        sum_sq = sum_sq + (rosary{i}.radius - r_min)^2;
    end
    sm = sum_sq / k;
end



%Função que calcula a menor distancia entre um ponto P e um triangulo ABC
%no espaço
function distancia = DistanciaPontoTriangulo(triangulo, ponto)
    %Vertices do triangulo
    vertice_A = triangulo(1, :);
    vertice_B = triangulo(2, :);
    vertice_C = triangulo(3, :);
    
    %Arestas do triangulo
    arestas_AB = vertice_B - vertice_A;
    arestas_AC = vertice_C - vertice_A;
    arestas_BC = vertice_C - vertice_B;
    
    %Vetor entre ponto e vertices do triangulo
    vetor_AP = ponto - vertice_A;
    vetor_BP = ponto - vertice_B;
    vetor_CP = ponto - vertice_C;

    %Ponto localizado na região do vertice A
    
    %Quando P avança na direção de AB a partir de A
    d1 = dot(arestas_AB, vetor_AP);
    %Quando P avança na direção de AC a partir de A
    d2 = dot(arestas_AC, vetor_AP);
    %Se P esta antes de A na direção AB
    % e P tambem esta antes de A na direção
    if d1 <= 0 && d2 <= 0
        distancia = norm(vetor_AP);
        return
    end

    %Ponto localizado na região do vertice B

    %Quando P avança na direção de AB a partir de B
    d3 = dot(arestas_AB, vetor_BP);
    %Quando P avança na direção de AC a partir de B
    d4 = dot(arestas_AC, vetor_BP);
    %Se P esta depois de B na direção AB
    % e se ponto não subiu na direção de AC mais do que avançou em AB
    if d3 >= 0 && d4 <= d3
        distancia = norm(vetor_BP);
        return
    end

    %Ponto localizado na região da Aresta AB
    
    %Verifica se P esta no mesmo lado da aresta AB que o triangulo
    %   vc > 0 -> Dentro do triangulo
    %   vc = 0 -> Sobre a aresta AB
    %   vc < 0 -> Fora, do outro lado
    vc = d1 * d4 - d3 * d2; %Determinante
    %Se P esta alinhado com a aresta AB, 
    %não esta antes de A e não passou por B
    if vc <= 0 && d1 >= 0 && d3 <= 0
        %Busca um v, tal que L(v) = A + v*AB
        %Calculo de peso quando ponto esta perto de B
        v = d1 / (d1 - d3);
        %Projeção de ponto na reta AB mais proximo de P
        proj_P = vertice_A + v * arestas_AB; %Representação da reta  
        distancia = norm(ponto - proj_P);
        return
    end

    %Ponto localizado na região do vertice C

    %Quando P avança na direção de AB a partir de C
    d5 = dot(arestas_AB, vetor_CP);
    %Quando P avança na direção de AC a partir de C
    d6 = dot(arestas_AC, vetor_CP);
    %Se P está "para frente" na direção de AC
    %e  P não está indo para dentro do triângulo
    if d6 >= 0 && d5 <= d6
        distancia = norm(vetor_CP);
        return
    end

    %Ponto localizado na região da aresta AC
    
    %P projeta na aresta AC
    vb = d5*d2 - d1*d6; %Determinante
    %Se P esta alinhado com a aresta AC, 
    %não esta antes de A e não passou por C
    if vb <= 0 && d2 >= 0 && d6 <= 0
        %Calculo de peso quando ponto esta perto de C
        w = d2 / (d2 - d6);
        proj_P = vertice_A + w * arestas_AC;
        distancia = norm(ponto - proj_P);
        return
    end

    %Ponto localizado na região da aresta BC
    
    va = d3 * d6 - d5 * d4;%Determinante
    %va <= 0 - Siginifica que P esta alinhado com BC ou fora do triangulo
    %(d4 - d3) >= 0 - P não esta voltado para A e esta indo na direção de C
    sentidoPpatindoB = d4 - d3;
    %(d5 - d6) >= 0 - Esta dentro da faixa da aresta BC
    sentidoPpatindoC = d5 - d6;
    if va <= 0 && sentidoPpatindoB >= 0 && sentidoPpatindoC >= 0
        w = sentidoPpatindoB / (sentidoPpatindoB + sentidoPpatindoC);
        proj_P = vertice_B + w * arestas_BC;
        distancia = norm(ponto - proj_P);
        return
    end

    %Ponto localizado no interior do triangulo

    normalizar = 1 / (va + vb + vc); 
    %Quando o ponto ta perto de B
    v = vb * normalizar;
    %Quando o ponto ta perto de C
    w = vc * normalizar;
    proj_P = vertice_A + arestas_AB * v + arestas_AC * w;
    distancia = norm(ponto - proj_P);
end

%Função otimizada para calcular a menor distancia só com os K triangulos mais
%proximos(Evita ter que calcular para todos os triangulos dos mapa)
function distancia = DistanciaKtriangulosProximos(ponto, vertices, faces, arvoreKdimensional, k)
    
    %Busca indice dos k triangulos mais proximos
    indice = knnsearch(arvoreKdimensional, ponto, 'K', k);

    %Definir menor distancia como infinito
    distancia = inf;

    for i = 1:length(indice)
        triangulo = vertices(faces(indice(i), :), :);
        distancia_auxiliar = DistanciaPontoTriangulo(triangulo, ponto);
        if distancia_auxiliar < distancia
            distancia = distancia_auxiliar;
        end
    end

end

%Função para calcular o raio da bolha 3D
function raio = calcular_raio(ponto, vertices, faces, arvoreKdimensional, k, z_minimo, z_maximo)
    %Se o ponto está fora da região vertical permitida, o raio é zero e
    %ele não é usado
    if ponto(3) < z_minimo || ponto(3) > z_maximo
        raio = 0;
        return
    end

    raio = DistanciaKtriangulosProximos(ponto, vertices, faces, arvoreKdimensional, k);
    %limitar raio máximo
    raio_maximo = 100;
    raio = min(raio, raio_maximo);

    %Garantia para não existir raio 0 
    raio = max(raio, 0.5);
end

function [latitude, longitude] = PosGPS(x, y, latitude_atual, longitude_atual)
    
    %O planeta é esferico então uso uma aproximação plana:
    %Circunferencia(metros)/360 = 40075000/360
    latitude_metros_por_grau = 111320; %111,32 km aproximadamente
    
    %C = 2*PI*R*cos(teta) -> cosd é o equivalente a cosseno em graus
    longitude_metros_por_grau = latitude_metros_por_grau * cosd(latitude_atual);
    
    latitude = latitude_atual + (y / latitude_metros_por_grau);
    longitude = longitude_atual + (x / longitude_metros_por_grau);
end

%Escala do mapa 3D do copelia no matlab
escala = 1; 

% Interface com RemoteApi do Coppelia
objetoAPI_remota = remApi('remoteApi');
% Fechar conexões antigas
objetoAPI_remota.simxFinish(-1);

id = objetoAPI_remota.simxStart('127.0.0.1', 19000, true, true, 5000, 5);

if id < 0
    disp('Falha ao tentar conectar o CoppeliaSim com o matlab')
    objetoAPI_remota.delete;
    return;
end

%Manipulador dos objetos na simulação
%Obs: No matlab ~ significa ignorar valor de saida

%Posição alvo do drone
[~, alvo] = objetoAPI_remota.simxGetObjectHandle(id, 'alvo', objetoAPI_remota.simx_opmode_blocking);
%Referencia de local de partida do drone
[~, ref_partida] = objetoAPI_remota.simxGetObjectHandle(id, 'partida', objetoAPI_remota.simx_opmode_blocking);
%Referencia do local de destino do drone
[~, ref_destino] = objetoAPI_remota.simxGetObjectHandle(id, 'destino', objetoAPI_remota.simx_opmode_blocking);

%Referência do mapa correto
[~, ref_mapa] = objetoAPI_remota.simxGetObjectHandle(id, 'UEFScompleto', objetoAPI_remota.simx_opmode_blocking);

%Posição do mapa (não do alvo)
[~, pos_mapa] = objetoAPI_remota.simxGetObjectPosition(id, ref_mapa, -1, objetoAPI_remota.simx_opmode_blocking);

%Obtendo a posição inicial do drone no ambiente 3D. Modo de operação
%bloqueante. O execução so continua após obeter retorno
[~, PosInicialDrone] = objetoAPI_remota.simxGetObjectPosition(id, alvo, -1, objetoAPI_remota.simx_opmode_blocking);

%Ler STL. As formas são aproximadas por triangulos
face_vertex = stlread('UEFScompleto.stl');

%Lista de vertices representados por pontos no espaço. Ex: Va = [Xa Ya Za]
vertices = face_vertex.Points;

%Lista de conectividade. Cada item representa um triangulo esses triangulos
%são fornecidos como um indice posicional de cada vertice na lista de
%pontos. Ex T1 = [Va, Vb, Vc]
faces = face_vertex.ConnectivityList;

%Calculo do centro geometrico dos triangulos a partir da media dos vertices
centroides = (vertices(faces(:, 1), :) + vertices(faces(:, 2), :) + vertices(faces(:, 3), :)) / 3;

%KD_Tree para busca rapida. 
arvoreKdimensional = KDTreeSearcher(centroides);

%Numero de vizinhos a serem pesquisados
k_vizinhos = 200;

%Limites do mapa 3D
z_minimo = 4;
z_maximo = 100;

%Margem segura de obstaculos
margem = 4;

expandir_bolha = @(ponto) calcular_raio(ponto, vertices, faces, arvoreKdimensional, k_vizinhos, z_minimo, z_maximo); 

% Interface 2D para selecionar destino do drone

figure;
trisurf(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
view(2)
axis equal
title("Selecione o local de destino")
hold on

ponto = (PosInicialDrone - pos_mapa) / escala;
partida = ponto; %Definindo posição inicial do drone como ponto de partida
%Plotando na representação superior 2D do mapa a localização do drone
plot(ponto(1), ponto(2), 'go', 'LineWidth', 2)

while true
    [x, y] = ginput(1); %Pega posição do mouse na figura ao clicar
    ponto = [x, y, ponto(3)];
    if expandir_bolha(ponto) > 0
        %Adicionando mais 8 metro de altura ao destino 
        destino = ponto + [0 0 8];
        plot(x, y, 'ro', 'LineWidth', 2)
        break
    end
end
hold off


%Aplicando os metodos

raio_minimo = 2;
K = 5;
bias = 0.05;

algorithms = {'PFM','GBPF','RBPF','HPF'};

time_vals = zeros(1,length(algorithms));
bubbles_vals = zeros(1,length(algorithms));
length_vals = zeros(1,length(algorithms));
sm_vals = zeros(1,length(algorithms));
success_vals = false(1,length(algorithms));

todos_rosarios = cell(1,length(algorithms));

for alg_idx = 1:length(algorithms)

    fprintf('\nExecutando %s\n',algorithms{alg_idx});

    switch algorithms{alg_idx}

        case 'PFM'

            [rosary,success,stats] = ...
                PFM_STL(...
                partida,...
                destino,...
                expandir_bolha,...
                raio_minimo,...
                K,...
                z_minimo,...
                z_maximo);

        case 'GBPF'

            [rosary,success,stats] = ...
                GBPF_STL(...
                partida,...
                destino,...
                expandir_bolha,...
                raio_minimo,...
                bias,...
                z_minimo,...
                z_maximo,...
                vertices);

        case 'RBPF'

            [rosary,success,stats] = ...
                RBPF_STL(...
                partida,...
                destino,...
                expandir_bolha,...
                raio_minimo,...
                K,...
                z_minimo,...
                z_maximo);

        case 'HPF'

            [rosary,success,stats] = ...
                HPF_STL(...
                partida,...
                destino,...
                expandir_bolha,...
                raio_minimo,...
                K,...
                z_minimo,...
                z_maximo);

    end

    success_vals(alg_idx)=success;

    if success

        sm = compute_safety_metric(rosary,raio_minimo);

        time_vals(alg_idx)=stats.time;
        bubbles_vals(alg_idx)=stats.num_bubbles;
        length_vals(alg_idx)=stats.path_length;
        sm_vals(alg_idx)=sm;

        todos_rosarios{alg_idx}=rosary;

    end

end

fprintf('\n');
fprintf('=============================================\n');
fprintf('RESULTADOS\n');
fprintf('=============================================\n');

fprintf('%-10s | %10s | %10s | %10s | %10s\n',...
'Algoritmo','Tempo','Bolhas','Path','SM');

for i=1:length(algorithms)

    if success_vals(i)

        fprintf('%-10s | %10.3f | %10d | %10.2f | %10.3f\n',...
        algorithms{i},...
        time_vals(i),...
        bubbles_vals(i),...
        length_vals(i),...
        sm_vals(i));

    else

        fprintf('%-10s | FALHA\n',algorithms{i});

    end
end

figure;
hold on;
axis equal;
view(3);

trisurf(...
faces,...
vertices(:,1),...
vertices(:,2),...
vertices(:,3),...
'FaceAlpha',0.4,...
'EdgeColor','none');

cores = {'r','g','b','m'};

for a=1:length(algorithms)

    rosary = todos_rosarios{a};

    if isempty(rosary)
        continue;
    end

    centros = zeros(length(rosary),3);

    for i=1:length(rosary)

        centros(i,:)=rosary{i}.center;

    end

    plot3(...
    centros(:,1),...
    centros(:,2),...
    centros(:,3),...
    cores{a},...
    'LineWidth',3);

end


objetoAPI_remota.simxFinish(id);
objetoAPI_remota.delete();