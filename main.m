function [mapa,xmin,ymin,zmin,res] = ...
    create_stl_voxel_map(stl_file,res)

%% ===================== LEITURA STL =====================

TR = stlread(stl_file);

V = TR.Points;
F = TR.ConnectivityList;

%% ===================== LIMITES =====================

xmin = min(V(:,1));
xmax = max(V(:,1));

ymin = min(V(:,2))+2;
ymax = max(V(:,2));

zmin = min(V(:,3));
zmax = max(V(:,3));

Nx = ceil((xmax-xmin)/res)+1;
Ny = ceil((ymax-ymin)/res)+1;
Nz = ceil((zmax-zmin)/res)+1;

fprintf('Grade: %d x %d x %d\n',Nx,Ny,Nz);

mapa = false(Ny,Nx,Nz);

%% ===================== VOXELIZAÇÃO DA SUPERFÍCIE =====================

for f = 1:size(F,1)

    P1 = V(F(f,1),:);
    P2 = V(F(f,2),:);
    P3 = V(F(f,3),:);

    L12 = norm(P2-P1);
    L13 = norm(P3-P1);
    L23 = norm(P3-P2);

    maior_lado = max([L12 L13 L23]);

    N = max(3,ceil(maior_lado/res)*2);

    for i = 0:N

        for j = 0:(N-i)

            u = i/N;
            v = j/N;
            w = 1-u-v;

            P = u*P1 + v*P2 + w*P3;

            ix = round((P(1)-xmin)/res)+1;
            iy = round((P(2)-ymin)/res)+1;
            iz = round((P(3)-zmin)/res)+1;

            if ix>=1 && ix<=Nx && ...
               iy>=1 && iy<=Ny && ...
               iz>=1 && iz<=Nz

                mapa(iy,ix,iz) = true;

            end

        end

    end

end

fprintf('Superficie voxelizada.\n');

%% ===================== FECHAMENTO =====================

SE = ones(3,3,3);

mapa = imclose(mapa,SE);

%% ===================== PREENCHIMENTO VERTICAL =====================

for ix = 1:Nx

    for iy = 1:Ny

        col = find(mapa(iy,ix,:));

        if numel(col) >= 2

            zmin_col = min(col);
            zmax_col = max(col);

            mapa(iy,ix,zmin_col:zmax_col) = true;

        end

    end

end

fprintf('Volume preenchido.\n');

%% ===================== LIMPEZA =====================

mapa = imfill(mapa,'holes');

mapa = bwareaopen(mapa,10);

fprintf('Mapa limpo.\n');

end

function D = compute_distance_map(map,res)

D = bwdist(map);

D = D*res;

end

function [ix,iy,iz] = ...
    world_to_voxel(P,xmin,ymin,zmin,res)

ix = round((P(1)-xmin)/res)+1;
iy = round((P(2)-ymin)/res)+1;
iz = round((P(3)-zmin)/res)+1;

end

function radius = expand_bubble_3d(...
    center,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min)

[ix,iy,iz] = ...
world_to_voxel(...
    center,...
    xmin,ymin,zmin,...
    res);

[Ny,Nx,Nz] = size(D);

if ix<1 || ix>Nx || ...
   iy<1 || iy>Ny || ...
   iz<1 || iz>Nz

    radius = 0;
    return

end

radius = D(iy,ix,iz);

if radius < r_min
    radius = 0;
end

end

function [rosary, success, stats] = PFM_3D( ...
    q_init,...
    q_goal,...
    map,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min,...
    K)

tic;

F = {};
Q = {};

%% Bolha inicial

r_init = expand_bubble_3d( ...
    q_init,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min);

if r_init <= 0

    rosary = {};
    success = false;

    stats.time = toc;
    stats.num_bubbles = 0;
    stats.path_length = inf;
    stats.success = false;
    stats.max_iter = 0;

    return

end

F{1} = struct( ...
    'center',q_init,...
    'radius',r_init,...
    'parent',0);

Q{1} = 1;

%% Busca

success = false;
goal_bubble_idx = [];

max_iter = 3000;

for iter = 1:max_iter

    if isempty(Q)
        break;
    end

    parent_idx = Q{1};
    Q(1) = [];

    parent = F{parent_idx};

    %% Eq. 2

    N = K*floor(parent.radius/r_min);
    N = max(1,min(N,30));

    for n = 1:N

        %% direção aleatória

        direction = randn(1,3);

        if norm(direction) < 1e-6
            continue;
        end

        direction = direction/norm(direction);

        %% ponto na superfície

        q_surface = ...
            parent.center + ...
            direction*parent.radius;

        %% converter para voxel

        [ix,iy,iz] = world_to_voxel( ...
            q_surface,...
            xmin,ymin,zmin,...
            res);

        %% limites

        if ix < 1 || ix > size(map,2) || ...
           iy < 1 || iy > size(map,1) || ...
           iz < 1 || iz > size(map,3)

            continue;

        end

        %% obstáculo

        if map(iy,ix,iz)

            continue;

        end

        %% dentro de outra bolha?

        inside = false;

        for i = 1:length(F)

            if norm(q_surface - F{i}.center) < ...
                    F{i}.radius - 0.1

                inside = true;
                break;

            end

        end

        if inside
            continue;
        end

        %% raio via mapa de distância

        r_new = expand_bubble_3d( ...
            q_surface,...
            D,...
            xmin,ymin,zmin,...
            res,...
            r_min);

        if r_new < r_min
            continue;
        end

        %% adiciona bolha

        new_idx = length(F)+1;

        F{new_idx} = struct( ...
            'center',q_surface,...
            'radius',r_new,...
            'parent',parent_idx);

        Q{end+1} = new_idx;

        %% goal

        if norm(q_surface-q_goal) <= r_new

            success = true;
            goal_bubble_idx = new_idx;
            break;

        end

    end

    if success
        break;
    end

end

%% Resultado

if success

    rosary = ...
        extract_rosary_3d( ...
        F,...
        goal_bubble_idx);

    stats.path_length = ...
        compute_path_length_3d(rosary);

else

    rosary = {};
    stats.path_length = inf;

end

stats.time = toc;
stats.num_bubbles = length(F);
stats.success = success;
stats.max_iter = iter;

end

function [rosary, success, stats] = GBPF_3D( ...
    q_init,...
    q_goal,...
    map,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min,...
    bias)

tic;

F = {};

%% Bolha inicial

r_init = expand_bubble_3d( ...
    q_init,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min);

if r_init <= 0

    rosary = {};
    success = false;

    stats.time = toc;
    stats.num_bubbles = 0;
    stats.path_length = inf;
    stats.success = false;
    stats.max_iter = 0;

    return

end

F{1} = struct( ...
    'center',q_init,...
    'radius',r_init,...
    'parent',0);

%% Limites do STL

Nx = size(map,2);
Ny = size(map,1);
Nz = size(map,3);

xmax = xmin + (Nx-1)*res;
ymax = ymin + (Ny-1)*res;
zmax = zmin + (Nz-1)*res;

%% Busca

max_iter = 3000;

success = false;
goal_bubble_idx = [];

for iter = 1:max_iter

    %% Goal bias

    if rand < bias

        q_aux = q_goal;

    else

        q_aux = [ ...
            xmin + rand*(xmax-xmin), ...
            ymin + rand*(ymax-ymin), ...
            zmin + rand*(zmax-zmin)];

    end

    %% nearest bubble

    nearest_idx = 1;

    nearest_dist = ...
        norm(F{1}.center - q_aux);

    for i = 2:length(F)

        d = norm(F{i}.center - q_aux);

        if d < nearest_dist

            nearest_dist = d;
            nearest_idx = i;

        end

    end

    parent = F{nearest_idx};

    %% nearest config

    direction = q_aux - parent.center;

    if norm(direction) < 1e-6

        direction = randn(1,3);

    end

    direction = ...
        direction / norm(direction);

    q_near = ...
        parent.center + ...
        direction*parent.radius;

    %% converte para voxel

    [ix,iy,iz] = ...
        world_to_voxel( ...
        q_near,...
        xmin,ymin,zmin,...
        res);

    %% limites

    if ix < 1 || ix > size(map,2) || ...
       iy < 1 || iy > size(map,1) || ...
       iz < 1 || iz > size(map,3)

        continue;

    end

    %% obstáculo

    if map(iy,ix,iz)

        continue;

    end

    %% interior de outra bolha

    inside = false;

    for i = 1:length(F)

        if norm(q_near - F{i}.center) < ...
                F{i}.radius - 0.1

            inside = true;
            break;

        end

    end

    if inside
        continue;
    end

    %% novo raio

    r_new = expand_bubble_3d( ...
        q_near,...
        D,...
        xmin,ymin,zmin,...
        res,...
        r_min);

    if r_new < r_min
        continue;
    end

    %% adiciona bolha

    new_idx = length(F)+1;

    F{new_idx} = struct( ...
        'center',q_near,...
        'radius',r_new,...
        'parent',nearest_idx);

    %% goal

    if norm(q_near-q_goal) <= r_new

        success = true;
        goal_bubble_idx = new_idx;
        break;

    end

end

%% Resultado

if success

    rosary = ...
        extract_rosary_3d( ...
        F,...
        goal_bubble_idx);

    stats.path_length = ...
        compute_path_length_3d(rosary);

else

    rosary = {};
    stats.path_length = inf;

end

stats.time = toc;
stats.num_bubbles = length(F);
stats.success = success;
stats.max_iter = iter;

end

function [rosary, success, stats] = RBPF_3D( ...
    q_init,...
    q_goal,...
    map,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min,...
    K)

tic;

F = {};
OpenList = [];

%% Bolha inicial

r_init = expand_bubble_3d( ...
    q_init,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min);

if r_init <= 0

    rosary = {};
    success = false;

    stats.time = toc;
    stats.num_bubbles = 0;
    stats.path_length = inf;
    stats.success = false;
    stats.max_iter = 0;

    return

end

F{1} = struct( ...
    'center',q_init,...
    'radius',r_init,...
    'parent',0);

OpenList = 1;

%% Busca

max_iter = 3000;

success = false;
goal_bubble_idx = [];

for iter = 1:max_iter

    if isempty(OpenList)
        break;
    end

    %% Roleta baseada no raio

    radii = zeros(length(OpenList),1);

    for i = 1:length(OpenList)

        radii(i) = ...
            F{OpenList(i)}.radius;

    end

    if sum(radii) <= 0
        break;
    end

    prob = radii ./ sum(radii);

    r = rand;

    chosen_idx = ...
        find(r <= cumsum(prob), ...
        1,'first');

    parent_idx = ...
        OpenList(chosen_idx);

    parent = F{parent_idx};

    %% Remove da lista aberta

    OpenList(chosen_idx) = [];

    %% Número de filhos

    N = K * floor(parent.radius/r_min);

    N = max(1,min(N,30));

    for n = 1:N

        %% direção aleatória

        direction = randn(1,3);

        if norm(direction) < 1e-6
            continue;
        end

        direction = ...
            direction / norm(direction);

        %% superfície

        q_surface = ...
            parent.center + ...
            direction*parent.radius;

        %% voxel

        [ix,iy,iz] = ...
            world_to_voxel( ...
            q_surface,...
            xmin,ymin,zmin,...
            res);

        %% limites

        if ix < 1 || ix > size(map,2) || ...
           iy < 1 || iy > size(map,1) || ...
           iz < 1 || iz > size(map,3)

            continue;

        end

        %% obstáculo

        if map(iy,ix,iz)

            continue;

        end

        %% já coberto?

        covered = false;

        for i = 1:length(F)

            if norm(q_surface - F{i}.center) < ...
                    F{i}.radius - 0.1

                covered = true;
                break;

            end

        end

        if covered
            continue;
        end

        %% raio

        r_new = expand_bubble_3d( ...
            q_surface,...
            D,...
            xmin,ymin,zmin,...
            res,...
            r_min);

        if r_new < r_min
            continue;
        end

        %% adiciona

        new_idx = length(F)+1;

        F{new_idx} = struct( ...
            'center',q_surface,...
            'radius',r_new,...
            'parent',parent_idx);

        OpenList(end+1) = new_idx;

        %% goal

        if norm(q_surface-q_goal) <= r_new

            success = true;
            goal_bubble_idx = new_idx;
            break;

        end

    end

    if success
        break;
    end

end

%% Resultado

if success

    rosary = ...
        extract_rosary_3d( ...
        F,...
        goal_bubble_idx);

    stats.path_length = ...
        compute_path_length_3d(rosary);

else

    rosary = {};
    stats.path_length = inf;

end

stats.time = toc;
stats.num_bubbles = length(F);
stats.success = success;
stats.max_iter = iter;

end

function [rosary, success, stats] = HPF_3D( ...
    q_init,...
    q_goal,...
    map,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min,...
    K)

tic;

F = {};

OpenList = struct( ...
    'idx',{},...
    'g',{},...
    'h',{},...
    'f',{});

%% Bolha inicial

r_init = expand_bubble_3d( ...
    q_init,...
    D,...
    xmin,ymin,zmin,...
    res,...
    r_min);

if r_init <= 0

    rosary = {};
    success = false;

    stats.time = toc;
    stats.num_bubbles = 0;
    stats.path_length = inf;
    stats.success = false;
    stats.max_iter = 0;

    return

end

F{1} = struct( ...
    'center',q_init,...
    'radius',r_init,...
    'parent',0);

%% Custos iniciais

g_init = r_init;

h_init = norm(q_init - q_goal);

f_init = g_init + h_init;

OpenList(1).idx = 1;
OpenList(1).g   = g_init;
OpenList(1).h   = h_init;
OpenList(1).f   = f_init;

%% Busca

max_iter = 3000;

success = false;
goal_bubble_idx = [];

for iter = 1:max_iter

    if isempty(OpenList)
        break;
    end

    %% Menor f

    f_values = [OpenList.f];

    [~,min_idx] = min(f_values);

    current = OpenList(min_idx);

    OpenList(min_idx) = [];

    parent = F{current.idx};

    %% Número de filhos

    N = K * floor(parent.radius/r_min);

    N = max(1,min(N,30));

    for n = 1:N

        %% Direção aleatória

        direction = randn(1,3);

        if norm(direction) < 1e-6
            continue;
        end

        direction = direction/norm(direction);

        %% Superfície da bolha

        q_surface = ...
            parent.center + ...
            direction*parent.radius;

        %% Coordenadas voxel

        [ix,iy,iz] = ...
            world_to_voxel( ...
            q_surface,...
            xmin,ymin,zmin,...
            res);

        %% Limites

        if ix < 1 || ix > size(map,2) || ...
           iy < 1 || iy > size(map,1) || ...
           iz < 1 || iz > size(map,3)

            continue;

        end

        %% Colisão

        if map(iy,ix,iz)

            continue;

        end

        %% Cobertura

        covered = false;

        for i = 1:length(F)

            if norm(q_surface - F{i}.center) < ...
                    F{i}.radius - 0.1

                covered = true;
                break;

            end

        end

        if covered
            continue;
        end

        %% Novo raio

        r_new = expand_bubble_3d( ...
            q_surface,...
            D,...
            xmin,ymin,zmin,...
            res,...
            r_min);

        if r_new < r_min
            continue;
        end

        %% Nova bolha

        new_idx = length(F)+1;

        F{new_idx} = struct( ...
            'center',q_surface,...
            'radius',r_new,...
            'parent',current.idx);

        %% Equação 7

        g_new = current.g + r_new;

        h_new = norm(q_surface - q_goal);

        f_new = g_new + h_new;

        OpenList(end+1).idx = new_idx;
        OpenList(end).g = g_new;
        OpenList(end).h = h_new;
        OpenList(end).f = f_new;

        %% Goal

        if norm(q_surface-q_goal) <= r_new

            success = true;
            goal_bubble_idx = new_idx;

            break;

        end

    end

    if success
        break;
    end

end

%% Resultado

if success

    rosary = ...
        extract_rosary_3d( ...
        F,...
        goal_bubble_idx);

    stats.path_length = ...
        compute_path_length_3d(rosary);

else

    rosary = {};
    stats.path_length = inf;

end

stats.time = toc;
stats.num_bubbles = length(F);
stats.success = success;
stats.max_iter = iter;

end

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

function visualize_all_algorithms( ...
    map,...
    q_init,...
    q_goal,...
    todos_rosarios,...
    xmin,ymin,zmin,...
    res,...
    algorithms)

figure( ...
    'Name','Comparacao dos Algoritmos', ...
    'Position',[100 100 1400 900]);

hold on;

%% Obstáculos

[y_obs,x_obs,z_obs] = ...
    ind2sub(size(map),find(map));

X = xmin + (x_obs-1)*res;
Y = ymin + (y_obs-1)*res;
Z = zmin + (z_obs-1)*res;

scatter3( ...
    X,...
    Y,...
    Z,...
    8,...
    Z,...
    'filled');

colormap(jet);
colorbar;

%% Cores dos algoritmos

cores = { ...
    'r',...
    'g',...
    'b',...
    'm'};

%% Caminhos

for a = 1:length(todos_rosarios)

    rosary = todos_rosarios{a};

    if isempty(rosary)
        continue;
    end

    centers = zeros(length(rosary),3);

    for i = 1:length(rosary)

        centers(i,:) = rosary{i}.center;

    end

    plot3( ...
        centers(:,1),...
        centers(:,2),...
        centers(:,3),...
        cores{a},...
        'LineWidth',3);

    plot3( ...
        centers(:,1),...
        centers(:,2),...
        centers(:,3),...
        'o',...
        'Color',cores{a},...
        'MarkerFaceColor',cores{a},...
        'MarkerSize',5);

end

%% Start

plot3( ...
    q_init(1),...
    q_init(2),...
    q_init(3),...
    'ks',...
    'MarkerFaceColor','g',...
    'MarkerSize',12);

%% Goal

plot3( ...
    q_goal(1),...
    q_goal(2),...
    q_goal(3),...
    'kd',...
    'MarkerFaceColor','b',...
    'MarkerSize',12);

xlabel('X');
ylabel('Y');
zlabel('Z');

title('PFM vs GBPF vs RBPF vs HPF');

grid on;
axis equal;

view(45,30);

camlight;
lighting gouraud;

legend({ ...
    'Obstacles',...
    algorithms{1},...
    algorithms{2},...
    algorithms{3},...
    algorithms{4},...
    'Start',...
    'Goal'},...
    'Location','best');

end

%MAIN
clear
clc
close all

%% STL

res = 1;

[map,xmin,ymin,zmin,res] = ...
    create_stl_voxel_map("UEFScompleto.stl",res);

D = compute_distance_map(map,res);

%% Configuração
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

% Interface 2D para selecionar destino do drone

% Interface 2D para selecionar destino do drone
TR = stlread('UEFScompleto.stl');

vertices = TR.Points;
faces = TR.ConnectivityList;
figure;
trisurf(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
view(2)
axis equal
title("Selecione o local de destino")
hold on
ponto = (PosInicialDrone) / escala;
q_init = ponto; %Definindo posição inicial do drone como ponto de partida
%Plotando na representação superior 2D do mapa a localização do drone
plot(ponto(1), ponto(2), 'go', 'LineWidth', 2)

while true
    [x, y] = ginput(1); %Pega posição do mouse na figura ao clicar
    ponto = [x, y, ponto(3)];
    
    q_goal = ponto;
    plot(x, y, 'ro', 'LineWidth', 2)
    break
   
end
hold off


r_min = 0.5;
K = 30;
bias = 0.25;

%% Algoritmos

algorithms = { ...
    'PFM',...
    'GBPF',...
    'RBPF',...
    'HPF'};

todos_rosarios = cell(length(algorithms),1);
todos_stats = cell(length(algorithms),1);

for a = 1:length(algorithms)

    fprintf('\n====================\n');
    fprintf('%s\n',algorithms{a});
    fprintf('====================\n');

    switch algorithms{a}

        case 'PFM'

            [rosary,success,stats] = ...
                PFM_3D( ...
                q_init,...
                q_goal,...
                map,...
                D,...
                xmin,ymin,zmin,...
                res,...
                r_min,...
                K);

        case 'GBPF'

            [rosary,success,stats] = ...
                GBPF_3D( ...
                q_init,...
                q_goal,...
                map,...
                D,...
                xmin,ymin,zmin,...
                res,...
                r_min,...
                bias);

        case 'RBPF'

            [rosary,success,stats] = ...
                RBPF_3D( ...
                q_init,...
                q_goal,...
                map,...
                D,...
                xmin,ymin,zmin,...
                res,...
                r_min,...
                K);

        case 'HPF'

            [rosary,success,stats] = ...
                HPF_3D( ...
                q_init,...
                q_goal,...
                map,...
                D,...
                xmin,ymin,zmin,...
                res,...
                r_min,...
                K);

    end

    todos_rosarios{a} = rosary;
    todos_stats{a} = stats;

    fprintf('Success: %d\n',stats.success);
    fprintf('Bubbles: %d\n',stats.num_bubbles);
    fprintf('Time: %.3f s\n',stats.time);

end

visualize_all_algorithms( ...
    map,...
    q_init,...
    q_goal,...
    todos_rosarios,...
    xmin,ymin,zmin,...
    res,...
    algorithms);