clc; clear; close all; % limpa tudo (console, variáveis e figuras)
%--------------------Carregar mapa 2D e processar imagem-------------------

latitudes = [-12.1980927, -12.1988903];
longitudes = [-38.9725052, -38.970194];
zoom = 17;
area = CarregarMapa(latitudes, longitudes, zoom);
area = im2double(area); % Transforma a imagem no range de 0 a 1
[R,G,B] = imsplit(area); % Separa os canais da imagem

% Intensidade dos canais RGB dos obstaculos do mapa transformado no range de
% 0 a 1
%cores_obstaculos = [[252 248 239]/255 ; [231 232 235]/255];
cores_obstaculos = [[107 107 107]/255 ; [170 170 170]/255 ; [151 151 151]/255 ; [143 143 143]/255 ...
     ; [146 146 146]/255 ; [134 134 134]/255 ; [163 163 163]/255 ; [121 121 121]/255 ; [99 99 99]/255 ; [128 128 128]/255];

% Tolerância para binarização da imagem
tol = 0.02;

img_binarizada = false(size(R));

% Binarização da imagem de acordo com  a tolerancia(Obstaculos serão transformados
% em 0 e areas livres em 1)
for i = 1:size(cores_obstaculos,1)
    img_binarizada = img_binarizada | ...
        ( abs(R - cores_obstaculos(i,1)) < tol & ...
          abs(G - cores_obstaculos(i,2)) < tol & ...
          abs(B - cores_obstaculos(i,3)) < tol );
end


% Inversão dos valores
livre = ~img_binarizada;

% Elemento estrutural no formato de esfera com raio 1. Usado para
% fechamento da imagem(Vai realçar o que for branco e tirar ruidos da imagem)
elemento_estrutural = strel("sphere", 1);
mapa_filtrado = imclose(livre, elemento_estrutural);

mapa_filtrado([1 end],:) = 0;
mapa_filtrado(:,[1 end]) = 0;

% Obstaculos
mapa2D = ~mapa_filtrado;

%----------------------------Geração de mapa 3D---------------------------
% Cria relevo continuo ao multiplicar a matriz por 5 e aplicar suavização
%gaussiana
altura_base = imgaussfilt(double(mapa2D)*10, 3);
altura_base = round(altura_base);

altura_base(altura_base < 1) = 0;

% Define altura maxima do mundo tomando como base a altura dos obstaculos e
% deixando mais 15 espaços livres acima
max_altura_mundo = max(altura_base(:)) + 15;

% Criação do mapa 3D
mapa3D = zeros(size(mapa2D,1), size(mapa2D,2), max_altura_mundo);

for i = 1:size(mapa2D,1)
    for j = 1:size(mapa2D,2)

        altura_obstaculos = altura_base(i,j);

        if altura_obstaculos > 0
            % Gerando obstaculos(1) no eixo z até a altura dos edificios
            mapa3D(i,j,1:altura_obstaculos) = 1;
        end

        % Acima disso permanece livre (0)
    end
end

%-----------------------Calcula do mapa de distância-----------------------

distancia3D = bwdist(mapa3D);

% Incentivo para o movimento para cima   
[yy,xx,zz] = ndgrid(1:size(mapa3D,1),1:size(mapa3D,2),1:size(mapa3D,3));
peso_z = 0.5; % Peso para favorecer a altura(Quanto maior mais incentivo)

% Adiciona um bonus no eixo z que diz para o algoritimo que a subida é mais
% vantajosa ou seja que para cima tem mais espaços livres caso contrario o
% algoritimo escolheria distancias maiores
distancia3D = distancia3D + peso_z * (zz / size(mapa3D,3));

distancia3D = distancia3D - 2; % Remove valores muito pequenos
distancia3D(distancia3D < 0) = 0;

figure;
imagesc(area);
colormap(gray)
axis equal
title("Clique na PARTIDA e depois no DESTINO")
hold on

%----------------------Selecionar os pontos de interesse-------------------

% PARTIDA
while true
    [x, y] = ginput(1);
    x = round(x); y = round(y);

    if mapa2D(y,x) == 0
        partida = [x y 2]; % Começa um pouco acima do chão
        plot(x,y,'go','LineWidth',2)
        break
    end
end

% DESTINO
while true
    [x, y] = ginput(1);
    x = round(x); y = round(y);

    if mapa2D(y,x) == 0
        destino = [x y 2];
        plot(x,y,'ro','LineWidth',2)
        break
    end
end

hold off

%--------------------Metodo de espuma probabilistica-----------------------

raio_minimo = 2; 
k = 10; % Numero de amostras(Exploração)
fator_z = 2.5; % Ganho na direção vertical(z) durante a geração de pontos

espuma = []; % Armazena todas as bolhas
fila = []; % Guarda indices das bolhas que vão gerar bolhas filhas

% Função inline que recebe uma posição x,y e retorna o valor de x e y no
% mapa de distância
expandir_bolha = @(x, y, z) distancia3D(y, x, z);
r_0 = expandir_bolha(partida(1), partida(2), partida(3));
espuma = [partida r_0 0]; % Adiciona a primeira bolha (Centro, raio, bolha pai)
fila = [1];

encontrou_destino = false;
indice_do_destino = - 1;

% Enquanto ouver bolhas para expandir
while ~isempty(fila)
    indice = fila(1); % Pega a primeira bolha da fila
    fila(1) = []; % Remove essa bolha da fila
    centro_da_bolha_pai = espuma(indice, 1:3); % Pega as colunas 1, 2 e 3 da linha da bolha
    raio_da_bolha_pai = espuma(indice, 4);

     % Define quantas novas bolhas serão geradas
    numero_de_novas_bolhas = ceil(k * (raio_da_bolha_pai/raio_minimo));

    for i = 1:numero_de_novas_bolhas
        theta = 2*pi*rand; % Gera um angulo aleatorio.
        phi = acos(2*rand - 1);
        
        %Deslocamento nas três dimenões
        dx = raio_da_bolha_pai*sin(phi)*cos(theta);
        dy = raio_da_bolha_pai*sin(phi)*sin(theta);
        dz = fator_z * raio_da_bolha_pai*cos(phi);
        
        % Nova posição candidata para bolhas filhas
        novo_ponto = round(centro_da_bolha_pai + [dx dy dz]);
        x = novo_ponto(1); y = novo_ponto(2); z = novo_ponto(3);
        
        % Ignorando regiões fora do mapa
        if x<1 || y<1 || z<1 || ...
           x>size(mapa3D,2) || ...
           y>size(mapa3D,1) || ...
           z>size(mapa3D,3)
            continue
        end

        % Ignorando obstaculos
        if mapa3D(y,x,z) == 1
            continue
        end

        % Se a lista de espumas não estiver vazia
        if ~isempty(espuma)
            % ponto ,Norma 2(euclidiana), 2-calcula a norma ao longo das linha
            distancias = vecnorm(espuma(:,1:3)-novo_ponto,2,2); 
            % Verificando se o centro esta localizado dentro de outra
            % bolha. Caso sim, ignora esse ponto. 
            if any(distancias < espuma(:,4))
                continue
            end
        end

        % Expandindo nova bolha
        raio_novo = expandir_bolha(x,y,z);
        
        % Se o raio maximo da bolha for menor do que o raio minimo ignora
        if raio_novo < raio_minimo
            continue
        end
        
        %Adicionando nova bolha(linha) a matriz de espuma
        espuma = [espuma; novo_ponto raio_novo indice];
        fila = [fila size(espuma,1)];
        
        % Verifica se o destino esta no raio da esfera
        if norm(novo_ponto - destino) <= raio_novo
            % cria bolha exatamente no destino
            raio_destino = expandir_bolha(destino(1), destino(2), destino(3));
            espuma = [espuma; destino raio_destino size(espuma,1)];
            encontrou_destino = true;
            indice_destino = size(espuma,1);
            break % Para caso esteja
        end
    end

    if encontrou_destino
        break
    end
end

%--------------Geração do caminho entre as bolhas--------------------------

caminho = [];

% Começa a gerar o caminho se foi encontrado
if encontrou_destino
    indice = indice_destino;

    while indice ~= 0
        % Monta o caminho
        caminho = [espuma(indice,1:3); caminho];
        indice = espuma(indice,5);
    end
    caminho = smoothdata(caminho, 1, 'gaussian', 5);
else
    disp("Destino não encontrado")
end


%-----------------------Vizualização do metodo-----------------------------
figure; hold on; axis equal;
view(3)

% ind2sub: converte índices lineares em coordenadas.(x->linha, y->coluna,
% z->altura) e find(mapa3D) pega todos os indices onde mapa3D = 1
[x,y,z] = ind2sub(size(mapa3D), find(mapa3D));
plot3(y,x,z,'ks','MarkerSize',3,'MarkerFaceColor','k')

for i = 1:size(espuma,1)
    [sx,sy,sz] = sphere(8);
    surf(sx*espuma(i,4)+espuma(i,1), ...
         sy*espuma(i,4)+espuma(i,2), ...
         sz*espuma(i,4)+espuma(i,3), ...
         'FaceAlpha',0.07,'EdgeColor','none');
end

if ~isempty(caminho)
    plot3(caminho(:,1), caminho(:,2), caminho(:,3), 'r','LineWidth',3)
end

plot3(partida(1),partida(2),partida(3),'go','LineWidth',2)
plot3(destino(1),destino(2),destino(3),'ro','LineWidth',2)

title('Espuma Probabilistica 3D')
xlabel('X'); ylabel('Y'); zlabel('Z');