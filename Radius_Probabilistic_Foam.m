clc; clear; close all; % limpa tudo (console, variáveis e figuras)

latitudes = [-12.1980927, -12.1988903];
longitudes = [-38.9725052, -38.970194];
zoom = 17;
area = CarregarMapa(latitudes, longitudes, zoom);
area = im2double(area); % Transforma a imagem no range de 0 a 1
imshow(area)
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

mapa_filtrado(1,: ) = 0;
mapa_filtrado(end, :) = 0;
mapa_filtrado(:, 1) = 0;
mapa_filtrado(:, end) = 0;

% A função bwdist calcula a distancia de cada pixel 0 até o proximo pixel 1 mais
% proximo. O mapa é filtrado para inverter os valores já que obstaculos são
% considerados 1 no nosso codigo
mapa_de_distancia = bwdist(~mapa_filtrado);

mapa_de_distancia(1,: ) = 0;
mapa_de_distancia(end, :) = 0;
mapa_de_distancia(:, 1) = 0;
mapa_de_distancia(:, end) = 0;

% Margem de segurança para bolha
margem_para_bolha = 3;

mapa_de_distancia = mapa_de_distancia - margem_para_bolha;
mapa_de_distancia(mapa_de_distancia < 0) = 0; % Evita numeros negativos.

figure
imshow(area)
title("Clique no mapa para selecionar o ponto de partida e o destino!!")
hold on

while true
    [x, y] = ginput(1); % ginput recebe clique do mouse
    x = round(x);
    y = round(y);

    if mapa_filtrado(y,x)
        partida = [x y];
        plot(x,y,"ro","LineWidth",2, "Color", [1 0 0])
        break
    else
        disp("Clique em uma área livre!!")
    end
end


while true
    [x, y] = ginput(1);
    x = round(x);
    y = round(y);

    if mapa_filtrado(y,x)
        destino = [x y];
        plot(x,y,"ro","LineWidth",2, "Color",[0 1 0])
        break
    else
        disp("Clique em uma área livre!!")
    end
end

hold off

% Probabilist Foam Method

raio_minimo = 2; 
k = 10; % Numero de amostras(Exploração)

espuma = []; % Armazena todas as bolhas
fila = []; % Guarda indices das bolhas que vão gerar bolhas filhas

% Função inline que recebe uma posição x,y e retorna o valor de x e y no
%mapa de distância
espandir_bolha = @(x, y) mapa_de_distancia(y, x);

r_0 = espandir_bolha(partida(1), partida(2));
espuma = [partida r_0 0]; % Adiciona a primeira bolha (Centro, raio, bolha pai)

fila = [1]; % Adiciona bolha a fila(Indice 1)

encontrou_destino = false;
indice_do_destino = - 1;

% Enquanto ouver bolhas para expandir
while ~isempty(fila)
    indice = fila(1); % Pega a primeira bolha da fila
    fila(1) = []; % Remove essa bolha da fila
    centro_da_bolha_pai = espuma(indice, 1:2); % Pega as colunas 1 e 2 da linha da bolha
    raio_da_bolha_pai = espuma(indice, 3);
    
    % Define quantas novas bolhas serão geradas
    numero_de_novas_bolhas = ceil(k * (raio_da_bolha_pai/raio_minimo));

    for i = 1:numero_de_novas_bolhas
        theta = 2*pi*rand; % Gera um angulo aleatorio.

        x_n_bolha = round(centro_da_bolha_pai(1) + raio_da_bolha_pai*cos(theta));
        y_n_bolha = round(centro_da_bolha_pai(2) + raio_da_bolha_pai*sin(theta));
        
        % Garante bolha dentro da area da imagem
        if x_n_bolha < 1 || y_n_bolha < 1 || ...
           x_n_bolha > size(mapa_filtrado, 2) || ...
           y_n_bolha > size(mapa_filtrado, 1)
            continue
        end
        
        % Garante bolha em espaço livre
        if ~mapa_filtrado(y_n_bolha, x_n_bolha)
            continue
        end

        ponto_valido = [x_n_bolha, y_n_bolha];

        if ~isempty(espuma) %Verifica se tem bolhas dentro da espuma
            distancias = vecnorm(espuma(:, 1:2) - ponto_valido, 2, 2); % ponto ,Norma 2(euclidiana), 2-calcula a norma ao longo das linhas
            if any(distancias < espuma(:, 3)) % Se o ponto esta dentro de qualquer bolha existente ignore
                continue
            end
        end

        raio_novo = espandir_bolha(x_n_bolha, y_n_bolha);
        
        % Descarta bolhas com raio menor que o minimo
        if raio_novo < raio_minimo
            continue
        end
        
        %Adicionando nova bolha(linha) a matriz de espuma
        espuma = [espuma; x_n_bolha y_n_bolha raio_novo indice];
        
        novo_indice = size(espuma, 1);
        fila = [fila novo_indice];

        if norm(centro_da_bolha_pai - destino) <= raio_novo
            disp("Destino localizado!!")
            encontrou_destino = true;

            espuma = [espuma; destino mapa_de_distancia(destino(2), destino(1)) novo_indice];
            indice_do_destino = size(espuma, 1);
            break
        end
    end

    if encontrou_destino
        break
    end
end

if ~encontrou_destino
    disp("Destino não foi localizado!!")
end

% Geração do caminho entre as bolhas

caminho = [];

indice = indice_do_destino;

if encontrou_destino
    while indice ~= 0
        caminho = [espuma(indice, 1:2); caminho]; %Adiciona ponto ao caminho
        indice = espuma(indice, 4); % Vai para o pai
    end
end

figure
imshow(area)
hold on

% Desenhar Bolhas
for i = 1:size(espuma,1)
    viscircles(espuma(i,1:2), espuma(i,3), ...
        'Color','b','LineWidth',0.3);
end

% Desenhar caminho
if ~isempty(caminho)
    plot(caminho(:,1), caminho(:,2), 'r','LineWidth',2)
end

plot(partida(1), partida(2),'go','LineWidth',2)
plot(destino(1), destino(2),'ro','LineWidth',2)

title("Probabilistic Foam Method")








