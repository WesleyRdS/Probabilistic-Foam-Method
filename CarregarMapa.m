function mapa = CarregarMapa(latitude, longitude, zoom)
% latitude => [minima maxima]
% longitude => [minuma maxima]
% zoom de cada tile

numero_tiles = 2^zoom; % Em um mapa do tipo slippy map cada nivel de zoom divide o mapa..
% em uma potência de base 2 elevado ao zoom tiles horizontal e verticalmente

% Projeção de Web Mercator(Transforma superficie esferica em um plano)
% Projeta latitude somando a tangente e a secante da latitude em radianos e
% depois normalizamos para ficar entre 0 e 1
latitude_em_tile = @(lat) (1 - log(tan(deg2rad(lat)) + sec(deg2rad(lat))) / pi)/ 2 * numero_tiles;

% Longitude varia de -180 a 180 totalizando 360 e normalizado para 1 depois
% mutiplica-se pelo numero de tiles na horizontal
longitude_em_tile = @(lon) (longitude + 180) / 360 * numero_tiles;

% Convertendo corenadas geograficas em indices de tile de X e Y
x_tiles = longitude_em_tile(longitude);
y_tiles = latitude_em_tile(latitude);

% Determina os indices inteiros dos tiles que precisa ser baixado
% floor(min(...)) => Primeiro tile que cobre a area
% ceil(max()) => Ultimo tile que cobre a area
x_minimo = floor(min(x_tiles));
x_maximo = ceil(max(x_tiles));
y_minimo = floor(min(y_tiles));
y_maximo = ceil(max(y_tiles));

tamanho_do_tile = 256; % Tamanho padrão de cada tile em pixels
largura = (x_maximo - x_minimo + 1) * tamanho_do_tile;
altura = (y_maximo - y_minimo + 1) * tamanho_do_tile;

%Cria uma imagem RGBvazia com o tamanho calculado
area_imagem = zeros(altura, largura, 3, 'uint8');

for x_i = x_minimo:x_maximo
    for y_i = y_minimo:y_minimo
        % Importa o tile do mapa sem textos
        url = sprintf('https://basemaps.cartocdn.com/light_nolabels/%d/%d/%d.png', zoom, x_i, y_i);
        tile = imread(url);
        %Melhora o contraste
        tile = imadjust(tile);

        % Se o tile for grayscale converte para RGB duplicando os dados
        if ndims(tile) == 2
            tile = repmat(tile, [1 1 3]);
        end

        deslocamento_de_x = (x_i - x_minimo) * tamanho_do_tile + 1;
        deslocamento_de_y = (y_i - y_minimo) * tamanho_do_tile + 1;
        
        % Calcula as posições na imagem final onde cada tile vai ser
        % colocado
        indice_de_y = deslocamento_de_y:deslocamento_de_y + tamanho_do_tile - 1;
        indice_de_x = deslocamento_de_x:deslocamento_de_x + tamanho_do_tile - 1;

        % Coloca o tile na posição correta da imagem
        mapa(indice_de_y, indice_de_x, :) = tile;
    end
end

% Inverte a imagem verticalmente porque no sistema de tiles o indice Y
% aumenta para baixo
mapa = flipud(mapa);

end




