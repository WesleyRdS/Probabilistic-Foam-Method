clc;
clear;
close all;

%% ===================== LEITURA STL =====================

TR = stlread('modulo3.stl');

V = TR.Points;
F = TR.ConnectivityList;

%% ===================== RESOLUÇÃO =====================

res = 1.0;   % metros

xmin = min(V(:,1));
xmax = max(V(:,1));

ymin = min(V(:,2));
ymax = max(V(:,2));

zmin = min(V(:,3));
zmax = max(V(:,3));

Nx = ceil((xmax-xmin)/res)+1;
Ny = ceil((ymax-ymin)/res)+1;
Nz = ceil((zmax-zmin)/res)+1;

fprintf('Grade: %d x %d x %d\n',Nx,Ny,Nz);

mapa = false(Ny,Nx,Nz);

%% ===================== RASTERIZAÇÃO DOS TRIÂNGULOS =====================

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

                mapa(iy,ix,iz)=true;

            end

        end

    end

end

fprintf('Superfície voxelizada.\n');

%% ===================== FECHAMENTO MORFOLÓGICO =====================

SE = ones(3,3,3);
mapa = imclose(mapa,SE);

%% ===================== PREENCHIMENTO VERTICAL =====================

for ix = 1:Nx

    for iy = 1:Ny

        col = find(mapa(iy,ix,:));

        if numel(col)>=2

            zmin_col = min(col);
            zmax_col = max(col);

            mapa(iy,ix,zmin_col:zmax_col)=true;

        end

    end

end

fprintf('Volume preenchido.\n');

%% ===================== VISUALIZAÇÃO =====================

[y,x,z] = ind2sub(size(mapa),find(mapa));

X = xmin + (x-1)*res;
Y = ymin + (y-1)*res;
Z = zmin + (z-1)*res;

figure;

scatter3( ...
    X,...
    Y,...
    Z,...
    5,...
    Z,...
    'filled');

axis equal;
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');

title('Volume voxelizado');

view(3);

%% ===================== DISTANCE TRANSFORM =====================

livre = ~mapa;

D = bwdist(mapa)*res;

figure;

xs = round(Nx/2);
ys = round(Ny/2);
zs = round(Nz/2);

slice(double(D),xs,ys,zs)

shading interp
colorbar
title('Mapa de distância')

%% ===================== EXEMPLO DE CONSULTA =====================

p = [0 0 10];

ix = round((p(1)-xmin)/res)+1;
iy = round((p(2)-ymin)/res)+1;
iz = round((p(3)-zmin)/res)+1;

if ix>=1 && ix<=Nx && ...
   iy>=1 && iy<=Ny && ...
   iz>=1 && iz<=Nz

    raio = D(iy,ix,iz);

    fprintf('Raio livre = %.3f m\n',raio);

end