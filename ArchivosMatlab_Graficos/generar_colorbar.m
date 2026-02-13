clc; clear; close all;

%% Configuración
rangoColores = [-90 -40]; % Rango consistente con los otros scripts
cmap = jet; % Mapa de colores usado anteriormente

%% Generar Figura solo con Colorbar
fig = figure('Name', 'Colorbar Horizontal', 'Color', 'w', ...
    'Units', 'normalized', 'Position', [0.3 0.4 0.4 0.15]); % Figura pequeña y ancha

% Crear ejes invisibles para alojar el colorbar
ax = axes('Position', [0.1 0.3 0.8 0.0], 'Visible', 'off');
clim(rangoColores);
colormap(cmap);

% Crear colorbar horizontal
c = colorbar('Location', 'north', 'Position', [0.1 0.4 0.8 0.3]);
c.Label.String = 'Potencia (dBm)';
c.Label.FontSize = 12;
c.FontSize = 10;
c.LineWidth = 1;
c.TickDirection = 'out';

% Guardar imagen
if ~exist('Resultados_Imagenes', 'dir')
    mkdir('Resultados_Imagenes');
end
saveas(fig, fullfile('Resultados_Imagenes', 'Colorbar_Horizontal.png'));
saveas(fig, fullfile('Resultados_Imagenes', 'Colorbar_Horizontal.fig'));

disp('Colorbar horizontal generado en Resultados_Imagenes/Colorbar_Horizontal.png');
