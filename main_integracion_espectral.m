clc; clear; close all;

%% --- CONFIGURACIÓN GENERAL ---
carpeta = 'ArchivosIQ';
if ~exist(carpeta, 'dir'), error(['La carpeta "' carpeta '" no existe.']); end
archivos = dir(fullfile(carpeta, '*.iq'));
if isempty(archivos), error('No se encontraron archivos .iq.'); end

disp('--- Integración Espectral (Todos los Métodos) ---');
for i = 1:length(archivos)
    fprintf('%d: %s (%.2f MB)\n', i, archivos(i).name, archivos(i).bytes/1024/1024);
end
fprintf('\n');

idx = input('Selecciona el número del archivo: ');
if isempty(idx) || idx < 1 || idx > length(archivos), error('Selección inválida.'); end
archivoSeleccionado = archivos(idx);

disp(' ');
inicio_pct = input('Inicio (%): '); if isempty(inicio_pct), inicio_pct=0; end
fin_pct = input('Fin (%): '); if isempty(fin_pct), fin_pct=100; end
if fin_pct - inicio_pct > 100, fin_pct = inicio_pct + 100; end

%% --- EJECUCIÓN EN LOTE ---
run_batch_mode = true;
batch_idx = idx;
batch_inicio_pct = inicio_pct;
batch_fin_pct = fin_pct;

results_batch = struct();

% 1. Método Directo
fprintf('\n>>> Ejecutando Método Directo...\n');
try
    main_directo;
catch ME
    warning('Error en Método Directo: %s', ME.message);
end

% 2. Método Indirecto
fprintf('\n>>> Ejecutando Método Indirecto...\n');
try
    main_indirecto;
catch ME
    warning('Error en Método Indirecto: %s', ME.message);
end

% 3. MUSIC
fprintf('\n>>> Ejecutando MUSIC...\n');
try
    main_music;
catch ME
    warning('Error en MUSIC: %s', ME.message);
end

% 4. AR (Yule-Walker)
fprintf('\n>>> Ejecutando AR...\n');
try
    main_ar;
catch ME
    warning('Error en AR: %s', ME.message);
end

% 5. CWT
fprintf('\n>>> Ejecutando CWT...\n');
try
    main_cwt;
catch ME
    warning('Error en CWT: %s', ME.message);
end

% 6. Superlet
fprintf('\n>>> Ejecutando Superlet...\n');
try
    main_Superlet;
catch ME
    warning('Error en Superlet: %s', ME.message);
end

%% --- VISUALIZACIÓN INTEGRADA ---
fprintf('\nGenerando Gráfico Integrado...\n');

figInt = figure('Name', 'Integración Espectral', 'Color', 'w', 'Units', 'normalized', 'Position', [0.3 0.05 0.4 0.9]);
tiledlayout(6, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

metodos = {'Directo', 'Indirecto', 'MUSIC', 'AR', 'CWT', 'Superlet'};
fields  = {'Directo', 'Indirecto', 'MUSIC', 'AR', 'CWT', 'Superlet'};
% Colores distintivos (RGB)
colores_list = [
    0 0.4470 0.7410;  % Azul
    0.8500 0.3250 0.0980; % Rojo anaranjado
    0.9290 0.6940 0.1250; % Amarillo ocre
    0.4940 0.1840 0.5560; % Violeta
    0.4660 0.6740 0.1880; % Verde
    0.6350 0.0780 0.1840  % Rojo vino
    ];

f_HI = 1420e6;

for k = 1:length(metodos)
    nexttile;
    metodo = metodos{k};
    field = fields{k};
    color_linea = colores_list(k, :);

    if isfield(results_batch, field)
        res = results_batch.(field);
        f_axis = res.f;
        P_axis = res.P;

        plot(f_axis/1e6, P_axis, 'Color', color_linea, 'LineWidth', 1.5);
        grid on;
        xline(f_HI/1e6, 'k--', 'LineWidth', 1); % Línea HI en negro para contraste

        % Ajuste de límites inteligente (Asegurar que se vea la línea HI)
        min_f = min(f_axis);
        max_f = max(f_axis);
        % Expandir un poco si HI está en el borde o fuera
        if f_HI < min_f || f_HI > max_f
            % Si HI está fuera, lo incluimos forzosamente
            min_f = min(min_f, f_HI - 10e3); % -10kHz margen
            max_f = max(max_f, f_HI + 10e3); % +10kHz margen
        end
        xlim([min_f max_f]/1e6);

        % Título y Etiquetas
        title(metodo, 'FontWeight', 'bold');
        if k == length(metodos)
            xlabel('Frecuencia (MHz)');
        else
            xticklabels({});
        end
        ylabel('Potencia (dBm)');
    else
        text(0.5, 0.5, ['Datos de ' metodo ' no disponibles'], ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        axis off;
    end
end

sgtitle(['Integración Espectral | ' archivoSeleccionado.name ' | Inicio: ' num2str(inicio_pct) '% - Fin: ' num2str(fin_pct) '%'], 'Interpreter', 'none');

% Guardar Figura
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
[~, name_no_ext, ~] = fileparts(archivoSeleccionado.name);
dir_salida = fullfile('Resultados_Imagenes', 'INTEGRACION');
if ~exist(dir_salida, 'dir'), mkdir(dir_salida); end

saveas(figInt, fullfile(dir_salida, ['Integracion_' name_no_ext '_' timestamp '.png']));
saveas(figInt, fullfile(dir_salida, ['Integracion_' name_no_ext '_' timestamp '.fig']));

disp('--- Integración Completada ---');
