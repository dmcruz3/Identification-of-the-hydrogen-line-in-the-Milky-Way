clc;
clear;
close all;

%% 1. SELECCIÓN DE ARCHIVO
disp('--- Visualizador Universal de Matrices Espectrales (.mat) ---');
start_path = fullfile(pwd, 'Resultados_Datos');
if ~exist(start_path, 'dir'), start_path = pwd; end

[file, path] = uigetfile(fullfile(start_path, '*.mat'), 'Selecciona un archivo .mat');
if isequal(file, 0)
    disp('Selección cancelada.');
    return;
end
fullpath = fullfile(path, file);
disp(['Cargando: ' fullpath]);

%% 2. CARGA DE DATOS
data_struct = load(fullpath);
vars = fieldnames(data_struct);

% Variables objetivo
F = [];
T = [];
P = [];

%% 3. DETECCIÓN AUTOMÁTICA DE VARIABLES
disp('Analizando variables...');

% --- Buscar FRECUENCIA (F) ---
% Posibles nombres: f, f_abs, F_abs, freq, F, F_stft, frequency, Freq_Save
nombres_f = {'Freq_Save', 'f_abs', 'F_abs', 'f', 'F', 'freq', 'frequency', 'F_stft'};
for i = 1:length(nombres_f)
    if isfield(data_struct, nombres_f{i})
        temp = data_struct.(nombres_f{i});
        if isvector(temp) && length(temp) > 1
            F = temp;
            disp([' > Frecuencia encontrada: ' nombres_f{i}]);
            break;
        end
    end
end

% --- Buscar TIEMPO (T) ---
% Posibles nombres: t, t_abs, T_abs, time, T, t_signal, T_stft, Accum_Time
nombres_t = {'Accum_Time', 't_abs', 'T_abs', 't', 'T', 'time', 't_signal', 'T_stft', 'time_vector'};
for i = 1:length(nombres_t)
    if isfield(data_struct, nombres_t{i})
        temp = data_struct.(nombres_t{i});
        if isvector(temp) && length(temp) > 1
            T = temp;
            disp([' > Tiempo encontrado: ' nombres_t{i}]);
            break;
        end
    end
end

% --- Buscar POTENCIA / ESPECTROGAM (P) ---
% Posibles nombres: P_dBm, S, Spec, P, data, power, spectrogram, Accum_P
nombres_p = {'Accum_P', 'P_dBm', 'S', 'Spec', 'P', 'data', 'power', 'spectrogram'};
for i = 1:length(nombres_p)
    if isfield(data_struct, nombres_p{i})
        temp = data_struct.(nombres_p{i});
        if ~isvector(temp) && ismatrix(temp) && numel(temp) > 100
            % Verificar dimensiones sean consistentes con F y T si ya existen
            dim_ok = true;
            if ~isempty(F)
                if size(temp, 1) ~= length(F) && size(temp, 2) ~= length(F)
                    dim_ok = false; % No coincide con Frecuencia
                end
            end

            if dim_ok
                P = temp;

                % Si es complejo (como 'S'), sacar modulo cuadrado y log
                if ~isreal(P)
                    disp('   (Convirtiendo datos complejos a Potencia dB)');
                    P = 10*log10(abs(P).^2 + eps);
                end

                disp([' > Matriz de Datos encontrada: ' nombres_p{i}]);
                break;
            end
        end
    end
end

% --- Fallback: Heurística por tamaño ---
if isempty(P)
    disp(' ! Buscando matriz más grande como P...');
    max_el = 0;
    for i = 1:length(vars)
        val = data_struct.(vars{i});
        if isnumeric(val) && ~isvector(val) && ismatrix(val)
            if numel(val) > max_el
                max_el = numel(val);
                P = val;
                if ~isreal(P), P = 10*log10(abs(P).^2 + eps); end
                disp([' ! Asumiendo P es: ' vars{i}]);
            end
        end
    end
end

if isempty(P)
    error('No se pudo encontrar una matriz de datos válida en el archivo.');
end

% --- Reconstrucción de Ejes si faltan ---
if isempty(F)
    disp(' ! Frecuencia no encontrada. Generando eje índice.');
    F = 1:size(P, 1);
end
if isempty(T)
    disp(' ! Tiempo no encontrado. Generando eje índice.');
    T = 1:size(P, 2);
end

% --- Asegurar Orientación ---
% Convención: Filas = Frecuencia, Columnas = Tiempo
% F debe tener longitud size(P, 1) y T size(P, 2)
if length(F) == size(P, 2) && length(T) == size(P, 1)
    disp(' ! Transponiendo P para coincidir con ejes encontrados...');
    P = P.';
end

% Verificar nuevamente dimensiones
if length(F) ~= size(P, 1) || length(T) ~= size(P, 2)
    % Intento de swap de ejes
    if length(F) == size(P, 2) && length(T) == size(P, 1)
        P = P.';
        disp(' ! Transpuesto P (segundo intento).');
    else
        warning('Las dimensiones de F y T no coinciden exactamente con P. Se intentará graficar igual.');
        % Truncar o ajustar
        min_f = min(length(F), size(P, 1));
        min_t = min(length(T), size(P, 2));
        F = F(1:min_f);
        T = T(1:min_t);
        P = P(1:min_f, 1:min_t);
    end
end

%% 4. VISUALIZACIÓN
rangoColores = [-90 -40]; % Default
% Auto-ajuste de rango si los datos son muy distintos
p_med = median(P(:));
if abs(p_med - mean(rangoColores)) > 30
    % Si la media está muy lejos, ajustar
    rangoColores = [p_med-20, p_med+30];
    disp([' ! Ajustando rango de colores a: ' num2str(rangoColores)]);
end

disp('Generando gráficas...');

% --- 2D Espectrograma ---
fig2D = figure('Name', ['2D: ' file], 'Color', 'w');
ax2D = axes;
surf(ax2D, F/1e6, T, P.', 'EdgeColor', 'none');
view(0, 90);
axis tight;
colormap(jet);
c = colorbar;
c.Label.String = 'Potencia (dBm/dB)';
clim(rangoColores);
xlabel('Frecuencia (MHz)');
ylabel('Tiempo (s)');
title(['Espectrograma 2D: ' file], 'Interpreter', 'none');

% --- 3D Waterfall ---
fig3D = figure('Name', ['3D: ' file], 'Color', 'w');
ax3D = axes;
surf(ax3D, F/1e6, T, P.', 'EdgeColor', 'none');
view(-45, 60);
axis tight;
grid on;
colormap(jet);
colorbar;
clim(rangoColores);
zlim(rangoColores);
xlabel('Frecuencia (MHz)');
ylabel('Tiempo (s)');
zlabel('Potencia (dBm)');
title(['Waterfall 3D: ' file], 'Interpreter', 'none');

disp('¡Visualización lista!');
