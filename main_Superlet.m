if ~exist('run_batch_mode', 'var')
    clc; clear; close all;
    run_batch_mode = false;
end

%% --- CONFIGURACIÓN GENERAL ---
fs = 1e6;
f_HI = 1420e6;
carpeta = 'ArchivosIQ';

% Configuración Superlet
metodo = 'faslt';
baseCycles = 3;
srord = [5 10];
df = 2000;
max_samples_slt = 5000;

% Configuración Visualización
anchoBandaDefault = 400e3;
rangoColores = [-70 -30];
alinearRuido = false;
nivelRuidoObjetivo = 0;
umbral_guardado_dBm = -50;
offset_calibracion = -130; % Calibración basada en 'RadioTelescopio'

%% 1. SELECCIÓN DE ARCHIVO
if ~exist(carpeta, 'dir'), error(['La carpeta "' carpeta '" no existe.']); end
archivos = dir(fullfile(carpeta, '*.iq'));
if isempty(archivos), error('No se encontraron archivos .iq.'); end

disp('--- Archivos encontrados (SUPERLET) ---');
for i = 1:length(archivos)
    fprintf('%d: %s (%.2f MB)\n', i, archivos(i).name, archivos(i).bytes/1024/1024);
end
fprintf('\n');

if ~run_batch_mode
    idx = input('Selecciona el número del archivo: ');
    if isempty(idx) || idx < 1 || idx > length(archivos), error('Selección inválida.'); end
else
    idx = batch_idx;
end
archivoSeleccionado = archivos(idx);

fullPath = fullfile(carpeta, archivoSeleccionado.name);
[~, name_no_ext, ~] = fileparts(archivoSeleccionado.name);
timestamp_inicio = datestr(now, 'yyyymmdd_HHMMSS');
% dir_base_img definido más abajo
dir_base_dat = fullfile('Resultados_Datos', ['SUPERLET_' name_no_ext]);
if ~exist(dir_base_dat, 'dir'), mkdir(dir_base_dat); end

fc = f_HI;

%% 2. LECTURA DE DATOS
disp(' ');
if ~run_batch_mode
    inicio_pct = input('Inicio (%): '); if isempty(inicio_pct), inicio_pct=0; end
    fin_pct = input('Fin (%): '); if isempty(fin_pct), fin_pct=100; end
else
    inicio_pct = batch_inicio_pct;
    fin_pct = batch_fin_pct;
end
if fin_pct - inicio_pct > 100, fin_pct = inicio_pct + 100; end
pct_elegido = fin_pct - inicio_pct;

folder_suffix = sprintf('_Inicio%d_Fin%d', inicio_pct, fin_pct);
dir_base_img = fullfile('Resultados_Imagenes', ['SUPERLET_' name_no_ext], [timestamp_inicio folder_suffix]);
if ~exist(dir_base_img, 'dir'), mkdir(dir_base_img); end

totalMuestrasComplejas = archivoSeleccionado.bytes / 4;
byteInicio = floor((inicio_pct / 100) * totalMuestrasComplejas) * 4;
muestrasComplejasLeer = floor((pct_elegido / 100) * totalMuestrasComplejas);
muestrasShortLeer = muestrasComplejasLeer * 2;
timeOffset = (byteInicio / 4) / fs;

f = fopen(fullPath, 'r');
if f < 0, error('No se pudo abrir el archivo.'); end
fseek(f, byteInicio, 'bof');
s = fread(f, muestrasShortLeer, 'short=>single');
fclose(f);

y = s(1:2:end) + 1i*s(2:2:end);
% Normalización ELIMINADA (dBm reales)
% m = max(abs(y));
% if m > 0, y = y / m; end

%% 3. CONFIGURACIÓN DE FRECUENCIA Y DIEZMADO
anchoBanda = anchoBandaDefault;
f_min_req = fc - anchoBanda/2;
min_fs_needed = anchoBanda * 1.2;
max_decim_bw = floor(fs / min_fs_needed);
if max_decim_bw < 1, max_decim_bw = 1; end

decim_mem = ceil(length(y) / max_samples_slt);
factor_d = min(decim_mem, max_decim_bw);
if factor_d < 1, factor_d = 1; end

if factor_d > 1
    fprintf('--> Diezmando por factor %d...\n', factor_d);
    y = y(1:factor_d:end);
    fs = fs / factor_d;
end

t_rel = (0:length(y)-1) / fs;
t_signal = t_rel + timeOffset;

% Shift para Superlet (solo frecuencias positivas)
B = anchoBanda;
f_shift = fc - f_min_req;
y_shifted = y .* exp(1i*2*pi*f_shift*t_rel(:));

fois = 0:df:B;
if fois(end) ~= B, fois = [fois B]; end
fois(fois == 0) = [];
f_rf = f_min_req + fois;

%% 4. PROCESAMIENTO (SUPERLET)
disp('--> Iniciando procesamiento por bloques...');
blockSize = 5000;
overlap_block = blockSize/2;
step_block = blockSize - overlap_block;
numBlocks = floor((length(y) - overlap_block) / step_block);

% Inicializar Figuras
if ~run_batch_mode
    fig2D = figure('Color','w', 'Name', 'Superlet Viewer HI (2D)'); ax2D = axes; hold(ax2D, 'on');
    xlabel('Frecuencia (MHz)'); ylabel('Tiempo (s)');
    title(sprintf('Superlet | %s', name_no_ext), 'Interpreter', 'none');
    axis tight; colormap(jet); view(0, 90); clim(ax2D, rangoColores); colorbar(ax2D);
else
    ax2D = []; fig2D = [];
end

if ~run_batch_mode
    fig3D = figure('Color', 'w', 'Name', 'Waterfall 3D - Superlet'); ax3D = axes; hold(ax3D, 'on');
    xlabel('Frecuencia (MHz)'); ylabel('Tiempo (s)'); zlabel('Potencia (dBm)');
    title(['Waterfall 3D | ' archivoSeleccionado.name], 'Interpreter', 'none');
    view(-45, 60); axis tight; colormap(jet); grid on; clim(ax3D, rangoColores); zlim(ax3D, rangoColores); colorbar(ax3D);
else
    ax3D = []; fig3D = [];
end

mask_roi = f_rf < fc;
f_roi = f_rf(mask_roi);
lista_final_eventos = [];
sum_spec = zeros(length(f_rf), 1);
count_spec = 0;

for k = 1:numBlocks
    idx_start = (k-1)*step_block + 1;
    idx_end = idx_start + blockSize - 1;

    y_block = y(idx_start:idx_end);
    t_block_rel = t_rel(idx_start:idx_end);
    t_block_abs = t_signal(idx_start:idx_end);

    y_block_shifted = y_block .* exp(1i*2*pi*f_shift*t_block_rel(:));
    x_in = y_block_shifted(:).';

    try
        slt = aslt(x_in, fs, fois, baseCycles, srord, 0);
    catch
        x_in = real(y_block_shifted(:).');
        slt = aslt(x_in, fs, fois, baseCycles, srord, 0);
    end

    S = abs(slt);
    % S es amplitud (magnitud del complejo), por lo tanto potencia es 20*log10(S) o 10*log10(S^2)
    S_dBm = 20*log10(S + eps('single')) + offset_calibracion;

    % Calcular y mostrar ruido siempre
    noiseFloor = median(S_dBm(:));
    fprintf('   > Nivel de Ruido Detectado (Bloque %d): %.2f dBm\n', k, noiseFloor);

    if alinearRuido
        S_dBm = S_dBm + (nivelRuidoObjetivo - noiseFloor);
    end

    sum_spec = sum_spec + sum(S_dBm, 2);
    count_spec = count_spec + size(S_dBm, 2);

    % Métricas Robustas (Solo ROI)
    if any(mask_roi)
        S_dBm_roi = S_dBm(mask_roi, :);
        m_bloque = calcularMetricas(f_roi, S_dBm_roi, t_block_abs, []);

        if ~isempty(m_bloque)
            if m_bloque.P_max_dBm > umbral_guardado_dBm
                lista_final_eventos = [lista_final_eventos, m_bloque];
            end
        end
    end

    % Graficar (Downsampling)
    puntos_plot = 200;
    step_plot = max(1, floor(length(t_block_abs) / puntos_plot));
    idx_plot = 1:step_plot:length(t_block_abs);
    t_plot = t_block_abs(idx_plot);
    S_dBm_plot = S_dBm(:, idx_plot);
    mask_plot_bw = abs(f_rf - fc) <= anchoBanda/2;

    if ~run_batch_mode
        surf(ax2D, f_rf(mask_plot_bw)/1e6, t_plot, S_dBm_plot(mask_plot_bw, :).', 'EdgeColor', 'none');
        surf(ax3D, f_rf(mask_plot_bw)/1e6, t_plot, S_dBm_plot(mask_plot_bw, :).', 'EdgeColor', 'none');
        drawnow limitrate;
    end
    if mod(k,10)==0, fprintf(' Bloque %d/%.0f (%.0f%%)\n', k, numBlocks, k/numBlocks*100); end
end

%% 5. RESULTADOS FINALES
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

if ~run_batch_mode
    saveas(fig2D, fullfile(dir_base_img, ['SUPERLET_2D_SinMarcadores_' name_no_ext '.png']));
    saveas(fig2D, fullfile(dir_base_img, ['SUPERLET_2D_SinMarcadores_' name_no_ext '.fig']));
    saveas(fig3D, fullfile(dir_base_img, ['SUPERLET_3D_SinMarcadores_' name_no_ext '.png']));
    saveas(fig3D, fullfile(dir_base_img, ['SUPERLET_3D_SinMarcadores_' name_no_ext '.fig']));
end

if ~isempty(lista_final_eventos)
    fprintf('\n--- RESULTADOS ESTADÍSTICOS (SUPERLET) ---\n');
    fprintf('Eventos Totales: %d\n', length(lista_final_eventos));

    if ~isempty(lista_final_eventos)
        fprintf('\n>>> MÉTRICAS ROBUSTAS (Promedio General) <<<\n');
        fprintf('  - Diferencia Pico-Ruido: %.2f dBm\n', mean([lista_final_eventos.Diferencia_Potencia_dBm]));
        fprintf('  - Delta Estabilidad:     %.2f dBm\n', mean([lista_final_eventos.Delta_Estabilidad_dBm]));
    end

    nombreExcel = fullfile(dir_base_dat, ['Reporte_Metricas_SUPERLET_' datestr(now, 'yyyymmdd_HHMMSS') '.xlsx']);
    registrarMetrica('Superlets', lista_final_eventos, nombreExcel);

    [~, idx_sort] = sort([lista_final_eventos.Diferencia_Potencia_dBm], 'descend');
    eventos_sorted = lista_final_eventos(idx_sort);

    for p=1:length(eventos_sorted)
        pk=eventos_sorted(p);
        if 1, c_m = 'm'; else, c_m = 'c'; end
        if ~run_batch_mode
            plot3(ax2D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, 200, 'v', 'MarkerSize', 10, 'MarkerFaceColor', c_m, 'MarkerEdgeColor', 'k');
            plot3(ax3D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, 200, 'v', 'MarkerSize', 10, 'MarkerFaceColor', c_m, 'MarkerEdgeColor', 'k');
        end
    end
end

if ~run_batch_mode
    saveas(fig2D, fullfile(dir_base_img, ['SUPERLET_2D_ConMarcadores_' name_no_ext '.png']));
    saveas(fig2D, fullfile(dir_base_img, ['SUPERLET_2D_ConMarcadores_' name_no_ext '.fig']));
    saveas(fig3D, fullfile(dir_base_img, ['SUPERLET_3D_ConMarcadores_' name_no_ext '.png']));
    saveas(fig3D, fullfile(dir_base_img, ['SUPERLET_3D_ConMarcadores_' name_no_ext '.fig']));
end

spec_profile = sum_spec / max(count_spec, 1);
if ~run_batch_mode
    figProfile = figure('Name', 'Perfil Promedio', 'Color', 'w');
    plot(f_rf/1e6, spec_profile, 'k', 'LineWidth', 1.2);
    grid on; xlabel('Frecuencia (MHz)'); ylabel('Amplitud Promedio (dBm)');
    xlim([f_min_req, f_min_req + anchoBanda]/1e6);
    xline(f_HI/1e6, 'r--', 'LineWidth', 1);
    title('Estimación Espectral (Superlet)');
    saveas(figProfile, fullfile(dir_base_img, ['SUPERLET_Perfil_' name_no_ext '.png']));
    saveas(figProfile, fullfile(dir_base_img, ['SUPERLET_Perfil_' name_no_ext '.fig']));
end

disp('Proceso Superlet completado.');

if run_batch_mode
    if ~exist('spec_profile', 'var')
        spec_profile = sum_spec / max(count_spec, 1);
    end
    results_batch.Superlet.f = f_rf;
    results_batch.Superlet.P = spec_profile;
end
