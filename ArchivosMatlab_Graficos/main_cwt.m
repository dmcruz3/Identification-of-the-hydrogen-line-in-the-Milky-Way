clc
clear
close all;

if ~exist('run_batch_mode', 'var')
    clc; clear; close all;
    run_batch_mode = false;
end

%% --- CONFIGURACIÓN GENERAL ---
fs = 1e6;
f_HI = 1420e6;
carpeta = 'ArchivosIQ';

% Configuración Visualización y Análisis
anchoBanda = 400e3;
% Rango Fijo de Visualización (MHz)
f_min_visual = 1419.8;
f_max_visual = 1420.2;
rangoColores = [-90 -40];
alinearRuido = true;
nivelRuidoObjetivo = -80;
umbral_guardado_dBm = -50;
offset_calibracion = -0; % Calibración basada en 'RadioTelescopio'

%% 1. SELECCIÓN DE ARCHIVO
if ~exist(carpeta, 'dir'), error(['La carpeta "' carpeta '" no existe.']); end
archivos = dir(fullfile(carpeta, '*.iq'));
if isempty(archivos), error('No se encontraron archivos .iq.'); end

disp('--- Archivos encontrados (CWT) ---');
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
dir_base_dat = fullfile('Resultados_Datos', ['CWT_' name_no_ext]);
if ~exist(dir_base_dat, 'dir'), mkdir(dir_base_dat); end
if ~exist(dir_base_dat, 'dir'), mkdir(dir_base_dat); end

fc = f_HI;

%% 2. LECTURA DE DATOS
disp(' ');
totalMuestras = archivoSeleccionado.bytes / 4;
duracionTotal = totalMuestras / fs;

if ~run_batch_mode
    disp(' ');
    fprintf('Duración Total del Archivo: %.2f segundos\n', duracionTotal);

    inicio_seg = input('Inicio (s): ');
    if isempty(inicio_seg), inicio_seg = 0; end

    fin_seg = input('Fin (s): ');
    if isempty(fin_seg), fin_seg = duracionTotal; end

    if fin_seg > duracionTotal, fin_seg = duracionTotal; end
    if inicio_seg < 0, inicio_seg = 0; end
    if fin_seg <= inicio_seg, fin_seg = duracionTotal; end

    duracion_elegida = fin_seg - inicio_seg;
    pct_usado = (duracion_elegida / duracionTotal) * 100;
    fprintf('--> Se analizará %.2f%% del archivo (%.2f s)\n', pct_usado, duracion_elegida);

    byteInicio = floor(inicio_seg * fs) * 4;
    muestrasLeer = floor(duracion_elegida * fs) * 2;

    folder_suffix = sprintf('_Inicio%.2fs_Fin%.2fs', inicio_seg, fin_seg);
    pct_display = pct_usado;
else
    inicio_pct = batch_inicio_pct;
    fin_pct = batch_fin_pct;
    if fin_pct - inicio_pct > 100, fin_pct = inicio_pct + 100; end
    pct_elegido = fin_pct - inicio_pct;

    byteInicio = floor((inicio_pct / 100) * totalMuestras) * 4;
    muestrasLeer = floor((pct_elegido / 100) * totalMuestras) * 2;

    folder_suffix = sprintf('_Inicio%d_Fin%d', inicio_pct, fin_pct);
    pct_display = pct_elegido;
end
dir_base_img = fullfile('Resultados_Imagenes', ['CWT_' name_no_ext], [timestamp_inicio folder_suffix]);
if ~exist(dir_base_img, 'dir'), mkdir(dir_base_img); end

% totalMuestras, byteInicio, muestrasLeer calculados arriba
timeOffset = (byteInicio / 4) / fs;

f = fopen(fullPath, 'r');
fseek(f, byteInicio, 'bof');
disp(['--> Leyendo ' num2str(pct_display) '%...']);
s = fread(f, muestrasLeer, 'short=>single');
fclose(f);

y = s(1:2:end) + 1i*s(2:2:end);
% Normalización ELIMINADA (dBm reales)
% m = max(abs(y));
% if m > 0, y = y / m; end

% --- OPTIMIZACIÓN DE MEMORIA ---
max_samples_cwt = 5000;
min_fs_needed = anchoBanda * 1.2;
max_decim_bw = floor(fs / min_fs_needed);
if max_decim_bw < 1, max_decim_bw = 1; end
decim_mem = ceil(length(y) / max_samples_cwt);
factor_d = min(decim_mem, max_decim_bw);

if factor_d > 1
    fprintf('Diezmando por factor %d (FS final: %.2f kHz).\n', factor_d, (fs/factor_d)/1e3);
    y = y(1:factor_d:end);
    fs = fs / factor_d;
end

t_signal = ((0:length(y)-1) / fs) + timeOffset;

%% 3. PROCESAMIENTO (CWT)
disp('--- Calculando CWT por Ventanas ---');
if exist('cwt', 'file') ~= 2, error('Wavelet Toolbox necesario.'); end

% Inicializar Figuras
% Límites de frecuencia para visualización centrada
f_min_vis = (fc - anchoBanda/2) / 1e6;
f_max_vis = (fc + anchoBanda/2) / 1e6;


if ~run_batch_mode
    fig3D = figure('Name', ' CWT (3D)', 'Color', 'w'); ax3D = axes; hold(ax3D, 'on');
    grid(ax3D, 'on'); xlabel(ax3D, 'Frecuencia (MHz)'); ylabel(ax3D, 'Tiempo (s)'); zlabel(ax3D, 'Magnitud (dBm)');
    title(sprintf(' CWT | %s', name_no_ext), 'Interpreter', 'none');
    clim(ax3D, rangoColores); zlim(ax3D, rangoColores); view(ax3D, -45, 60); colormap(fig3D, jet); colorbar(ax3D);
    xlim(ax3D, [f_min_visual f_max_visual]);
else
    ax3D = []; fig3D = [];
end

blockSize = 5000;
overlap_block = blockSize/2;
step_block = blockSize - overlap_block;
numBlocks = floor((length(y) - overlap_block) / step_block);

lista_final_eventos = [];
sum_spec = [];
count_spec = 0;
f_profile = [];

for k = 1:numBlocks
    idx_start = (k-1)*step_block + 1;
    idx_end = idx_start + blockSize - 1;

    y_block = y(idx_start:idx_end);
    t_block = t_signal(idx_start:idx_end);

    fprintf('Bloque %d/%d (%.1f%%)...\n', k, numBlocks, (k/numBlocks)*100);

    [cfs, f] = cwt(y_block, fs, 'VoicesPerOctave', 16);

    puntos_plot = 200;
    step_plot = max(1, floor(length(t_block) / puntos_plot));
    idx_plot = 1:step_plot:length(t_block);
    t_plot_vals = t_block(idx_plot);

    mag_for_metrics = [];   % Magnitud para cálculo de métricas
    f_for_metrics = [];     % Frecuencia para cálculo de métricas

    if ndims(cfs) == 3
        % Señal Compleja
        % Lado Izquierdo (f < fc): cfs(:,:,2) -> Frecuencias Negativas
        mag_left = 20*log10(abs(cfs(:,:,2))) + offset_calibracion;
        f_left = fc - f;

        % Lado Derecho (f > fc): cfs(:,:,1) -> Frecuencias Positivas
        mag_right = 20*log10(abs(cfs(:,:,1))) + offset_calibracion;
        f_right = fc + f;

        % Unificar Datos
        f_combined = [f_left; f_right];
        mag_combined = [mag_left; mag_right];

        % Ordenar
        [f_sorted, idx_sort] = sort(f_combined);
        mag_sorted = mag_combined(idx_sort, :);

        % Calcular Ruido
        noise_val = median(mag_sorted(:));
        fprintf('   > Nivel de Ruido Detectado (Bloque %d): %.2f dBm\n', k, noise_val);

        if alinearRuido
            offset = nivelRuidoObjetivo - noise_val;
            mag_sorted = mag_sorted + offset;
        end

        f_axis = f_sorted;
        mag_for_metrics = mag_sorted;
        f_for_metrics = f_axis;

        % Graficar - Filtrar por ancho de banda para no cargar gráfica
        % Usamos anchoBanda efectivo centrado en fc.
        % f_axis está en Hz.
        % Graficar - Usar los límites visuales definidos por el usuario
        % f_axis está en Hz. f_min_visual/f_max_visual en MHz.
        idx_vis = f_axis >= f_min_visual*1e6 & f_axis <= f_max_visual*1e6;
        f_vis = f_axis(idx_vis);

        if ~run_batch_mode && any(idx_vis)
            % surf(ax2D, f_vis/1e6, t_plot_vals, mag_sorted(idx_vis, idx_plot).', 'EdgeColor', 'none');
            surf(ax3D, f_vis/1e6, t_plot_vals, mag_sorted(idx_vis, idx_plot).', 'EdgeColor', 'none');
        end

    else
        % Señal Real
        mag = 20*log10(abs(cfs)) + offset_calibracion;
        % Calcular y mostrar ruido siempre
        noise_val = median(mag(:));
        fprintf('   > Nivel de Ruido Detectado (Bloque %d): %.2f dBm\n', k, noise_val);

        if alinearRuido
            offset = nivelRuidoObjetivo - noise_val;
            mag = mag + offset;
        end
        f_axis = fc + f;

        mag_for_metrics = mag;
        f_for_metrics = f_axis;

        % Graficar
        % Graficar - Usar los límites visuales definidos por el usuario
        idx_vis = f_axis >= f_min_visual*1e6 & f_axis <= f_max_visual*1e6;
        f_vis = f_axis(idx_vis);

        if ~run_batch_mode && any(idx_vis)
            surf(ax2D, f_vis/1e6, t_plot_vals, mag(idx_vis, idx_plot).', 'EdgeColor', 'none');
            surf(ax3D, f_vis/1e6, t_plot_vals, mag(idx_vis, idx_plot).', 'EdgeColor', 'none');
        end
    end

    % --- MÉTRICAS ROBUSTAS CWT ---
    % Tomamos el espectrograma completo del bloque
    % Ajuste: Convertir a lineal (Amplitud^2 o Potencia) para calcularMetricas.m
    % mag_for_metrics está en dBm (20*log10) o dB (20*log10 o 10*log10).
    % Asumimos que queremos magnitud de potencia lineal para SNR.
    P_lin_metrics = 10.^(mag_for_metrics./10);

    m_bloque = calcularMetricas(f_for_metrics, P_lin_metrics, t_block, []);

    if ~isempty(m_bloque)
        if m_bloque.P_max_dBm > umbral_guardado_dBm
            lista_final_eventos = [lista_final_eventos, m_bloque];
        end
    end

    % Acumular 1D (usando el centroide del bloque)
    % --- CÁLCULO DE PERFIL (Linear Average) ---
    % Convertimos a lineal si no lo está (CWT usa dBm)
    % mag_for_metrics ya está en dBm (20*log10), así que volvemos a amplitud lineal
    Amp_lin = 10.^(mag_for_metrics./20); % O 10^(dB/20) para amplitud, o 10^(dB/10) para potencia.
    % La métrica visual suele ser potencia. Usaremos 10^(dB/10) asumiendo que queremos Promedio de Potencia.
    P_lin = 10.^(mag_for_metrics./10);

    current_sum = sum(P_lin, 2);

    if isempty(sum_spec)
        sum_spec = current_sum;
        f_profile = f_for_metrics;
        count_spec = size(mag_for_metrics, 2);
    elseif size(current_sum, 1) == size(sum_spec, 1)
        sum_spec = sum_spec + current_sum;
        count_spec = count_spec + size(mag_for_metrics, 2);
    end

    clear cfs mag_neg mag;
    drawnow limitrate;
end

%% 4. RESULTADOS FINALES
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

if ~run_batch_mode
    % Asegurar límites antes de guardar "final de analisis"
    if ~isempty(ax3D) && isvalid(ax3D)
        xlim(ax3D, [f_min_visual f_max_visual]);
    end
    % saveas(fig2D, fullfile(dir_base_img, ['CWT_2D_SinMarcadores_' name_no_ext '.png']));
    % saveas(fig2D, fullfile(dir_base_img, ['CWT_2D_SinMarcadores_' name_no_ext '.fig']));
    saveas(fig3D, fullfile(dir_base_img, ['CWT_3D_SinMarcadores_' name_no_ext '.png']));
    saveas(fig3D, fullfile(dir_base_img, ['CWT_3D_SinMarcadores_' name_no_ext '.fig']));
end

if ~isempty(lista_final_eventos)
    fprintf('\n--- RESULTADOS ESTADÍSTICOS (CWT) ---\n');
    fprintf('Eventos Totales: %d\n', length(lista_final_eventos));

    if ~isempty(lista_final_eventos)
        fprintf('\n>>> MÉTRICAS ROBUSTAS (Promedio General) <<<\n');
        fprintf('  - Diferencia Pico-Ruido: %.2f dBm\n', mean([lista_final_eventos.Diferencia_Potencia_dBm]));
        fprintf('  - Delta Estabilidad:     %.2f dBm\n', mean([lista_final_eventos.Delta_Estabilidad_dBm]));
    end

    nombreExcel = fullfile(dir_base_dat, ['Reporte_Metricas_CWT_' datestr(now, 'yyyymmdd_HHMMSS') '.xlsx']);
    registrarMetrica('CWT', lista_final_eventos, nombreExcel);

    % Marcadores
    [~, idx_sort] = sort([lista_final_eventos.Diferencia_Potencia_dBm], 'descend');
    eventos_sorted = lista_final_eventos(idx_sort);

    for p = 1:length(eventos_sorted)
        pk = eventos_sorted(p);
        z_mark = pk.P_max_dBm;
        if 1, c='m'; else, c='c'; end
        if ~run_batch_mode
            plot3(ax3D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, z_mark, 'v', 'MarkerFaceColor',c,'MarkerEdgeColor','k');
            % plot3(ax2D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, 200, 'v', 'MarkerFaceColor',c,'MarkerEdgeColor','k');
        end
    end
end


if ~run_batch_mode
    % saveas(fig2D, fullfile(dir_base_img, ['CWT_2D_ConMarcadores_' name_no_ext '.png']));
    % saveas(fig2D, fullfile(dir_base_img, ['CWT_2D_ConMarcadores_' name_no_ext '.fig']));
    saveas(fig3D, fullfile(dir_base_img, ['CWT_3D_ConMarcadores_' name_no_ext '.png']));
    saveas(fig3D, fullfile(dir_base_img, ['CWT_3D_ConMarcadores_' name_no_ext '.fig']));
end

if ~isempty(sum_spec)
    spec_profile_lin = sum_spec / max(count_spec, 1);
    spec_profile = 10*log10(spec_profile_lin);
    if ~run_batch_mode
        figProfile = figure('Name', 'Perfil Espectral Promedio', 'Color', 'w');
        [f_sorted, idx_s] = sort(f_profile);
        spec_sorted = spec_profile(idx_s);
        plot(f_sorted/1e6, spec_sorted, 'k', 'LineWidth', 1.2);
        grid on; xlabel('Frecuencia (MHz)'); ylabel('Amplitud Promedio (dBm)');
        xline(fc/1e6, 'r--', 'LineWidth', 1);
        title('Estimación Espectral (CWT)');
        xlim([f_min_vis f_max_vis]);
        saveas(figProfile, fullfile(dir_base_img, ['CWT_Perfil_' name_no_ext '.png']));
        saveas(figProfile, fullfile(dir_base_img, ['CWT_Perfil_' name_no_ext '.fig']));
    end
end

disp('Proceso CWT Completado.');

if run_batch_mode
    if ~exist('spec_profile', 'var')
        spec_profile_lin = sum_spec / max(count_spec, 1);
        spec_profile = 10*log10(spec_profile_lin);
    end
    results_batch.CWT.f = f_profile;
    results_batch.CWT.P = spec_profile;
end
