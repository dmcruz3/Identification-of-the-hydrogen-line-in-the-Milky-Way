if ~exist('run_batch_mode', 'var')
    clc; clear; close all;
    run_batch_mode = false;
end

%% --- CONFIGURACIÓN GENERAL ---
fs = 1e6;
f_HI = 1420e6;
carpeta = 'ArchivosIQ';

% Configuración STFT
nfft = 1024;
window_size = 128;
overlap = 64;
nombre_ventana = 'Blackman';

umbral_guardado_dBm = -75;
offset_calibracion = -145;
anchoBanda = 400e3;
rangoColores = [-90 -60];
alinearRuido = false;
nivelRuidoObjetivo = 0;


%% 1. SELECCIÓN DE ARCHIVO
if ~exist(carpeta, 'dir'), error(['La carpeta "' carpeta '" no existe.']); end
archivos = dir(fullfile(carpeta, '*.iq'));
if isempty(archivos), error('No se encontraron archivos .iq.'); end

disp('--- Archivos encontrados ---');
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
% dir_base_img se definirá después de obtener los porcentajes
dir_base_dat = fullfile('Resultados_Datos', ['DIRECTO_' name_no_ext]);
if ~exist(dir_base_dat, 'dir'), mkdir(dir_base_dat); end

fc = f_HI;

%% 2. LECTURA DE DATOS
if ~run_batch_mode
    disp(' ');
    inicio_pct = input('Inicio (%): '); if isempty(inicio_pct), inicio_pct=0; end
    fin_pct = input('Fin (%): '); if isempty(fin_pct), fin_pct=100; end
else
    inicio_pct = batch_inicio_pct;
    fin_pct = batch_fin_pct;
end
pct_elegido = fin_pct - inicio_pct;

folder_suffix = sprintf('_Inicio%d_Fin%d', inicio_pct, fin_pct);
dir_base_img = fullfile('Resultados_Imagenes', ['DIRECTO_' name_no_ext], [timestamp_inicio folder_suffix]);
if ~exist(dir_base_img, 'dir'), mkdir(dir_base_img); end

totalMuestras = archivoSeleccionado.bytes / 4;
byteInicio = floor((inicio_pct / 100) * totalMuestras) * 4;
muestrasLeer = floor((pct_elegido / 100) * totalMuestras) * 2;
timeOffset = (byteInicio / 4) / fs;

f = fopen(fullPath, 'r');
fseek(f, byteInicio, 'bof');
disp(['--> Leyendo ' num2str(pct_elegido) '%...']);
s = fread(f, muestrasLeer, 'short=>single');
fclose(f);

y = s(1:2:end) + 1i*s(2:2:end);
% Normalización ELIMINADA para mantener unidades relativas/reales (dBm)
% m = max(abs(y));
% if m > 0
%     y = y / m;
% end
t_signal = ((0:length(y)-1) / fs) + timeOffset;

%% 3. PROCESAMIENTO (STFT)
disp('--- Iniciando STFT por Bloques ---');

blockSize = 5000;
overlap_block = blockSize/2;
step_block = blockSize - overlap_block;
numBlocks = floor((length(y) - overlap_block) / step_block);

switch lower(nombre_ventana)
    case 'blackman', win = blackman(window_size);
    case 'flattop', win = flattopwin(window_size);
    otherwise, win = hamming(window_size);
end

% Inicializar Figuras
if ~run_batch_mode
    fig2D = figure('Name', 'Espectro Directo 2D', 'Color', 'w'); ax2D = axes; hold(ax2D,'on');
    xlabel('Frecuencia (MHz)'); ylabel('Tiempo (s)');
    title(sprintf('Analisis Directo | %s', archivoSeleccionado.name), 'Interpreter', 'none');
    clim(ax2D, rangoColores); colormap(jet); colorbar; view(0,90);

    fig3D = figure('Name', 'Waterfall Directo 3D', 'Color', 'w'); ax3D = axes; hold(ax3D,'on');
    xlabel('Frecuencia (MHz)'); ylabel('Tiempo (s)');
    zlabel('Potencia (dBm)');
    title(sprintf('Waterfall 3D | %s', archivoSeleccionado.name), 'Interpreter', 'none');
    view(-45, 60); grid on; clim(ax3D, rangoColores); zlim(ax3D, rangoColores); colormap(jet); colorbar;
end

lista_final_eventos = [];

% Eje FFT
f = (-nfft/2 : nfft/2-1) * (fs / nfft);
f_abs = fc + f;

sum_spec = zeros(nfft, 1);
count_spec = 0;

for k = 1:numBlocks
    idx_start = (k-1)*step_block + 1;
    idx_end = idx_start + blockSize - 1;

    y_block = y(idx_start:idx_end);
    t_start_blk = t_signal(idx_start);

    [S, F_stft, T_stft] = spectrogram(y_block, win, overlap, nfft, fs, 'centered');

    T_abs = T_stft + t_start_blk;
    F_abs = fc + F_stft;

    % --- Modo Relativo (dBm) ---
    P_dBm = 10*log10(abs(S).^2 + eps) + offset_calibracion;

    % Calcular y mostrar ruido siempre
    ruido_est = median(P_dBm(:));
    fprintf('   > Nivel de Ruido Detectado (Bloque %d): %.2f dBm\n', k, ruido_est);

    if alinearRuido
        P_dBm = P_dBm + (nivelRuidoObjetivo - ruido_est);
    end

    sum_spec = sum_spec + sum(P_dBm, 2);
    count_spec = count_spec + size(P_dBm, 2);

    % --- CÁLCULO DE MÉTRICAS ROBUSTAS ---
    % Solo analizamos la parte izquierda (f < 0 localmente, < fc absoluto) si se desea HI
    mask_roi = F_stft < 0;
    P_roi = P_dBm(mask_roi, :);
    f_roi_curr = F_abs(mask_roi);

    m_bloque = calcularMetricas(f_roi_curr, P_roi, T_abs, []);

    if ~isempty(m_bloque)
        % Filtro Previo (Nivel > Umbral) para no llenar de ruido
        if m_bloque.P_max_dBm > umbral_guardado_dBm
            lista_final_eventos = [lista_final_eventos, m_bloque];
        end
    end

    % Graficar (Downsampling)
    idx_vis = abs(F_stft) <= (anchoBanda/2);
    F_vis = F_abs(idx_vis);
    P_vis = P_dBm(idx_vis, :);

    step_plot = max(1, floor(length(T_abs)/50));
    idx_t_plot = 1:step_plot:length(T_abs);
    T_plot = T_abs(idx_t_plot);
    P_plot = P_vis(:, idx_t_plot);

    if ~run_batch_mode
        surf(ax2D, F_vis/1e6, T_plot, P_plot.', 'EdgeColor', 'none');
        surf(ax3D, F_vis/1e6, T_plot, P_plot.', 'EdgeColor', 'none');
        drawnow limitrate;
    end

    if mod(k, 10) == 0
        %fprintf('Bloque %d/%d (%.1f%%)\n', k, numBlocks, 100*k/numBlocks);
    end
end

%% 4. RESULTADOS FINALES
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

if ~run_batch_mode
    saveas(fig2D, fullfile(dir_base_img, ['Directo_2D_SinMarcadores_' name_no_ext '.png']));
    saveas(fig2D, fullfile(dir_base_img, ['Directo_2D_SinMarcadores_' name_no_ext '.fig']));
    saveas(fig3D, fullfile(dir_base_img, ['Directo_3D_SinMarcadores_' name_no_ext '.png']));
    saveas(fig3D, fullfile(dir_base_img, ['Directo_3D_SinMarcadores_' name_no_ext '.fig']));

    if ~isempty(lista_final_eventos)
        fprintf('\n--- RESULTADOS ESTADÍSTICOS ---\n');
        fprintf('Eventos Totales Detectados: %d\n', length(lista_final_eventos));

        % Reporte de Métricas Promedio (Todos los eventos)
        if ~isempty(lista_final_eventos)
            fprintf('\n>>> MÉTRICAS ROBUSTAS (Promedio General) <<<\n');
            fprintf('  - Diferencia Pico-Ruido: %.2f dBm\n', mean([lista_final_eventos.Diferencia_Potencia_dBm]));
            fprintf('  - Potencia Máxima:       %.2f dBm\n', mean([lista_final_eventos.P_max_dBm]));
            fprintf('  - Piso de Ruido:         %.2f dBm\n', mean([lista_final_eventos.P_ruido_dBm]));
            fprintf('  - Potencia Promedio:     %.2f dBm\n', mean([lista_final_eventos.P_promedio_dBm]));
            fprintf('  - Delta Estabilidad:     %.2f dBm\n', mean([lista_final_eventos.Delta_Estabilidad_dBm]));
        end

        % Reporte de Métricas Promedio (Estables)


        % Marcadores en Gráficas
        [~, idx_sort] = sort([lista_final_eventos.Diferencia_Potencia_dBm], 'descend');
        eventos_sorted = lista_final_eventos(idx_sort);

        for p_idx = 1:length(eventos_sorted)
            pk = eventos_sorted(p_idx);
            z_mark = min(max(pk.P_max_dBm, rangoColores(1)), rangoColores(2));

            % Color único para todos, ya no validamos estabilidad
            color_mark = 'm';


            plot3(ax3D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, z_mark, 'v', 'MarkerSize', 8, ...
                'MarkerFaceColor', color_mark, 'MarkerEdgeColor','k');
            plot3(ax2D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, 200, 'v', 'MarkerSize', 8, ...
                'MarkerFaceColor', color_mark, 'MarkerEdgeColor','k');

            if p_idx <= 5
                text(ax3D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, z_mark, sprintf(' Dif: %.1f', pk.Diferencia_Potencia_dBm), ...
                    'Color', color_mark, 'FontSize',8, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 1);
            end
        end

        nombreExcel = fullfile(dir_base_dat, ['Reporte_Metricas_Directo_' datestr(now, 'yyyymmdd_HHMMSS') '.xlsx']);
        registrarMetrica('Directo', lista_final_eventos, nombreExcel);
    end

    saveas(fig2D, fullfile(dir_base_img, ['Directo_2D_ConMarcadores_' name_no_ext '.png']));
    saveas(fig2D, fullfile(dir_base_img, ['Directo_2D_ConMarcadores_' name_no_ext '.fig']));
    saveas(fig3D, fullfile(dir_base_img, ['Directo_3D_ConMarcadores_' name_no_ext '.png']));
    saveas(fig3D, fullfile(dir_base_img, ['Directo_3D_ConMarcadores_' name_no_ext '.fig']));

    % Perfil 1D
    figProfile = figure('Name', 'Perfil Promedio', 'Color', 'w');
    spec_profile = sum_spec / max(count_spec, 1);
    plot(f_abs/1e6, spec_profile, 'k', 'LineWidth', 1.2);
    grid on; xlabel('Frecuencia (MHz)');
    ylabel('Amplitud Promedio (dBm)');
    xline(f_HI/1e6, 'r--', 'LineWidth', 1);
    title('Estimación Espectral (Directo)');
    saveas(figProfile, fullfile(dir_base_img, ['Directo_Perfil_' name_no_ext '.png']));
    saveas(figProfile, fullfile(dir_base_img, ['Directo_Perfil_' name_no_ext '.fig']));
end

if run_batch_mode
    % Recalcular spec_profile por si no entró al bloque de ploteo
    if ~exist('spec_profile', 'var')
        spec_profile = sum_spec / max(count_spec, 1);
    end
    results_batch.Directo.f = f_abs;
    results_batch.Directo.P = spec_profile;
end

disp('Proceso Directo completado.');
