clc
clear all;
close all;

if ~exist('run_batch_mode', 'var')
    clc; clear; close all;
    run_batch_mode = false;
end

%% --- CONFIGURACIÓN GENERAL ---
fs = 1e6;
f_HI = 1420e6;
carpeta = 'ArchivosIQ';

% Configuración Correlación
nfft = 1024;
max_lag = 37;
nombre_ventana_lag = 'bartlett';
window_len = 128;
overlap = window_len/2;

% Configuración Visualización
anchoBanda = 400e3;
rangoColores = [-90 -40];
alinearRuido = false;
umbral_guardado_dBm = -70;
offset_calibracion = -120; % Calibración basada en 'RadioTelescopio'

%% 1. SELECCIÓN DE ARCHIVO
if ~exist(carpeta, 'dir'), error(['La carpeta "' carpeta '" no existe.']); end
archivos = dir(fullfile(carpeta, '*.iq'));
if isempty(archivos), error('No se encontraron archivos .iq.'); end

disp('--- Archivos (INDIRECTO) ---');
for i=1:length(archivos), fprintf('%d: %s \n', i, archivos(i).name); end
fprintf('\n');

if ~run_batch_mode
    idx = input('Selección: ');
    if isempty(idx) || idx < 1 || idx > length(archivos), error('Selección inválida.'); end
else
    idx = batch_idx;
end
archivoSeleccionado = archivos(idx);
fullPath = fullfile(carpeta, archivoSeleccionado.name);

[~, name_no_ext, ~] = fileparts(archivoSeleccionado.name);
timestamp_inicio = datestr(now, 'yyyymmdd_HHMMSS');
dir_base_dat = fullfile('Resultados_Datos', ['INDIRECTO_' name_no_ext]);
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
pct = fin_pct-inicio_pct;

folder_suffix = sprintf('_Inicio%d_Fin%d', inicio_pct, fin_pct);
dir_base_img = fullfile('Resultados_Imagenes', ['INDIRECTO_' name_no_ext], [timestamp_inicio folder_suffix]);
if ~exist(dir_base_img, 'dir'), mkdir(dir_base_img); end

totalB = archivoSeleccionado.bytes;
byteInit = floor((inicio_pct/100)*totalB/4)*4;
nRead = floor((pct/100)*totalB/4)*2;
timeOff = (byteInit/4)/fs;

f = fopen(fullPath,'r');
fseek(f,byteInit,'bof');
s = fread(f,nRead,'short=>single');
fclose(f);

y = s(1:2:end)+1i*s(2:2:end);
% Normalización ELIMINADA (dBm reales)
% m = max(abs(y));
% if m>0, y=y/m; end
t_signal = (0:length(y)-1)/fs + timeOff;

%% 3. PROCESAMIENTO (INDIRECTO)
disp('--- Iniciando Método Indirecto ---');
blockSize = 5000;
overlap_block = blockSize/2;
step_block = blockSize - overlap_block;
numBlocks = floor((length(y) - overlap_block) / step_block);

w_lag = bartlett(2*max_lag+1);

% Figuras: Se crean al final para evitar bloqueo
if ~run_batch_mode
    fig2D=figure('Name', 'Espectro Indirecto 2D', 'Color','w'); ax2D=axes; hold(ax2D,'on');
    title(sprintf('Indirecto 2D | %s', name_no_ext), 'Interpreter', 'none');
    view(0,90); clim(ax2D,rangoColores); colormap(jet); colorbar; xlabel('Frecuencia (MHz)'); ylabel('Tiempo (s)');

    fig3D=figure('Name', 'Waterfall Indirecto 3D', 'Color','w'); ax3D=axes; hold(ax3D,'on');
    title(sprintf('Waterfall Indirecto | %s', name_no_ext), 'Interpreter', 'none');
    view(-45,60); grid on; clim(ax3D,rangoColores); zlim(ax3D,rangoColores); colormap(jet); colorbar; xlabel('Frecuencia (MHz)'); ylabel('Tiempo (s)'); zlabel('Potencia (dBm)');
else
    ax2D=[]; fig2D=[]; ax3D=[]; fig3D=[];
end

lista_final_eventos = [];

f_vec = (-nfft/2:nfft/2-1)*(fs/nfft);
f_abs = fc + f_vec;

sum_spec = zeros(nfft, 1);
count_spec = 0;

% Acumuladores para visualización final
Accum_Freq = f_abs/1e6; % Vector de frecuencias (constante)
Accum_Time = [];        % Se irá concatenando
Accum_P    = [];        % Se irá concatenando

for k=1:numBlocks
    idx_s = (k-1)*step_block+1;
    idx_e = idx_s+blockSize-1;
    y_blk = y(idx_s:idx_e);
    t_start = t_signal(idx_s);

    [segments, ~] = buffer(y_blk, window_len, overlap, 'nodelay');
    [~, num_seg] = size(segments);
    t_seg_abs = t_start + ((0:num_seg-1)*(window_len-overlap) + window_len/2)/fs;

    P_blk = zeros(nfft, num_seg);

    for i=1:num_seg
        seg = segments(:,i);
        [R, ~] = xcorr(seg, max_lag, 'biased');
        R_w = R .* w_lag;
        P = fftshift(abs(fft(R_w, nfft)));
        P_blk(:,i) = P;
    end

    P_dBm = 10*log10(P_blk + eps) + offset_calibracion;

    % Calcular y mostrar ruido siempre
    r_est = median(P_dBm(:));
    fprintf('   > Nivel de Ruido Detectado (Bloque %d): %.2f dBm\n', k, r_est);

    if alinearRuido
        P_dBm = P_dBm + (0 - r_est);
    end

    sum_spec = sum_spec + sum(P_dBm, 2);
    count_spec = count_spec + size(P_dBm, 2);

    % Métricas Robustas
    mask_roi = f_vec < 0;
    P_roi = P_dBm(mask_roi,:);
    f_roi_curr = f_abs(mask_roi);

    m_bloque = calcularMetricas(f_roi_curr, P_roi, t_seg_abs, []);

    if ~isempty(m_bloque)
        if m_bloque.P_max_dBm > umbral_guardado_dBm
            lista_final_eventos = [lista_final_eventos, m_bloque];
        end
    end

    % Acumular datos para gráfico final
    mask_vis = abs(f_vec) <= anchoBanda/2;
    % Transponemos P_dBm para que sea (Tiempo x Frecuencia) si usamos surf(X,Y,Z) donde Z es Matriz
    % O mantenemos (Frecuencia x Tiempo) y usamos 'shading flat'.
    % Normalmente waterfall espera Z(y,x).
    % Aquí acumularemos tal cual sale: Freq x Time

    % Para simplificar la concatenación temporal:
    Accum_Time = [Accum_Time, t_seg_abs];
    Accum_P    = [Accum_P, P_dBm(mask_vis,:)];

    if mod(k,10)==0, fprintf(' Bloque %d/%.0f (%.0f%%)\n', k, numBlocks, k/numBlocks*100); end
end

%% 4. RESULTADOS FINALES Y EXCEL (Prioritario)
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% --- GUARDADO DE DATOS CRUDOS (.MAT) ---
fprintf('Guardando datos crudos en .mat...\n');
% Eje Frecuencia recortado usado para visualización
mask_vis_final = abs(f_vec) <= anchoBanda/2;
Freq_Save = f_abs(mask_vis_final);

nombreMat = fullfile(dir_base_dat, ['Datos_Procesados_Indirecto_' timestamp '.mat']);
try
    save(nombreMat, 'Accum_P', 'Accum_Time', 'Freq_Save', 'lista_final_eventos', 'archivoSeleccionado', '-v7.3');
    fprintf('Datos guardados exitosamente en: %s\n', nombreMat);
catch ME
    fprintf('Error al guardar .mat: %s\n', ME.message);
end
if ~isempty(lista_final_eventos)
    fprintf('\n--- RESULTADOS ESTADÍSTICOS (Indirecto) ---\n');
    fprintf('Eventos Totales: %d\n', length(lista_final_eventos));

    if ~isempty(lista_final_eventos)
        fprintf('\n>>> MÉTRICAS ROBUSTAS (Promedio General) <<<\n');
        fprintf('  - Diferencia Pico-Ruido: %.2f dBm\n', mean([lista_final_eventos.Diferencia_Potencia_dBm]));
        fprintf('  - Potencia Máxima:       %.2f dBm\n', mean([lista_final_eventos.P_max_dBm]));
        fprintf('  - Delta Estabilidad:     %.2f dBm\n', mean([lista_final_eventos.Delta_Estabilidad_dBm]));
    end

    [~, idx] = sort([lista_final_eventos.Diferencia_Potencia_dBm], 'descend');
    lista_final_eventos = lista_final_eventos(idx);

    % Generar Excel antes de graficar
    nombreExcel = fullfile(dir_base_dat, ['Reporte_Metricas_Indirecto_' datestr(now, 'yyyymmdd_HHMMSS') '.xlsx']);
    registrarMetrica('Indirecto', lista_final_eventos, nombreExcel);
end

% --- GRAFICADO FINAL ---
if ~run_batch_mode
    fprintf('Generando gráficos finales...\n');
    % Eje Frecuencia recortado
    mask_vis = abs(f_vec) <= anchoBanda/2;
    Freq_Plot = f_abs(mask_vis)/1e6;

    % Downsampling visual si es muy grande (opcional, por ahora directo)
    % surf(X, Y, Z) -> X=Frec, Y=Time, Z=Potencia(Frec,Time)'

    surf(ax2D, Freq_Plot, Accum_Time, Accum_P.', 'EdgeColor','none');
    surf(ax3D, Freq_Plot, Accum_Time, Accum_P.', 'EdgeColor','none');

    % Marcadores (dependen de lista_final_eventos ordenada)
    if ~isempty(lista_final_eventos)
        for p=1:length(lista_final_eventos)
            pk=lista_final_eventos(p);
            z_mark = min(max(pk.P_max_dBm, rangoColores(1)), rangoColores(2));
            if 1, c='m'; else, c='c'; end
            plot3(ax3D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, z_mark, 'v', 'MarkerFaceColor',c,'MarkerEdgeColor','k');
            plot3(ax2D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, 200, 'v', 'MarkerFaceColor',c,'MarkerEdgeColor','k');
            if p <= 5
                text(ax3D, pk.Freq_Hz/1e6, pk.Tiempo_Max_Seg, z_mark, sprintf(' Dif: %.1f', pk.Diferencia_Potencia_dBm), 'Color',c, 'BackgroundColor', 'w', 'Margin', 1);
            end
        end
    end

    drawnow;
end

%% GUARDADO DE IMAGENES
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

if ~run_batch_mode
    saveas(fig2D, fullfile(dir_base_img, ['INDIRECTO_2D_SinMarcadores_' name_no_ext '.png']));
    saveas(fig2D, fullfile(dir_base_img, ['INDIRECTO_2D_SinMarcadores_' name_no_ext '.fig']));
    saveas(fig3D, fullfile(dir_base_img, ['INDIRECTO_3D_SinMarcadores_' name_no_ext '.png']));
    saveas(fig3D, fullfile(dir_base_img, ['INDIRECTO_3D_SinMarcadores_' name_no_ext '.fig']));
end


return;

if ~run_batch_mode
    saveas(fig2D, fullfile(dir_base_img, ['INDIRECTO_2D_ConMarcadores_' name_no_ext '.png']));
    saveas(fig2D, fullfile(dir_base_img, ['INDIRECTO_2D_ConMarcadores_' name_no_ext '.fig']));
    saveas(fig3D, fullfile(dir_base_img, ['INDIRECTO_3D_ConMarcadores_' name_no_ext '.png']));
    saveas(fig3D, fullfile(dir_base_img, ['INDIRECTO_3D_ConMarcadores_' name_no_ext '.fig']));
end

spec_profile = sum_spec / max(count_spec, 1);
if ~run_batch_mode
    figProfile = figure('Name', 'Perfil Promedio', 'Color', 'w');
    plot(f_abs/1e6, spec_profile, 'k', 'LineWidth', 1.2);
    grid on; xlabel('Frecuencia (MHz)'); ylabel('Amplitud Promedio (dBm)');
    xline(f_HI/1e6, 'r--', 'LineWidth', 1);
    title('Estimación Espectral (Indirecto)');
    saveas(figProfile, fullfile(dir_base_img, ['INDIRECTO_Perfil_' name_no_ext '.png']));
    saveas(figProfile, fullfile(dir_base_img, ['INDIRECTO_Perfil_' name_no_ext '.fig']));
end

disp('Proceso Indirecto completado.');

if run_batch_mode
    if ~exist('spec_profile', 'var')
        spec_profile = sum_spec / max(count_spec, 1);
    end
    results_batch.Indirecto.f = f_abs;
    results_batch.Indirecto.P = spec_profile;
end
