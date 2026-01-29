function metrics = calcularMetricas(f, P_dBm, t, config)
% calcularMetricas - Calcula métricas robustas para detección de HI.
%
% Entradas:
%   f      : Vector de frecuencias (Hz)
%   P_dBm   : Matriz o Vector de Potencia (dBm). Si es matriz: [Frec x Tiempo]
%   t      : Vector de tiempo (s)
%   config : Estructura con parámetros (opcional, por compatibilidad)
%
% Salida:
%   metrics : Vector de estructuras con los eventos detectados.
%             Campos: Freq_Hz, Tiempo_Max_Seg, P_max_dBm, P_ruido_dBm,
%                     Diferencia_Potencia_dBm, Delta_Estabilidad_dBm,
%                     Validacion_Estabilidad (0/1).

metrics = [];

% 0. Validaciones básicas
if isempty(P_dBm) || isempty(f), return; end

% Asegurar dimensiones (Power debe ser vector para picos instantáneos o matriz)
% Si es matriz, operamos sobre el perfil "Max Hold" o promedio según lógica.
% Para detección robusta, usaremos el Máximo en el bloque.

if ismatrix(P_dBm) && size(P_dBm, 2) > 1
    % Si entra un bloque (espectrograma), colapsamos a perfil máximo para detección
    [P_perfil, idx_t_max] = max(P_dBm, [], 2); % Max sobre tiempo
else
    P_perfil = P_dBm(:);
    idx_t_max = ones(size(P_perfil));
end

% 1. Encontrar Pico Máximo y Piso de Ruido Global
[P_max_val, idx_peak] = max(P_perfil);
freq_peak = f(idx_peak);

% Estimación robusta del ruido (Mediana de todo el espectro)
% Excluimos la zona inmediata del pico para no contaminar el ruido
ancho_excl = 20e3; % +/- 20 kHz alrededor del pico
mask_noise = abs(f - freq_peak) > ancho_excl;

if any(mask_noise)
    P_ruido_val = median(P_perfil(mask_noise));
else
    P_ruido_val = median(P_perfil);
end

% 2. Métrica Principal: Diferencia Pico - Ruido
Diferencia_dBm = P_max_val - P_ruido_val;

% 3. Análisis de Estabilidad (Solo si es matriz y tenemos dimensión tiempo)
Delta_Estabilidad = 0;
% Es_Estable = 1; % Eliminado

tiempo_evento = t(1); % Default

if ismatrix(P_dBm) && size(P_dBm, 2) > 1
    % Recuperamos la traza temporal en la frecuencia del pico
    traza_temporal = P_dBm(idx_peak, :);

    % Delta Estabilidad: Variación (Max - Min) en esa frecuencia durante el bloque
    % Un transitorio tendrá un Delta muy alto. Una señal continua, bajo.
    Delta_Estabilidad = max(traza_temporal) - min(traza_temporal);

    % Criterio de Validación: Delta < 6 dBm (definido por usuario)
    % (Logica de transietoriedad eliminada por solicitud del usuario)
    % if Delta_Estabilidad < 6
    %     Es_Estable = 1;
    % else
    %     Es_Estable = 0; % Transitorio
    % end

    % Tiempo exacto del máximo
    tiempo_evento = t(idx_t_max(idx_peak));
end

% 4. Empaquetar Resultado
% Solo reportamos el pico dominante del bloque/segmento para robustez.

m_struct = struct();
m_struct.Freq_Hz = freq_peak;
m_struct.Tiempo_Max_Seg = tiempo_evento;

m_struct.P_max_dBm = P_max_val;
m_struct.P_ruido_dBm = P_ruido_val;
m_struct.P_promedio_dBm = mean(P_perfil); % Métrica auxiliar

m_struct.Diferencia_Potencia_dBm = Diferencia_dBm; % A MAXIMIZAR
m_struct.Delta_Estabilidad_dBm = Delta_Estabilidad; % A MINIMIZAR (<6)
% m_struct.Validacion_Estabilidad = Es_Estable; % Eliminado

metrics = [metrics, m_struct];
end
