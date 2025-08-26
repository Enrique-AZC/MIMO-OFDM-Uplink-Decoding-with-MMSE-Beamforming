%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Proyecto: Decodificación MIMO-OFDM Uplink con Beamforming MMSE
%
% Grupo: Enrique Alcalá-Zamora Castro
%        Andrea Bañón Triguero
%        José Palenzuela Porcel
%        Heiner Fernando Buitrago
%        Silvia Daniela Arechúa Batalla
%        
%
% Descripción general:
% Este script realiza la decodificación de señales MIMO-OFDM en el uplink
% utilizando beamforming MMSE. A partir de las señales recibidas en la BS
% (10 antenas), se estiman los canales, se aplica ecualización espacial 
% y se reconstruyen los símbolos QAM transmitidos por 3 usuarios. 
% Finalmente, se reconstruyen las señales de audio originales en formato 
% WAV de 8 y 16 bits.

% Estructura del procesamiento:
%
% 1. Definición de parámetros del sistema
%
% 2. Construcción de las secuencias piloto ortogonales:
%    - Se usa la DFT (fft(eye)) para generar pilotos ortogonales
%    - Se seleccionan las primeras K filas para los usuarios

% 3. Carga de señales recibidas en la BS:
%    - xmimo: matriz M×Nt, donde cada fila es la señal recibida por una antena
%    - Se dividen en nB bloques de coherencia

% 4. Demodulación OFDM y procesamiento por bloque:
%    - Se elimina el CP y se aplica FFT para obtener Xcoh
%    - Se separan símbolos piloto (Xpil) y símbolos de datos (Xdat)
%    - Se estima el canal H(f) por subportadora usando LS: H = Yp * pinv(P)
%    - Se aplica beamforming MMSE: W = (H'H + (1/SNR)*I)^(-1) * H'
%    - Se decodifican los datos: ŝ = W * Xdat

% 5. Decisión QAM-16:
%    - Se demodulan los símbolos complejos a enteros (0..15)
%    - Se guarda en 'simbolos_output.mat' como matriz W4rec (K × N)

% 6. Reconstrucción de audio:
%    - Cada 2 símbolos QAM (4 bits cada uno) se combinan en un byte (8 bits)
%    - Se generan archivos WAV en:
%       • Formato PCM8 (8-bit unsigned)
%       • Formato PCM16 (16-bit signed con mapeo lineal)
%    - Se realiza clamping para evitar saturación de audio

% Resultado:
% - El sistema logra separar correctamente las señales de los 3 usuarios,
%   estimar el canal, aplicar beamforming MMSE y reconstruir el audio.
% - Se obtienen archivos: simbolos_output.mat, user1/2/3_PCM8.wav, user1/2/3_PCM16.wav

% NOTA:
% Aunque no se ha usado la función 'qam16_to_audio8bit', su funcionalidad
% ha sido implementada directamente mediante ensamblado de nibbles (LSB, MSB).
% Esto garantiza una reconstrucción fiel del audio sin distorsiones.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc; close all;

%% 1. Parámetros del sistema MIMO-OFDM
M       = 10;                   % Número de antenas en la estación base (BS)
K       = 3;                    % Número de usuarios móviles (MS)
Nfft    = 256;                  % Número de subportadoras OFDM
Lcp     = 4;                    % Longitud del prefijo cíclico (CP)
Tc      = 10;                   % Número total de símbolos OFDM por bloque de coherencia
Tp      = 4;                    % Número de símbolos piloto por bloque
Tu      = Tc - Tp;              % Número de símbolos de datos por bloque
Es      = 10;                   % Energía promedio por símbolo (QAM-16)
SNR_dB  = 8;                    % SNR en dB
SNR_lin = 10^(SNR_dB/10);       % SNR en escala lineal

% Generación de pilotos ortogonales mediante DFT
pilots = sqrt(Es)*fft(eye(Tp)); % Tp × Tp matriz de pilotos
Pmat   = pilots(1:K,:);         % Se seleccionan las primeras K filas → K × Tp

%% 2. Carga de señales y preasignación
load('signals_input.mat','xmimo');  % Carga señales M×Nt (M=10, Nt=nº muestras)
Nt      = size(xmimo,2);            % Número total de muestras temporales
Tcoh    = Tc*(Nfft+Lcp);            % Muestras por bloque de coherencia
nB      = floor(Nt/Tcoh);           % Número de bloques de coherencia completos

% Reservas de matrices para almacenamiento
H_est    = zeros(M,K,Nfft,nB);      % Canal estimado H para cada bloque
symb_est = zeros(K,Nfft,Tu,nB);     % Símbolos estimados tras beamforming

%% 3. Procesamiento por bloque de coherencia
for b = 1:nB
  % 3.1: Extrae el bloque de coherencia b
  i0   = (b-1)*Tcoh + 1;
  xcoh = xmimo(:, i0:i0+Tcoh-1);    % M × (Nfft+Lcp) × Tc

  % 3.2: Demodulación OFDM (quita CP y aplica FFT)
  Xcoh = demod_ofdm_block(xcoh,Nfft,Lcp); % M × (Nfft×Tc)
  X3   = reshape(Xcoh, M, Nfft, Tc);      % M × Nfft × Tc

  % 3.3: Separación de pilotos y datos
  Xpil = X3(:,:,1:Tp);      % M × Nfft × Tp (símbolos piloto)
  Xdat = X3(:,:,Tp+1:end);  % M × Nfft × Tu (símbolos de datos)

  % 3.4: Estimación de canal H(f) por subportadora
  for f = 1:Nfft
    Yp = squeeze(Xpil(:,f,:));        % M × Tp
    H_est(:,:,f,b) = Yp * pinv(Pmat); % M × K (estimación LS)
  end

  % 3.5: Beamforming MMSE para cada subportadora
  for f = 1:Nfft
    Hf = squeeze(H_est(:,:,f,b));       
    Wf = (Hf'*Hf + (1/SNR_lin)*eye(K)) \ Hf';  % K × M (MMSE weights)
    for t = 1:Tu
      symb_est(:,f,t,b) = Wf * Xdat(:,f,t);    % Decodificación: ŝ = W·x
    end
  end
end

%% 4. Decisión QAM-16 y guardado
SY    = reshape(symb_est, K, []);  % K × Nsímbolos (vector por usuario)
W4rec = qamdemod(SY, 16, 'UnitAveragePower', false); % Decisión QAM-16
save('simbolos_output.mat','W4rec');  % Guarda resultados para evaluación SER

%% 5. Reconstrucción de audio para cada usuario
fs = 8000; % Frecuencia de muestreo del audio (Hz)

for k = 1:K
  cwk = W4rec(k,:);  % Codewords QAM-16 del usuario k

  % Verifica paridad (número par de símbolos)
  if mod(numel(cwk),2)
    cwk(end+1) = 0;  % Añade un cero si es impar
  end

  cw_pairs = reshape(cwk,2,[]).';  % n × 2: [LSB MSB] por fila

  % Combina los nibbles en bytes de 8 bits
  LSB = uint8(cw_pairs(:,1));
  MSB = uint8(cw_pairs(:,2));
  audio8 = bitshift(MSB,4) + LSB;  % Combina como: byte = (MSB<<4) | LSB

  % ===== OPCIÓN 1: Exportar WAV en 8-bit unsigned =====
  fn8 = sprintf('user%d_PCM8.wav',k);
  audiowrite(fn8, audio8, fs, 'BitsPerSample', 8);
  fprintf('Escrito %s\n', fn8);

  % ===== OPCIÓN 2: Exportar WAV en 16-bit signed =====
  % Reescala de [0,255] a [-32768,32767]
  audio16 = int16(double(audio8) * (65535/255) - 32768);
  % Clamping para evitar valores fuera del rango válido
  audio16 = max(min(audio16,32767), -32768);

  fn16 = sprintf('user%d_PCM16.wav',k);
  audiowrite(fn16, audio16, fs, 'BitsPerSample', 16);
  fprintf('Escrito %s\n', fn16);
end