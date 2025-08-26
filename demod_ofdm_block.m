function Xcoh=demod_ofdm_block(xcoh,M,L)

% Inputs:
%   xcoh: signal recibida (una fila por antena)
%   M: no. de subportadoras
%   L: intervalo de guarda
%   Tcoh: duracion (en no. de simbolos OFDM) del intervalo de coherencia

% Output Xcoh: bloque COH de simbolos en tiempo-frecuencia

N=M+L;      %Duracion del simbolo OFDM con intervalo de guarda
[Narr,aux]=size(xcoh); %No. antenas receptoras, Tcoh*N
Tcoh=aux/N; %Duracion (en no. de simbolos OFDM) del intervalo COH
Xcoh=zeros(Narr,M*Tcoh);
for nsimb=0:Tcoh-1 
  % Demodular OFDM
  xsimb=xcoh(:,nsimb*N+L+1:(nsimb+1)*N); %Quitar intervalos de guarda
  Xsimb=fft(xsimb.').'; %Pasar simbolo a dominio freq
  Xcoh(:,nsimb*M+1:(nsimb+1)*M)=Xsimb;
end