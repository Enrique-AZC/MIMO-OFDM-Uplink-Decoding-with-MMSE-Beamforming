# MIMO-OFDM Uplink Decoding with MMSE Beamforming

## Overview  
This project implements a **MIMO-OFDM uplink decoding system** with **Minimum Mean Square Error (MMSE) beamforming** at the base station (BS).  
The system simulates a **10-antenna BS** receiving uplink transmissions from **3 mobile users**, each transmitting QAM-16 modulated signals.  

The full pipeline includes **pilot-based channel estimation**, **MMSE equalization**, **symbol demodulation**, and **audio reconstruction** from the decoded streams. The results show **excellent separation of users, accurate channel estimation, and faithful audio recovery** in both PCM8 and PCM16 formats.

---

## Authors  
- Enrique Alcalá-Zamora Castro  
- Andrea Bañón Triguero   

---

## Objectives  

- Decode uplink transmissions in a **MIMO-OFDM system** using MMSE beamforming.  
- Estimate frequency-selective channels using **orthogonal pilot sequences**.  
- Apply **per-subcarrier MMSE equalization** to separate multiple users.  
- Recover **QAM-16 symbols** and reconstruct original audio signals.  
- Validate performance through **bit/symbol error rate (BER/SER)** and audio quality.  

---

## Methodology  

### System Parameters  
- **Base Station (BS):** \( M = 10 \) antennas  
- **Users (MS):** \( K = 3 \)  
- **Subcarriers:** \( N_{FFT} = 256 \)  
- **Cyclic Prefix (CP):** \( L_{cp} = 4 \)  
- **Block length:** \( T_c = 10 \) OFDM symbols per coherence block  
  - Pilots: \( T_p = 4 \)  
  - Data: \( T_u = T_c - T_p = 6 \)  
- **Modulation:** QAM-16  
- **SNR:** 8 dB  

### Processing Pipeline  

1. **Pilot Design**  
   - Orthogonal pilots generated via **DFT matrix**: `fft(eye(Tp))`  
   - First \(K\) rows assigned to users.  

2. **Signal Reception**  
   - Matrix `xmimo` stores BS antenna signals.  
   - Partitioned into **coherence blocks**.  

3. **OFDM Demodulation**  
   - Remove CP, apply FFT per antenna.  
   - Separate **pilot symbols** and **data symbols**.  

4. **Channel Estimation**  
   - For each subcarrier:  
     \[
     \hat{H}(f) = Y_p \cdot P^\dagger
     \]  
     where \(Y_p\) are received pilots, \(P\) is pilot matrix.  

5. **MMSE Beamforming**  
   - For each subcarrier:  
     \[
     W_f = (H_f^H H_f + \frac{1}{SNR} I)^{-1} H_f^H
     \]  
   - Apply to received data:  
     \[
     \hat{s} = W_f \cdot X_{dat}
     \]  

6. **QAM-16 Demodulation**  
   - Convert complex symbols to integers \([0..15]\).  
   - Store in `simbolos_output.mat` for performance evaluation.  

7. **Audio Reconstruction**  
   - Pack QAM symbols into **bytes** (2 symbols per byte).  
   - Export per-user audio as:  
     - **PCM8 (unsigned, 8-bit)**  
     - **PCM16 (signed, 16-bit, rescaled)**  

---

## Results  

- **Channel Estimation:** Accurate per-subcarrier LS estimation.  
- **Beamforming:** MMSE successfully separates users despite interference.  
- **Demodulation:** Robust QAM-16 decoding with very low SER at SNR = 8 dB.  
- **Audio Quality:**  
  - Clean reconstruction of user audio streams.  
  - No audible distortion after PCM8 and PCM16 reconstruction.  
  - Clamping ensures valid dynamic range in PCM16.  

The system achieves **excellent overall performance**: multiuser separation, reliable symbol detection, and high-quality audio recovery.

---

## Technical Relevance  

This project demonstrates:  
- **Advanced multi-antenna signal processing** for MIMO uplink.  
- **OFDM baseband processing pipeline** with pilot-based channel estimation.  
- **Optimal linear receivers (MMSE)** in multiuser environments.  
- **End-to-end validation** from digital communications theory to perceptual audio recovery.  

---

## Future Work  

- Implement alternative beamformers: **Zero-Forcing (ZF)**, **Maximum Ratio Combining (MRC)**.  
- Evaluate under **time-varying channels** with Doppler.  
- Extend to **higher-order modulations** (e.g., QAM-64, QAM-256).  
- Perform **BER/SER statistical analysis** vs. SNR.  
- Simulate **massive MIMO (M > 100 antennas)** scenarios.  
- Incorporate **FEC coding (e.g., LDPC, Turbo)** for error resilience.  

---

## Requirements  

- **MATLAB R2021a+** or **GNU Octave 6.0+**  
- Input file: `signals_input.mat` containing `xmimo` (received signals).  
- Auxiliary functions: `demod_ofdm_block.m`, `qamdemod`.  

---

## License  

This work is released for academic and research purposes. Proper attribution is required.  
