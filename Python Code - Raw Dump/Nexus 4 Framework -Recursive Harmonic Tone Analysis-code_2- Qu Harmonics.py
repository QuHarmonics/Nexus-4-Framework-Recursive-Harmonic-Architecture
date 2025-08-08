import numpy as np
import scipy.io.wavfile as wavfile

# Parameters
sample_rate = 44100  # Hz
duration_per_tone = 0.5  # seconds
phi = 0.61803398875
H = 0.35
F = 1.0
base_freq = 80.0  # base frequency in Hz

# Generate tones based on R values
tones = []
for i in range(100):
    t = phi * i
    R = np.exp(H * F * t)
    freq = base_freq * R
    time_array = np.linspace(0, duration_per_tone, int(sample_rate * duration_per_tone), endpoint=False)
    tone_wave = 0.5 * np.sin(2 * np.pi * freq * time_array)
    tones.append(tone_wave)

# Concatenate all tones
audio_signal = np.concatenate(tones)

# Normalize to 16-bit PCM range
audio_signal_int16 = np.int16(audio_signal / np.max(np.abs(audio_signal)) * 32767)

# Save to WAV file
wav_path = "d:\\harmonic_recursive_tones.wav"
wavfile.write(wav_path, sample_rate, audio_signal_int16)

wav_path
