import numpy as np
import matplotlib.pyplot as plt
from FIRWIN import FIRWIN, PI

# %% Window TEST
N = 100

w1 = FIRWIN.Window.Rectangular(N)
w2 = FIRWIN.Window.Hanning(N)
w3 = FIRWIN.Window.Hamming(N)
w4 = FIRWIN.Window.Blackman(N)
w5 = FIRWIN.Window.Kaiser(N, 5)

plt.figure(figsize=(12.8, 7.2))
def Draw_Win(w, N:int, window_name):
    plt.plot(w, linewidth = 5, label=f"{window_name}")
plt.title("Window Function")
Draw_Win(w1, N, "Rectangular")
Draw_Win(w2, N, "Hanning")
Draw_Win(w3, N, "Hamming")
Draw_Win(w4, N, "Blackman")
Draw_Win(w5, N, "Kaiser")
plt.legend()
plt.tight_layout()

# %% Design Examples
Rp = 0.1 # ripple [dB]
As = 50  # decay (passband v.s. stopband)

# LPF & HPF
Fs = 8000 # sampling frequency
fs = 2000 # stopband frequency
fp = 1500 # passband frequency

delta_f = fs - fp

LPF = FIRWIN(fp, delta_f, Rp, As, Fs, "Lowpass", draw=True)
h = LPF.design()
# Freq_response(h)

HPF = FIRWIN(fp, delta_f, Rp, As, Fs, "Highpass", draw=True)
h = HPF.design()
# Freq_response(h)

# BPF & BRJF
ws1 = 0.1*PI * 8e3  # stopband1 angular frequency
wp1 = 0.25*PI * 8e3 # passband1 angular frequency
ws2 = 0.7*PI * 8e3  # stopband2 angular frequency
wp2 = 0.55*PI * 8e3 # passband2 angular frequency

fp1 = wp1 / (2*PI)
fp2 = wp2 / (2*PI)

delta_w1 = wp1 - ws1
delta_w2 = ws2 - wp2
delta_w = min(delta_w1, delta_w2) # Choose the minimum trasition band
delta_f = delta_w / (2*PI)

BPF = FIRWIN([fp1, fp2], delta_f, Rp, As, Fs, type="Bandpass", draw=True)
h = BPF.design()
# Freq_response(h)

BRJF = FIRWIN([fp1, fp2], delta_f, Rp, As, Fs, type="Bandreject", draw=True)
h = BRJF.design()
# Freq_response(h)

plt.show()