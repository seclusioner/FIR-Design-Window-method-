from scipy.fftpack import fft,ifft
import numpy as np
import matplotlib.pyplot as plt

PI = np.pi

class Window():
    def __init__(self):
        pass

    def Rectangular(N):
        return [1] * N

    def Hanning(N):
        return [.5 + .5 * np.cos(2 * PI * n / (N - 1)) for n in range(N)]

    def Hamming(N):
        return [.54 - .46 * np.cos(2*np.pi*n / (N-1)) for n in range(N)]

    def Blackman(N):
        return [.42 - .5 * np.cos(2*np.pi*n / (N-1)) + .08 * np.cos(4*np.pi*n / (N-1)) for n in range(N)]

# %% FIR window LPF
Fs = 8000
fs = 2000 / Fs
fp = 1500 / Fs
delta_f = fs - fp
fp += delta_f / 2
wp = 2 * np.pi * fp

k = 3.3
N = int(np.ceil(k / delta_f))

hd=np.zeros(N)
w=Window.Hamming(N)
for n in range(N):
    if n==(N-1)/2:
        hd[n] = 2*fp
    else:
        hd[n] = 2*fp*np.sin((n-(N-1)/2)*wp) / ((n-(N-1)/2)*wp)
    
h = hd * w
plt.figure(figsize=(12.8, 7.2))
plt.subplot(2,1,1)
plt.stem(h)
H = fft(h, 8192)
plt.subplot(2,1,2)
plt.plot(20 * np.log10(np.abs(H[0:int(len(H) / 2)])))
plt.tight_layout()

# %% FIR window HPF
Fs = 10
fs = 2.5 / Fs
fp = 2 / Fs
delta_f = fs - fp
fp += delta_f / 2

k = 3.3
N = int(np.ceil(k / delta_f))

wp = 2 * np.pi * fp

hd=np.zeros(N)
w=np.zeros(N)
for n in range(N):
    w[n] = 0.54 - 0.46*np.cos(2*np.pi*n / (N-1))
    if n==(N-1)/2:
        hd[n] = 1-2*fp
    else:
        hd[n] = -2*fp*np.sin((n-(N-1)/2)*wp) / ((n-(N-1)/2)*wp)
    
h = hd * w

plt.figure(figsize=(12.8, 7.2))
plt.subplot(2,1,1)
plt.stem(h)
H = fft(h, 8192)
plt.subplot(2,1,2)
plt.plot(20 * np.log10(np.abs(H[0:int(len(H) / 2)])))
plt.tight_layout()

# %% FIR window BPF
ws1 = 0.1*PI
wp1 = 0.25*PI
ws2 = 0.7*PI
wp2 = 0.55*PI

delta_w1 = wp1 - ws1
delta_w2 = ws2 - wp2
delta_w = min(delta_w1, delta_w2)
delta_f = delta_w / (2*PI)

N = np.ceil(5.5 / delta_f).astype(np.int) + 1
# print(n)

fp1 = wp1 / (2*PI)
fp2 = wp2 / (2*PI)

hd = np.zeros(N)
w  = np.zeros(N)
for n in range(N):
    w[n] = 0.42 - 0.5 * np.cos(2 * PI * n / (N-1)) + 0.08 * np.cos(4*PI*n / (N-1))
    if n == (N-1)/2:
        hd[n] = 2 * (fp2 - fp1)
    else:
        hd[n] = 2 * fp2 * np.sin(wp2 * (n - (N-1) / 2)) / ((n - (N-1) / 2) * wp2) - 2 * fp1 * np.sin(wp1 * (n - (N-1) / 2)) / ((n - (N-1) / 2) * wp1)

h = hd * w

plt.figure(figsize=(12.8, 7.2))
plt.subplot(2,1,1)
plt.stem(h)
H = fft(h, 8192)
plt.subplot(2,1,2)
plt.plot(20 * np.log10(np.abs(H[0:int(len(H) / 2)])))

# %% FIR window BRJF
ws1 = 0.1*PI
wp1 = 0.25*PI
ws2 = 0.7*PI
wp2 = 0.55*PI

delta_w1 = wp1 - ws1
delta_w2 = ws2 - wp2
delta_w = min(delta_w1, delta_w2)
delta_f = delta_w / (2*PI)

N = np.ceil(5.5 / delta_f).astype(np.int) + 1
# print(n)

fp1 = wp1 / (2*PI)
fp2 = wp2 / (2*PI)

hd = np.zeros(N)
w  = np.zeros(N)
for n in range(N):
    w[n] = 0.42 - 0.5 * np.cos(2 * PI * n / (N-1)) + 0.08 * np.cos(4*PI*n / (N-1))
    if n==(N-1)/2:
        hd[n] = 1 - 2*(fp2-fp1)
    else:
        hd[n] = -(2 * fp2 * np.sin(wp2 * (n - (N-1) / 2)) / ((n - (N-1) / 2) * wp2) - 2 * fp1 * np.sin(wp1 * (n - (N-1) / 2)) / ((n - (N-1) / 2) * wp1))
    
h = hd * w
plt.figure(figsize=(12.8, 7.2))
plt.subplot(2,1,1)
plt.stem(h)
H = fft(h, 8192)
plt.subplot(2,1,2)
plt.plot(20 * np.log10(np.abs(H[0:int(len(H) / 2)])))
plt.tight_layout()

plt.show()