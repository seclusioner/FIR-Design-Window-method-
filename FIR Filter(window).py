import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft

PI = np.pi

def Freq_response(h):
    plt.figure(figsize=(12.8, 7.2))
    plt.subplot(2,1,1)
    plt.title("Impulse Response")
    plt.stem(h)
    H = fft(h, 8192)
    plt.subplot(2,1,2)
    plt.title("Frequency Response (log scale)")
    plt.plot(20 * np.log10(np.abs(H[0:int(len(H) / 2)])), linewidth = 5)
    plt.tight_layout()
        
class FIRWIN():
    class Window():
        def __init__(self):
            pass

        def Rectangular(N):
            return [1] * N

        def Hanning(N):
            return [.5 - .5 * np.cos(2 * PI * n / (N - 1)) for n in range(N)]

        def Hamming(N):
            return [.54 - .46 * np.cos(2*np.pi*n / (N-1)) for n in range(N)]

        def Blackman(N):
            return [.42 - .5 * np.cos(2*np.pi*n / (N-1)) + .08 * np.cos(4*np.pi*n / (N-1)) for n in range(N)]
        
        def Kaiser(N, beta):
            from scipy.special import i0
            return [i0(beta * np.sqrt(1 - ((2.0 * n) / (N - 1) - 1)**2)) / i0(beta) for n in range(N)]

    def __init__(self, fp, df, Rp, As, Fs, type="Lowpass", draw=False):
        if isinstance(fp, int):
            self.fp   = fp / Fs
        elif isinstance(fp, list):
            self.fp1  = fp[0] / Fs
            self.fp2  = fp[1] / Fs
        self.df   = df / Fs
        self.Rp   = Rp
        self.As   = As
        self.Fs   = Fs
        self.type = type
        self.draw = draw
        
    def design(self):
        '''
        Window          Transition      Passband        Main/Side Lobe      Stopband
        Function        Width (Hz)      Ripple (dB)     Relation (dB)       Attenuation (dB)
        -------------------------------------------------------------------------------------
        Rectangular     0.9 / N         0.7416          13                  21
        
        Hanning         3.1 / N         0.0546          31                  44

        Hamming         3.3 / N         0.0194          41                  53

        Blackman        5.5 / N         0.0017          57                  74

               ~        2.93 / N        0.0274          -- (beta = 4.54)    50
        Kaiser ~        4.32 / N        0.00275         -- (beta = 6.76)    70
               ~        5.71 / N        0.000275        -- (beta = 8.96)    90
        '''
        
        # Choose Window
        if self.As < 20: # Rectangular
            n = self.getN(0.9) 
            n = n if n % 2 == 1 else self.getN(0.9) + 1
            w = self.Window.Rectangular(n)
            
        elif self.As < 43: # Hanning
            n = self.getN(3.1)
            n = n if n % 2 == 1 else self.getN(3.1) + 1
            w = self.Window.Hanning(n)
            
        elif self.As < 53: # Hamming
            n = self.getN(3.3)
            n = n if n % 2 == 1 else self.getN(3.3) + 1
            w = self.Window.Hamming(n)
            
        elif self.As < 73: # Blackman
            n = self.getN(5.5)
            n = n if n % 2 == 1 else self.getN(5.5) + 1
            w = self.Window.Blackman(n)
        else:
            n = self.getN(5.71)
            n = n if n % 2 == 1 else self.getN(5.71) + 1
            w = self.Window.Kaiser(n, 8.96)
                
        '''
        type  of  filter               ideal  impulse  response
          hd[n], n!=0                           hd[0]
        ------------------------------------------------------------
         Lowpass     2fc*sin[nwc]/nwc           2fc
        
         Highpass    -2fc*sin[nwc]/nwc          1-2fc
        
         Bandpass    2f2*sin[nwc2]/nwc2         2(f2-f1)
                     - 2f1*sin[nwc1]/nwc1
        
         Bandreject  2f1*sin[nwc1]/nwc1         1-2(f2-f1)
                     - 2f2*sin[nwc2]/nwc2
        ------------------------------------------------------------
        '''
        
        if self.type == "Lowpass":
            hd = self.hd_lpf(n)
            
        elif self.type == "Highpass":
            hd = self.hd_hpf(n)
            
        elif self.type == "Bandpass":
            hd = self.hd_bpf(n)
            
        elif self.type == "Bandreject":
            hd = self.hd_brjf(n)
        else:
            print("Input is an unsupported type")
            return 0
        
        h = hd * np.array(w)
        
        if(self.draw):
            self.Filter_spec(h, Fs)
            
        return h

    def getN(self, k): # N = k / df (order of filter)
        return math.ceil(k / self.df)
  
    def hd_lpf(self, N):
        fp = self.fp + self.df / 2
        wp = 2 * PI * fp        
        return [(2 * fp) if n == (N-1)//2 else (2 * fp * np.sin((n - (N-1) // 2) * wp) / ((n - (N-1) // 2) * wp)) for n in range(N)]

    def hd_hpf(self, N):
        fp = self.fp + self.df / 2
        wp = 2 * PI * fp        
        return [(1 - 2 * fp) if n == (N-1)//2 else -(2 * fp * np.sin((n - (N-1) // 2) * wp) / ((n - (N-1) // 2) * wp)) for n in range(N)]

    def hd_bpf(self, N):
        wp1 = 2 * PI * self.fp1
        wp2 = 2 * PI * self.fp2
        return [(2 * (self.fp2 - self.fp1)) if n == (N-1)//2 else (2 * self.fp2 * np.sin(wp2 * (n - (N-1) / 2)) / ((n - (N-1) / 2) * wp2) - 2 * self.fp1 * np.sin(wp1 * (n - (N-1) / 2)) / ((n - (N-1) / 2) * wp1)) for n in range(N)]

    def hd_brjf(self, N):
        wp1 = 2 * PI * self.fp1
        wp2 = 2 * PI * self.fp2
        return [(1 - 2 * (self.fp2 - self.fp1)) if n == (N-1)//2 else (-(2 * self.fp2 * np.sin(wp2 * (n - (N-1) / 2)) / ((n - (N-1) / 2) * wp2) - 2 * self.fp1 * np.sin(wp1 * (n - (N-1) / 2)) / ((n - (N-1) / 2) * wp1))) for n in range(N)]

    def Filter_spec(self, h: list, sampling_rate: int):
        J= 0 + 1j
        def DTFT(x: list, N: int):
            Xjw = []
            w = np.linspace(0, np.pi, N)
            for i in range(N):
                ans = 0
                for n in range(len(x)):
                    ans += x[n] * np.exp(-J*w[i]*n)
                Xjw.append(ans)

            return Xjw
    
        N=len(h)
        n=np.arange(len(h))
        plt.figure(figsize=(12, 7.2))
        plt.suptitle('Filter Response')
        plt.subplot(321)
        plt.stem(n, h)
        plt.title('impulse response')
        plt.tight_layout()

        DATA=5000
        Xjw=DTFT(h, DATA)
        n=np.arange(len(Xjw))
        delta_f = (n/DATA)*(sampling_rate//2)

        plt.subplot(322)
        plt.plot(delta_f, np.abs(Xjw))
        plt.title('frequency response(DTFT)')
        plt.ylabel('magnitude |Xj\u03C9|') # unicode
        plt.xlabel('f(Hz)')
        plt.tight_layout()
        plt.subplot(323)
        plt.title('frequency response(DTFT)')
        plt.ylabel('magnitude(dB)')
        plt.xlabel('f(Hz)')
        plt.plot(delta_f, 20*np.log(np.abs(Xjw)))
        plt.tight_layout()

        # phase response --------------------------------
        plt.subplot(324)
        plt.title("Phase Response")
        plt.phase_spectrum(h)
        plt.tight_layout()

        X=fft(h)
        N=len(X)
        n=np.arange(N)
        delta_f = n*sampling_rate/N
        plt.subplot(325)
        plt.title('frequency response(FFT)')
        plt.ylabel('magnitude |X(\u03C9)|')
        plt.xlabel('f(Hz)')
        plt.xlim(0, sampling_rate//2)
        plt.stem(delta_f, np.abs(X))

        u = np.ones(len(h))
        y = np.convolve(u, h)
        n=np.arange(len(y))
        plt.subplot(326)
        plt.xlim(0, len(y)//2)
        plt.title('Step Response')
        plt.stem(n, y)
        plt.tight_layout()

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
