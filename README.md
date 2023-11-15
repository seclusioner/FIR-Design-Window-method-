# FIR Design (Window Method)

In this section, we implement a Python FIR (Finite Impulse Response) filter design class using the window method. We use a design table based on the DSP textbook to construct this class.

## Window Functions
We consider the following window functions:
- Rectangular
- Hanning
- Hamming
- Blackman
- Kaiser

### Filter Design Parameters

| Window        | Transition  | Passband | Main/Side Lobe   | Stopband          |
|---------------|-------------|----------|------------------|-------------------|
| Function      | Width (Hz)  | Ripple   | Relation (dB)    | Attenuation (dB)  |
| Rectangular   | 0.9 / N     | 0.7416   | 13               | 21                |
| Hanning       | 3.1 / N     | 0.0546   | 31               | 44                |
| Hamming       | 3.3 / N     | 0.0194   | 41               | 53                |
| Blackman      | 5.5 / N     | 0.0017   | 57               | 74                |
| Kaiser        | 2.93 / N    | 0.0274   | <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>β</mi></math> = 4.54        | 50               |
|               | 4.32 / N    | 0.00275  | <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>β</mi></math> = 6.76        | 70               |
|               | 5.71 / N    | 0.000275 | <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>β</mi></math> = 8.96        | 90               |

## Filter Types
We consider the following filter types:
- Lowpass (LPF)
- Highpass (HPF)
- Bandpass (BPF)
- Bandreject (BRJF)

### Ideal Impulse Responses

| Type of Filter | Ideal Impulse Response               |                       |
|----------------|--------------------------------------|-----------------------|
|                | hd[n], n <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mo>≠</mo></math> 0 | hd[0] |
| Lowpass        | Formula 1                            | 2fc               |
| Highpass       | Formula 2                            | 1-2fc             |
| Bandpass       | Formula 3                            | 2(f2-f1)          |
| Bandreject     | Formula 4                            | 1-2(f2-f1)        |

## Formulas


Formula 1:

$$ 2fc \cdot \frac{\sin(nwc)}{nwc}  $$

Formula 2:

$$ -2fc \cdot \frac{\sin(nwc)}{nwc} $$

Formula 3:

$$ 2f_{2}  \cdot \frac{\sin(nwc_{2})}{nwc2}-2f_{1} \cdot \frac{\sin(nwc_{1})}{nwc_{1}} $$

Formula 4:

$$ 2f_{1} \cdot \frac{\sin(nwc_{1})}{nwc_{1}}-2f2 \cdot \frac{\sin(nwc2)}{nwc2} $$

## Application

$$ y[n] = \sum_{k=0}^{N} h[k] \cdot x[n-k] $$

<hr>
<br>
<div align="right">
    <p><strong> Authors: D.S. 、 Yuyu378 </strong></p>
</div>

