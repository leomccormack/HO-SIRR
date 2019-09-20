# HO-SIRR

Higher-order Spatial Impulse Response Rendering (HO-SIRR) is a rendering method, which can generate output loudspeaker array impulse responses (IRs) from input spherical harmonic (Ambisonic) IRs. The method makes assumptions regarding the composition of the sound-field and extracts parameters over time, which allows the HOSIRR renderer to map in the input to the output in an adaptive and more informed manner when compared to linear methods; such as Ambisonics. 

The idea is that you then convolve a monophonic source with this loudspeaker array IR, and it is not only reproduced in the direction relative to the source/receiver during the capturing of the Ambisonic IRs, but it should also exhibit all of the spatial characteristics of the captured space.

**(Please note that the code is currently a "work in progress". Feedback is most welcome!)**

![](SIRR_vs_Ambi_vs_Reference.png)

## Getting Started

The code is reliant on the following Matlab libraries:
* [Spherical-Harmonic-Transform](https://github.com/polarch/Spherical-Harmonic-Transform)
* [Higher-Order-Ambisonics](https://github.com/polarch/Higher-Order-Ambisonics)
* [Vector-Base-Amplitude-Panning](https://github.com/polarch/Vector-Base-Amplitude-Panning)

The HO-SIRR script may be configured, for example, as:
```
pars.order = 3; % transform order of input signals  
pars.fs = 48e3; 
pars.ls_dirs_deg = ???; % your loudspeaker directions, [azi elev], in degrees

% Specify windowing sizes for multi-resolution STFT, for example:
pars.multires_winsize = [1024 128 16]; 
pars.multires_xovers = [500 2e3]; % 1024 up to 500Hz, then 128 up to 2kHz, then 16 past 2kHz (Note HOSIRR uses 50% window overlap)

% Diffuse rendering parameters
pars.RENDER_DIFFUSE = 1;     % 0: disable diffuse rendering, 1: enable
pars.decorrelationType = 'noise'; % {'phase','noise'}, decorrelation via convolution with 'noise', or via randomising the 'phase'
pars.maxDiffuseAnalysis_Hz = 6e3; % frequency up to which to estimate the diffuseness parameter 
pars.alpha_diff = 0.975;     % minimum diffuseness averaging coefficient (one-pole filter)

% Optionally, the first highest N peaks of the response may be isolated
% and panned with broad-band DoA estimates, which can reduce timbral 
% colourations of the method
pars.BROADBAND_DIRECT = 1;   % 0: disabled, 1: enabled
pars.nBroadbandPeaks = 1;    % number of peaks to isolate 
```

The input Ambisonic signals may be rendered for the loudspeaker set-up as:

```
pars.chOrdering = 'ACN'; % 'ACN', or 'WXYZ'  
pars.normScheme = 'N3D'; % 'N3D', or 'SN3D'
shir = audioread( ??? ) % load input audio
[lsir, lsir_ndiff, lsir_diff, pars, analysis] = HOSIRR(shir, pars);
% lsir       - output loudspeaker impulse responses
% lsir_ndiff - output loudspeaker impulse responses, direct stream only
% lsir_diff  - output loudspeaker impulse responses, diffuse stream only
% pars       - output parameters
% analysis   - analysed parameters stored during rendering
```

## Developers

* **Leo McCormack** - Matlab and algorithm design (contact: leo.mccormack@aalto.fi)
* **Archontis Politis** - Matlab and algorithm design
* **Ville Pulkki** - Matlab and algorithm design

## License

This code is provided under the [BSD 3-clause license](https://opensource.org/licenses/BSD-3-Clause). 

## References 

[1] McCormack, L., Politis, A., Scheuregger, O., and Pulkki, V. (2019). "**Higher-order processing of spatial impulse responses**".
In Proceedings of the 23rd International Congress on Acoustics, 9--13 September 2019 in Aachen, Germany.

[2] Politis, A. and Pulkki, V., 2016. "**Acoustic intensity, energy-density and diffuseness estimation in a directionally-constrained region**". 
arXiv preprint arXiv:1609.03409.

[3] Merimaa, J. and Pulkki, V., 2005. "**Spatial impulse response rendering I: Analysis and synthesis**". 
Journal of the Audio Engineering Society, 53(12), pp.1115-1127.

[4] Pulkki, V. and Merimaa, J., 2006. "**Spatial impulse response rendering II: Reproduction of diffuse sound and listening tests**". 
Journal of the Audio Engineering Society, 54(1/2), pp.3-20.

[5] Favrot, S. and Buchholz, J.M., 2010. "**LoRA: A loudspeaker-based room auralization system. Acta Acustica united with Acustica**", 96(2),  pp.364-375.
