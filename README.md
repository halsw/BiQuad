# Bi-quadratic filters for Arduino and Teensy
A library with template classes for the biquad and one pole filters

The biquad filters are the preferred solution for the implementation of high order IIR filters because of their stability and replacement for FIR filters because of their speed and memory usage.

In the library there are 4 distinct implementations of the biquad filters, the direct forms I & II, and the transposed direct forms I & II   


## Example file
The included .ino file provides an example usage and also helps to measure execution time and performance of the filters. The performance is tested by comparing the filter output to a noisy fast changing signal that is recovered through low pass filtering and can also be plotted through the Arduino IDE Serial Plotter.

To compile the .ino file you must first install the **FixedPoints** library

**PLOT** defined parameter controls the program output as 
0  print execution time of all forms
1  plot direct form I floating point implementation
2  plot direct form I fixed point implementation
3  plot direct form II floating point implementation
4  plot direct form II fixed point implementation
5  plot transposed direct form I floating point implementation
6  plot transposed direct form I fixed point implementation
7  plot transposed direct form II floating point implementation
8  plot transposed direct form II fixed point implementation
9  plot cascaded direct form I floating point implementation
10 plot cascaded direct form I fixed point implementation
11 plot cascaded direct form II floating point implementation
12 plot cascaded direct form II fixed point implementation
13 plot cascaded transposed direct form I floating point implementation
14 plot cascaded transposed direct form I fixed point implementation
15 plot cascaded transposed direct form II floating point implementation
16 plot cascaded transposed direct form II fixed point implementation

**TFLOAT** defined parameter controls the type used for the floating point implementation

**TFIXED** defined parameter controls the type used for the floating point implementation

**CASCADE** defined parameter controls the number of biquads to cascade, the order of the IIR filter is two times this parameter

**PLOT_FLAGS** defined parameter controls what gets plotted

You may also want to change the test signal by changing its Fourier components **in const double amp\[10]** (phase is not included but his gives it a periodic steep transition) and also the noise by changing the **NOISE** defined parameter

## Filter Coefficients
There are many ways to calculate the coefficients. The library included the biquadCoefs() function that uses predefined filters from  the "Cookbook formulae for audio EQ biquad filter coefficients"

The filter coefficients may also be calculated from S-domain transfer functions of analog filters 
        * 1 with    (1+2.z^-1+z^-2).(1-cos(w0))
        * S with    (1-z^-2).sin(w0)
        * S^2 with   (1-2.z^-1+z^-2).(1+cos(w0))
        * 1+S^2 with 2.(1-2.cos(w0).z^-1+z^-2)
where w0 is the resonance angular frequency w0=2.π.frequency.sampling period  

Another way is to factorize the Z-domain transfer function H(z) = A(z)/B(z) nominator and denominator and split them to biquads Hi(z)=Ai(z)/Bi(z)=g.(1-z0).(1-z1)/{(1-p0).(1-p1)}, for conjugate poles it is more convenient to use the polar notation of radius *R* and angle *θ*. In this case angle *θ* is considered the resonance frequency w0 and *R* controls the Q factor of the filter.

There are also two functions in the library, one is bandwith2Q() that converts bandwidth (between 3dB cutoff points) to Q factor and the other shelfSlope2Q() that converts shelving slope to Q factor. The shelf slope is the slope of high and low shelf filters and has a maximum of S=1. 