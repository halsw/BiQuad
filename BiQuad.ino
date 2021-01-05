/*
 * This file is part of the BiQuad library
 * Usage: Provide an example use of the library
 *        and test performance of biquad filters
 * Dependencies: FixedPoints library
 * 
 * Version 1.0.1
 * Developed by Evan https://github.com/halsw
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "BiQuad.h"
#include <limits.h> 
#include <FixedPoints.h>
#include <FixedPointsCommon.h>

//Define the types used for floating point and fixed point calculations
#define TFLOAT double
#define TFIXED SQ3x4

//Define the cycles run for profiling the filters
#define CYCLES 1000

//Define the plot values to output 0=profiling 1=fltDF1 2=fxdDF1 3=fltDF2 4=fxdDF2 5=fltTF1 6=fxdTF1 7=fltTF2 8=fxdTF2 9-16=cascaded versions
#define PLOT 16

//Define the number of filters to cascade (if PLOT>8)
#define CASCADE 3

//Define what to plot 0x1=Signal 0x2=Noisy Signal 0x4=Filter Output 0x8=Error(x100) 0x10=Noise 0x20=Ececution time(ms) 0x40=Load
#define PLOT_FLAGS 0x1D

//Define the noise standard deviation
#define NOISE 2.0

//Define the period in ms
#define PERIOD_MS 5

#define PROFILE(X,T) double e##X, a##X;\ 
  static double x##X = 0.0, m##X = 0.0, v##X = 0.0;\
  T i##X = noisy;\  
  unsigned long t##X = micros();\
  T r##X = X.filter(i##X);\
  t##X = micros() - t##X;\
  x##X += (t##X - x##X)/sample;\
  e##X = static_cast<double>(r##X) - value - m##X;\
  a##X = e##X / sample;\
  if (sample>1) v##X -= v##X / ( sample - 1);\
  v##X += e##X * a##X;\
  m##X += a##X

#define PPRINT(X) Serial.print(#X ": ");\
  Serial.print(x##X);\
  Serial.print("us error:");\ 
  Serial.println(sqrt(v##X),5);\
  m##X = v##X = x##X = 0.0;

#define EXECMS(X)   unsigned long t0 = micros();\
  do {X;} while (0);\
  float execms = (micros() - t0) / 1000.0

BQFilterDFI<TFLOAT> float_direct_form1(1.0, 2.0, 1.0, 2.0, 1.0);
BQFilterDFI<TFIXED> fixed_direct_form1(1.0, 2.0, 1.0, 2.0, 1.0);
BQFilterDFII<TFLOAT> float_direct_form2(1.0, 2.0, 1.0, 2.0, 1.0);
BQFilterDFII<TFIXED> fixed_direct_form2(1.0, 2.0, 1.0, 2.0, 1.0);
BQFilterTFI<TFLOAT> float_xposed_form1(1.0, 2.0, 1.0, 2.0, 1.0);
BQFilterTFI<TFIXED> fixed_xposed_form1(1.0, 2.0, 1.0, 2.0, 1.0);
BQFilterTFII<TFLOAT> float_xposed_form2(1.0, 2.0, 1.0, 2.0, 1.0);
BQFilterTFII<TFIXED> fixed_xposed_form2(1.0, 2.0, 1.0, 2.0, 1.0);
BQFilter<TFLOAT> *float_cascade;
BQFilter<TFIXED> *fixed_cascade;

const double amp[10] = {1.0, 0.988, 0.951, 0.891, 0.809, 0.707, 0.588, 0.454, 0.309, 0.156}; // Signal generator amplitudes
const double f0 = 50.0 / PERIOD_MS; //The signal generator base frequency
const double dt = PERIOD_MS / 1000.0; //Time increment

unsigned int wait() {
  static unsigned int load = 0;
  static unsigned long loopMS=PERIOD_MS;
  unsigned long now=millis();
  if (now>loopMS) {
    load = 100<<8;
    loopMS = now + PERIOD_MS;
  } else {
    load += (25600 - (loopMS-now)*(25600/PERIOD_MS) - load)>>3;
    while (loopMS>millis());
    loopMS += PERIOD_MS;
  }
  return(load);
}

void setup() {
  Serial.begin(115200);
  randomSeed(analogRead(0)); //assuming A0 is not connected
  #if PLOT == 0
  Serial.print("Biquad filters example, please wait ");
  Serial.print((unsigned long)CYCLES*PERIOD_MS/1000);
  Serial.println("s between updates...");
  #elif PLOT == 9
  float_cascade = createBQFilter<TFLOAT>(BQ_DIRECT1, BQ_LOPASS, dt, 10*f0, 0.0, 0.0, CASCADE);
  #elif PLOT == 10
  fixed_cascade = createBQFilter<TFIXED>(BQ_DIRECT1, BQ_LOPASS, dt, 10*f0, 0.0, 0.0, CASCADE);
  #elif PLOT == 11
  float_cascade = createBQFilter<TFLOAT>(BQ_DIRECT2, BQ_LOPASS, dt, 10*f0, 0.0, 0.0, CASCADE);
  #elif PLOT == 12
  fixed_cascade = createBQFilter<TFIXED>(BQ_DIRECT2, BQ_LOPASS, dt, 10*f0, 0.0, 0.0, CASCADE);
  #elif PLOT == 13
  float_cascade = createBQFilter<TFLOAT>(BQ_TRANSPOSED1, BQ_LOPASS, dt, 10*f0, 0.0, 0.0, CASCADE);
  #elif PLOT == 14
  fixed_cascade = createBQFilter<TFIXED>(BQ_TRANSPOSED1, BQ_LOPASS, dt, 10*f0, 0.0, 0.0, CASCADE);
  #elif PLOT == 15
  float_cascade = createBQFilter<TFLOAT>(BQ_TRANSPOSED2, BQ_LOPASS, dt, 10*f0, 0.0, 0.0, CASCADE);
  #elif PLOT == 16
  fixed_cascade = createBQFilter<TFIXED>(BQ_TRANSPOSED2, BQ_LOPASS, dt, 10*f0, 0.0, 0.0, CASCADE);
  #endif
}

void loop() {
  static int sample = 1;
  int i, load;
  double f, noisy, value, recover;
  static double t=0;
  f = TWO_PI*f0*t;
  value = 0.0;
  for (i=0; i<10; i++) value += amp[i]*sin(i*f); 
  noisy = value+NOISE*((LONG_MAX>>1)-random(LONG_MAX))/LONG_MAX;
  load = wait(); //Sample at exact intervals
  #if PLOT == 1
  EXECMS(recover=static_cast<double>(float_direct_form1.filter(value)));
  #elif PLOT == 2
  EXECMS(recover=static_cast<double>(fixed_direct_form1.filter(value)));
  #elif PLOT == 3
  EXECMS(recover=static_cast<double>(float_direct_form2.filter(value)));
  #elif PLOT == 4
  EXECMS(recover=static_cast<double>(fixed_direct_form2.filter(value)));
  #elif PLOT == 5
  EXECMS(recover=static_cast<double>(float_xposed_form1.filter(value)));
  #elif PLOT == 6
  EXECMS(recover=static_cast<double>(fixed_xposed_form1.filter(value)));
  #elif PLOT == 7
  EXECMS(recover=static_cast<double>(float_xposed_form2.filter(value)));
  #elif PLOT == 8
  EXECMS(recover=static_cast<double>(fixed_xposed_form2.filter(value)));
  #elif (PLOT == 9) || (PLOT == 11) || (PLOT == 13) || (PLOT == 15)
  EXECMS(recover=static_cast<double>(float_cascade->filter(value)));
  #elif (PLOT == 10) || (PLOT == 12) || (PLOT == 14) || (PLOT == 16)
  EXECMS(recover=static_cast<double>(fixed_cascade->filter(value)));
  #endif
  t+= PERIOD_MS/1000.0; //update simulation time
  #if PLOT == 0
  PROFILE(float_direct_form1, TFLOAT);
  PROFILE(fixed_direct_form1, TFIXED);
  PROFILE(float_direct_form2, TFLOAT);
  PROFILE(fixed_direct_form2, TFIXED);
  PROFILE(float_xposed_form1, TFLOAT);
  PROFILE(fixed_xposed_form1, TFIXED);
  PROFILE(float_xposed_form2, TFLOAT);
  PROFILE(fixed_xposed_form2, TFIXED);
  if (sample++ > CYCLES) {
    sample = 1;
    Serial.print("Load(%):");
    Serial.println(load/256.0);
    PPRINT(float_direct_form1);
    PPRINT(fixed_direct_form1);
    PPRINT(float_direct_form2);
    PPRINT(fixed_direct_form2);
    PPRINT(float_xposed_form1);
    PPRINT(fixed_xposed_form1);
    PPRINT(float_xposed_form2);
    PPRINT(fixed_xposed_form2);
  }
  #else
  #if (PLOT_FLAGS & 0x1) != 0
  Serial.print(value,6);
  Serial.print(" ");
  #endif  
  #if (PLOT_FLAGS & 0x2) != 0
  Serial.print(noisy,6);
  Serial.print(" ");
  #endif  
  #if (PLOT_FLAGS & 0x4) != 0
  Serial.print(recover,6);
  Serial.print(" ");
  #endif  
  #if (PLOT_FLAGS & 0x08) != 0
  Serial.print(100.0*(recover-value),3);
  Serial.print(" ");
  #endif
  #if (PLOT_FLAGS & 0x10) != 0
  Serial.print(noisy-value,3);
  Serial.print(" ");
  #endif
  #if (PLOT_FLAGS & 0x20) != 0
  Serial.print(execms,3);
  Serial.print(" ");
  #endif
  #if (PLOT_FLAGS & 0x40) != 0
  Serial.print(load/25600.0);
  #endif
  Serial.println(" ");
  #endif
}
