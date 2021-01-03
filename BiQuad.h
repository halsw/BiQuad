/*
 * This file is part of the BiQuad library
 * Usage: A template library for the implementation
 *        of biqudratic filters for arduino/teensy
 * Version 1.0
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
 * 
 * Classes:
 *   BQFilter     the one pole filter first order filter, not a biquadratic one but maybe useful for odd order cascaded filters
 *   BQFilterDFI  the biquadratic filter direct form 1
 *   BQFilterDFII the biquadratic filter direct form 2
 *   BQFilterTFI  the biquadratic filter transposed direct form 1
 *   BQFilterTFII the biquadratic filter transposed direct form 2
 *    
 * Functions: (these are seperate from classes so that they won't compile if not needed)
 *   void firstOrdCoefs() Calculate one pole filter coeffcicients in float array coefs[5] based on design parameters 
 *   BQFilter<typedef> *createFOFilter() Create one pole filter based on design parameters (its ether low or high pass)
 *   float bandwith2Q() convert bandwith to Q factor, for use in biquadCoefs()
 *   float shelfSlope2Q() convert shelving slope to Q factor, for use in biquadCoefs()
 *   void biquadCoefs() Calculate biquad filter coeffcicients in float array coefs[5] based on design parameters
 *   BQFilter<typedef>* createBQFilter() Create biquad filter based on design parameters
 *   void configBQFilter() Change on the fly biquad filter coefficients based on design parameters
 *   
 * Transfer functions:
 *                         b0
 * First order H(z) = ------------
 *                     1 + a1.z^-1
 *                     
 *                    b0 + b1.z^-1 + b2^z^-2
 * Biquadratic H(z) = ----------------------
 *                     1 + a1.z^-1 + a2^z^-2
 *                     
 * Conversion from S-domain to Z-domain (by substituting in analog filter transfer functions one can compute the biquad coefficients b0,b1,b2,a1,a2                   
 *                    1        (1+2.z^-1+z^-2).(1-cos(w0))
 *                    S        (1-z^-2).sin(w0)
 *                   S^2       (1-2.z^-1+z^-2).(1+cos(w0))
 *                  1+S^2      2.(1-2.cos(w0).z^-1+z^-2) 
 *                  
 * The coefficients may also be calculated by factoring the transfer function nominator and denominator and splitting them ro biquad pairs
 */
 
#ifndef BIQUAD_H
#define BIQUAD_H

//****** The biquad filter types for biquadCoefs() *********
#define BQ_LOPASS 0
#define BQ_HIPASS 1
#define BQ_LOSHELF 2
#define BQ_HISHELF 3
#define BQ_BANDPASS 4
//constant skirt gain S, use shelfSlope2Q() for S to Q conversion
#define BQ_BANDPASS_SKIRT 5
//0db peak gain
#define BQ_BANDPASS_0DB 6
#define BQ_NOTCH 7
#define BQ_INOTCH 8
#define BQ_ALLPASS 9
//**********************************************************

//***** The biquad filter forms for createBQFilter() *******
#define BQ_FIRSTORDER 0
#define BQ_DIRECT1 1
#define BQ_DIRECT2 2
#define BQ_TRANSPOSED1 3
#define BQ_TRANSPOSED2 4
//**********************************************************

template <class BQ> class BQFilterDFI;
template <class BQ> class BQFilterDFII;
template <class BQ> class BQFilterTFI;
template <class BQ> class BQFilterTFII;

template <class BQ>
class BQFilter {
  friend class BQFilterDFI<BQ>;
  friend class BQFilterDFII<BQ>;
  friend class BQFilterTFI<BQ>;
  friend class BQFilterTFII<BQ>;
  protected:
    char nam[8];
    static int registry;    
    BQFilter *cascd;  
    BQ b0;
    BQ b1;
    BQ b2;
    BQ a1;
    BQ a2;
    BQ y[2] = {0.0};    
  public:
    BQFilter (BQ b0, BQ b1, BQ b2=0.0, BQ a1=0.0, BQ a2=0.0, char* fname=NULL, BQFilter<BQ> *cascd=NULL) {
      this->cascd = cascd;
      this->b0 = b0;
      this->b1 = b1;
      this->b2 = b2;
      this->a1 = a1;
      this->a2 = a2;
    }
    BQFilter (BQ *c, char* fname=NULL, BQFilter<BQ> *cascd=NULL) {
      this->cascd = cascd;
      setCoefs(c);
    }
    virtual BQ filter( BQ value) {
      y[0] = value * b0 - y[0] * b1;
      if (cascd) return cascd->filter(y[0]);
      return y[0];
    }    

    //Get a cascaded filter 
    BQFilter* cascade(unsigned char order=0) {
      BQFilter* f=this;
      while (order-- && f) f=f->cascd;
      return(f);
    }    

    //Get the number of cascaded filters
    unsigned char cascaded() {
      unsigned char i=1;
      BQFilter* f=this;
      while (f=f->cascd) i++;
      return(i);
    }    

    //Set the filter coefficients
    void setCoefs(BQ* c, unsigned char order=0) {
      BQFilter* f=cascade(order);
      f->b0 = *c++;
      f->b1 = *c++;
      f->b2 = *c++;
      f->a1 = *c++;
      f->a2 = *c;
    }

    //Set cascaded filter
    inline void setCascaded(BQFilter* cascd, unsigned char order=0) {
      cascade(order)->cascd=cascd;
    }

    //Get cascaded filter
    inline BQFilter* getCascaded(unsigned char order=0) {
      return cascade(order)->cascd;
    }

    //Set the name of the filter
    inline void setName(char* nam, unsigned char order=0) {
      strncpy(cascade(order)->nam,8,nam);
    }

    //Set the name of the filter from a String type
    void setName(String nam, unsigned char order=0) {
      setName(const_cast<char*>(nam.c_str()), order);
    }

    //Get the name of the filter
    const char* getName(unsigned char order=0) {
      return(cascade(order)->nam);
    }

    //Get the filter gain (the gain is not in Db to avoid exp() function that may not be available for <BQ> type) 
    BQ getGain(bool cascaded=false) {
      BQFilter* f=this;
      BQ g = 1.0;
      do {
        if (f->b2 == 0.0 && f->a1 == 0.0 && f->a2 == 0.0)
          g *= f->b0 / (1.0 + f->b1);
        else
          g *= (f->b0 + f->b1 + f->b2) / (1.0 + f->a1 + f->a2);
        f = f->cascd;
      } while ( f && cascaded);    
      return; 
    }
  };

template <class BQ> int BQFilter<BQ>::registry = 0;

template <class BQ>
class BQFilterDFII: public BQFilter<BQ> {
  public:
    BQFilterDFII (BQ b0, BQ b1, BQ b2, BQ a1, BQ a2, char* fname=NULL, BQFilter<BQ> *cascd=NULL) :
      BQFilter<BQ> (b0, b1, b2, a1, a2, fname, cascd) {};
    BQFilterDFII (BQ *c, char* fname=NULL, BQFilter<BQ> *cascd=NULL) : BQFilter<BQ> (c, fname, cascd) {};
    virtual BQ filter( BQ value) {
      value -= this->a1*this->y[0] + this->a2*this->y[1];
      value = this->b0*value + this->b1*this->y[0] + this->b2*this->y[1];
      this->y[1] = this->y[0];
      this->y[0] = value;
      if (this->cascd) return this->cascd->filter(value);
      return value; 
    }    
  };

template <class BQ>
class BQFilterDFI: public BQFilter<BQ> {
  friend class BQFilterTFI<BQ>;
  friend class BQFilterTFII<BQ>;
  protected:
    BQ x[2]= {0.0};    
  public:
    BQFilterDFI (BQ b0, BQ b1, BQ b2, BQ a1, BQ a2, char* fname=NULL, BQFilter<BQ> *cascd=NULL) :
      BQFilter<BQ> (b0, b1, b2, a1, a2, fname, cascd) {};
    BQFilterDFI (BQ *c, char* fname=NULL, BQFilter<BQ> *cascd=NULL) : BQFilter<BQ> (c, fname, cascd) {};
    virtual BQ filter( BQ value) {
      BQ out =  this->b0*value + this->b1*this->x[0] + this->b2*this->x[1] - this->a1*this->y[0] - this->a2*this->y[1];
      this->x[1] = this->x[0];
      this->x[0] = value;
      this->y[1] = this->y[0];
      this->y[0] = out;
      if (this->cascd) return this->cascd->filter(out);
      return out; 
    }    
  };

template <class BQ>
class BQFilterTFI: public BQFilterDFI<BQ> {
  public:
    BQFilterTFI (BQ b0, BQ b1, BQ b2, BQ a1, BQ a2, char* fname=NULL, BQFilter<BQ> *cascd=NULL) :
      BQFilterDFI<BQ> (b0, b1, b2, a1, a2, fname, cascd) {};
    BQFilterTFI (BQ *c, char* fname=NULL, BQFilter<BQ> *cascd=NULL) : BQFilterDFI<BQ> (c, fname, cascd) {};
    virtual BQ filter( BQ value) {
      BQ out;
      value += this->x[1];
      this->x[1] = this->x[0] - this->a1*value;
      this->x[0] = -this->a2*value;
      out = this->b0*value + this->y[1];
      this->y[1] = this->y[0] + this->b1*value;
      this->y[0] = this->b2*value;
      if (this->cascd) return this->cascd->filter(out);
      return out; 
    }    
  };

template <class BQ>
class BQFilterTFII: public BQFilterDFI<BQ> {
  public:
    BQFilterTFII (BQ b0, BQ b1, BQ b2, BQ a1, BQ a2, char* fname=NULL, BQFilter<BQ> *cascd=NULL) :
      BQFilterDFI<BQ> (b0, b1, b2, a1, a2, fname, cascd) {};
    BQFilterTFII (BQ *c, char* fname=NULL, BQFilter<BQ> *cascd=NULL) : BQFilterDFI<BQ> (c, fname, cascd) {};
    virtual BQ filter( BQ value) {
      BQ out = this->b0*value + this->y[1];
      this->y[1] = this->y[0];
      this->y[0] = this->x[1] + this->b1*value - this->a1*out;
      this->x[1] = this->x[0];
      this->x[0] = this->b2*value - this->a2*out;
      if (this->cascd) return this->cascd->filter(out);
      return out; 
    }    
  };

//Calculate first order filter coeffcicients in float array coefs[5]
void firstOrdCoefs(float* coefs, unsigned char type, float period, float frequency) {
  float a1, b0, w0 = TWO_PI*frequency*period;
  switch (type % 2) {
    case BQ_LOPASS: {
      a1 = -exp(-w0);
      b0 = 1.0 + a1;
      break;
    }
    case BQ_HIPASS: {
      a1 = exp(w0 - PI);
      b0 = 1.0 - a1;
      break;
    }
  }  
  *coefs++ = b0;
  *coefs++ = a1;
  *coefs++ = 0.0;
  *coefs++ = 0.0;
  *coefs = 0.0;
}

//Create first order filter based on design parameters
template <class BQ>
BQFilter<BQ> *createFOFilter(unsigned char type, float period, float frequency) {
  float coefs[5];
  BQ coeff[5];
  int i;
  firstOrdCoefs(coefs, type, period, frequency);
  for (i=0; i<5; i++) coeff[i] = coefs[i];
  return new BQFilter<BQ>(coeff);
}

//Convert bandwith to Q factor
float bandwith2Q(float bw, float period, float frequency) {
  const float w0 = TWO_PI*frequency*period;
  return 0.5/sinh(M_LN2/2.0*bw*w0/sin(w0));
}

//Convert shelving slope to Q factor
float shelfSlope2Q(float S, float gainDb) {
  const float A = pow(10, gainDb / 40.0);
  return 1.0/sqrt((A+1.0/A)*(1/S-1.0)+2.0);
}

//Calculate biquad filter coeffcicients in float array coefs[5]
void biquadCoefs(float* coefs, unsigned char type, float period, float frequency, float Q=0.707, float gainDb=0.0) {
  const float A = pow(10, gainDb / 40.0);
  const float w0 = TWO_PI*frequency*period;
  const float s = sin(w0);
  const float c = cos(w0);
  const float a = s/(2.0*Q);
  const float b = 2.0*sqrt(A);
  float b0, b1, b2, a0, a1, a2;
  switch (type) {
    case BQ_LOPASS: {
      b1 = 1.0 - c;
      b0 = b2 = b1 / 2.0;
      a0 = 1.0 + a;
      a1 =-2.0*c;
      a2 = 1.0 - a;
      break;
    }
    case BQ_HIPASS: {
      b1 = - 1.0 - c;
      b0 = b2 = b1 / -2.0;
      a0 = 1.0 + a;
      a1 = -2.0*c;
      a2 = 1.0 - a;
      break;
    }
    case BQ_LOSHELF: {
       b0 =  A*((A+1.0) - (A-1.0)*c + b*a);
       b1 =  2.0*A*((A-1.0) - (A+1.0)*c);
       b2 =  A*((A+1.0) - (A-1,0)*c - b*a);
       a0 =  (A+1.0) + (A-1.0)*c + b*a;
       a1 = -2.0*((A-1.0) + (A+1.0)*c);
       a2 =  (A+1.0) + (A-1.0)*c - b*a;
      break;
    }
    case BQ_HISHELF: {
       b0 =  A*((A+1.0) - (A-1.0)*c + b*a);
       b1 = -2.0*A*((A-1.0) - (A+1.0)*c);
       b2 =  A*((A+1.0) - (A-1,0)*c - b*a);
       a0 =  (A+1.0) - (A-1.0)*c + b*a;
       a1 =  2.0*((A-1.0) - (A+1.0)*c);
       a2 =  (A+1.0) - (A-1.0)*c - b*a;
      break;
    }
    case BQ_BANDPASS: {
       b0 = Q*a;
       b1 =  0.0;
       b2 = -b0;
       a0 =  1.0 + A;
       a1 = -2.0*c;
       a2 =  1.0 - A;
      break;
    }
    case BQ_BANDPASS_SKIRT: {
       b0 = b2 = -s / 2.0;
       b1 =  0.0;
       a0 =  1.0 + A;
       a1 = -2.0*c;
       a2 =  1.0 - A;
      break;
    }
    case BQ_BANDPASS_0DB: {
       b0 =  A;
       b1 =  0.0;
       b2 = -A;
       a0 =  1.0 + A;
       a1 = -2.0*c;
       a2 =  1.0 - A;
      break;
    }  
    case BQ_NOTCH: {
       b0 = b2 = 1.0;
       b1 = a1 = -2.0*c;
       a0 =  1.0 + A;
       a2 =  1.0 - A;
      break;
    }
    case BQ_INOTCH: {
       b0 =  1.0 + a*A;
       b1 = a1 = -2.0*c;
       b2 =  1.0 - a*A;
       a0 =  1.0 + a/A;
       a2 =  1.0 - a/A;
      break;
    }
    case BQ_ALLPASS: {
       b0 = a2 = 1.0 - a;
       b1 = a1 = -2.0*c;
       b2 = a0 = 1.0 + a;
      break;
    }
    default: {
       a0 = 1.0;
       b0 = b1 = b2 = a1 = a2 = 0.0;
      break;
    }
  }
  *coefs++ = b0 / a0;
  *coefs++ = b1 / a0;
  *coefs++ = b2 / a0;
  *coefs++ = a1 / a0;
  *coefs = a2 / a0;
}

//Create biquad filter based on design parameters
template <class BQ>
BQFilter<BQ> *createBQFilter(unsigned char form, unsigned char type, float period, float frequency, float Q=0.707, float gainDb=0.0, int cascade=0) {
  float coefs[5];
  BQ coeff[5];
  float pole;
  int i;
  if (form == BQ_FIRSTORDER) return createFOFilter<BQ>(type, period, frequency);
  if (cascade) {
    BQFilter<BQ> *f, *r;
    pole = HALF_PI/(cascade<<1);
    for (i=0; i<cascade; i++) { 
      Q = 0.5 / cos((1+2*i)*pole); 
      f = createBQFilter<BQ>(form, type, period, frequency, Q, gainDb, 0);
      if ( i )
        r->setCascaded(f,i-1);
      else
        r = f;
    }
    return(r);  
  }
  biquadCoefs(coefs, type, period, frequency, Q, gainDb);
  for (i=0; i<5; i++) coeff[i] = coefs[i];
  switch (form) {
    case BQ_DIRECT1: return new BQFilterDFI<BQ>(coeff);
    case BQ_DIRECT2: return new BQFilterDFII<BQ>(coeff);
    case BQ_TRANSPOSED1: return new BQFilterTFI<BQ>(coeff);
    case BQ_TRANSPOSED2: return new BQFilterTFII<BQ>(coeff);
  }
  return NULL;    
}

//Configure a biquad filter based on design parameters
template <class BQ>
void configBQFilter(BQFilter<BQ> *f, bool cascadeConfig, unsigned char type, float period, float frequency, float Q=0.707, float gainDb=0.0) {
  float coefs[5];
  BQ coeff[5];
  int n;
  float pole;
  if (cascadeConfig) {
    n = f->cascaded();
    pole = HALF_PI/(n<<1);
    n = 0;
  }
  do {
    biquadCoefs(coefs, type, period, frequency, cascadeConfig?0.5 / cos((1+2*n++)*pole):Q, gainDb);
    for (int i=0; i<5; i++) coeff[i] = coefs[i];
    f->setCoefs(coeff);
    f=f->getCascaded();
  } while (f && cascadeConfig); 
}

void designPolesZeroes(float* coefs, float zero1, float zero2, bool zeroesConjugate, float pole1, float pole2, bool polesConjugate) {
  // find coefficients by evalucating (z - z1)*(z - z2)/(z - p1)*(z - p2)  
}
#endif
