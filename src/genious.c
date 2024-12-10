#include <math.h>
#include "csdl.h"
#include "csoundCore.h"

#include "genious.h"

static int32_t lorenztable(FGDATA *ff, FUNC *ftp)
{
  /*
    p3 = size
    p4 = GEN routine
    p5 = choose axis to use values of; 0 = x, 1 = y, 2 = z
    p6 = min output; if p6 and p7 == 0 don't scale
    p7 = max output; if p6 and p7 == 0 don't scale
    p8 = norm; 0 == normalize, -1 == don't normalize
    p9 = stepsize
    p10 = x start
    p11 = y start
    p12 = z start
    p13 = sigma; 0 == default -> 10
    p14 = rho; 0 == default -> 28 
    p15 = beta; 0 == default -> 8./3.
    p16 = time delta; 0 == default -> 0.001
  */
  
  /* function table */
  MYFLT *fp = ftp->ftable;
  MYFLT sigma, rho, beta, time;
  
  /* default arguments */
  if (ff->e.p[13] == 0){
    sigma = 10;
  }
  else {
    sigma = ff->e.p[13];
  } 
    
  if (ff->e.p[14] == 0){
    rho = 28;
  }
  else {
    rho = ff->e.p[14];
  }
    
  if (ff->e.p[15] == 0){
    beta = 8.0/3.0;
  }
  else {
    beta = ff->e.p[15];
  }
     
  if (ff->e.p[16] == 0){
    time = 0.001;
  }
  else {
    time = ff->e.p[16];
  }

  /* p-field arguments */
  int32_t ft_size = ftp->flen;
  MYFLT min_out = ff->e.p[6];
  MYFLT max_out = ff->e.p[7];
  MYFLT step_size = ff->e.p[9];
  MYFLT x = ff->e.p[10];
  MYFLT y = ff->e.p[11];
  MYFLT z = ff->e.p[12];
  MYFLT norm = ff->e.p[8];  

  MYFLT xx, yy;
  /* write raw data to ft */
  for (int i = 0; i < ft_size; i++)
    {    
      if (ff->e.p[5] == 0){
	fp[i] = x;
      } else if (ff->e.p[5] == 1){
	fp[i] = y;
      } else if (ff->e.p[5] == 2){
	fp[i] = z;
      }
      
      for (int j = 0; j < step_size; j++) {
	xx = (sigma * (y - x))       * time;
	yy = (x     * (rho - z) - y) * time;
	z  += ((x * y) - (beta * z))  * time;
	x  += xx;
	y  += yy;	  
      } 
    }

  /* scale data of ft */
  if((min_out != 0.0) || (max_out != 0.0))
    scale_array(fp, min_out, max_out, ft_size);

  if(norm!=0.0)
    ff->e.p[4] = -1;
      
  return OK;
}

static int32_t thomastable(FGDATA *ff, FUNC *ftp)
{
  /*
    p3 = size
    p4 = GEN routine
    p5 = choose axis to use values of; 0 = x, 1 = y, 2 = z
    p6 = min output; if p6 and p7 == 0 don't scale
    p7 = max output; if p6 and p7 == 0 don't scale
    p8 = norm; 0 == normalize, -1 == don't normalize
    p9 = stepsize
    p10 = x start
    p11 = y start
    p12 = z start
    p13 = b constant value; 0 == default to 0.9
    p14 = time delta; 0 == default -> 0.001
  */
  
  /* function table */
  MYFLT *fp = ftp->ftable;
  MYFLT b, time;
  
  /* default arguments */
  if (ff->e.p[13] == 0){
    b = 0.9;
  } else {
    b = ff->e.p[13];
  } 
  if (ff->e.p[14] == 0){
    time = 0.001;
  } else {
    time = ff->e.p[14];
  }

  /* p-field arguments */
  int32_t ft_size = ftp->flen;
  MYFLT min_out = ff->e.p[6];
  MYFLT max_out = ff->e.p[7];
  MYFLT step_size = ff->e.p[9];
  MYFLT x = ff->e.p[10];
  MYFLT y = ff->e.p[11];
  MYFLT z = ff->e.p[12];
  MYFLT norm = ff->e.p[8];

  MYFLT xx, yy;
  
  /* write raw data to ft from the thomas attractor
     https://en.wikipedia.org/wiki/Thomas'_cyclically_symmetric_attractor
   */
  for (int i = 0; i < ft_size; i++)
    {    
      if (ff->e.p[5] == 0){
	fp[i] = x;
      } else if (ff->e.p[5] == 1){
	fp[i] = y;
      } else if (ff->e.p[5] == 2){
	fp[i] = z;
      }
      
      for (int j = 0; j < step_size; j++) {
	xx =  ((-b * x) + sin(y)) * time;
	yy =  ((-b * y) + sin(z)) * time;
	z  +=  ((-b * z) + sin(x)) * time;
	x  += xx;
	y  += yy;
      } 
    }

  /* scale data of ft */
  if((min_out != 0.0) || (max_out != 0.0))
    scale_array(fp, min_out, max_out, ft_size);

  if(norm!=0.0)
    ff->e.p[4] = -1;
      
  return OK;
}

static NGFENS localfgens[] = {
  { "lorenz", lorenztable},
  { "thomas", thomastable},
  { NULL, NULL}
};


FLINKAGE_BUILTIN(localfgens)
