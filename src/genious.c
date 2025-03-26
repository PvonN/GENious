#include "../include/genious.h"
#include "../include/genious_helper.h"

static int32_t lorenztable(FGDATA *ff, FUNC *ftp)
{
  /*
    p3 = size
    p4 = GEN routine
    p5 = choose axis to use values of; 0 = x, 1 = y, 2 = z
    p6 = min output; if p6 and p7 == 0 don't scale
    p7 = max output; if p6 and p7 == 0 don't scale
    p8 = norm; 0 == normalize, -1 == don't normalize
    p9 = x start
    p10 = y start
    p11 = z start
    p12 = sigma; 0 == default -> 10
    p13 = rho; 0 == default -> 28 
    p14 = beta; 0 == default -> 8./3.
    p15 = time delta; 0 == default -> 0.001
    p16 = stepsize; 0 == default -> 1
  */
  
  /* function table */
  MYFLT *fp = ftp->ftable;
  MYFLT sigma, rho, beta, time, step_size;
  
  /* default arguments */
  if (ff->e.p[12] == 0){
    sigma = 10;
  }
  else {
    sigma = ff->e.p[12];
  } 
    
  if (ff->e.p[13] == 0){
    rho = 28;
  }
  else {
    rho = ff->e.p[13];
  }
    
  if (ff->e.p[14] == 0){
    beta = 8.0/3.0;
  }
  else {
    beta = ff->e.p[14];
  }
     
  if (ff->e.p[15] == 0){
    time = 0.001;
  }
  else {
    time = ff->e.p[15];
  }

  if (ff->e.p[16] == 0){
    step_size = 1;
  }
  else {
    step_size = ff->e.p[16];
  }
  
  /* p-field arguments */
  int32_t ft_size = ftp->flen;
  MYFLT min_out = ff->e.p[6];
  MYFLT max_out = ff->e.p[7];
  MYFLT norm = ff->e.p[8];  
  MYFLT x = ff->e.p[9];
  MYFLT y = ff->e.p[10];
  MYFLT z = ff->e.p[11];


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
		xx = x + time * (sigma * (y - x));
		yy = y + time * (x * (rho - z) - y);
		z  = z + time * ((x * y) - (beta * z));
		x  = xx;
		y  = yy;	  
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
    p9 = x start
    p10 = y start
    p11 = z start
    p12 = b constant value; 0 == default to 0.2
    p13 = time delta; 0 == default -> 0.001
    p14 = stepsize; 0 == default -> 1
  */
  
  /* function table */
  MYFLT *fp = ftp->ftable;
  MYFLT b, time, step_size;
  
  /* default arguments */
  if (ff->e.p[12] == 0){
    b = 0.2;
  } else {
    b = ff->e.p[12];
  }
  
  if (ff->e.p[13] == 0){
    time = 0.01;
  } else {
    time = ff->e.p[13];
  }
  
  if (ff->e.p[14] == 0){
    step_size = 1;
  } else {
    step_size = ff->e.p[14];
  }
  
  
  /* p-field arguments */
  int32_t ft_size = ftp->flen;
  MYFLT min_out = ff->e.p[6];
  MYFLT max_out = ff->e.p[7];
  MYFLT norm = ff->e.p[8];
  MYFLT x = ff->e.p[9];
  MYFLT y = ff->e.p[10];
  MYFLT z = ff->e.p[11];


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
		xx = x + time * ((-b * x) + sin(y));
		yy = y + time * ((-b * y) + sin(z));
		z = z + time * ((-b * z) + sin(x));
		x = xx;
		y = yy;
      } 
    }

  /* scale data of ft */
  if((min_out != 0.0) || (max_out != 0.0))
    scale_array(fp, min_out, max_out, ft_size);

  if(norm!=0.0)
    ff->e.p[4] = -1;
      
  return OK;
}

static int32_t dadrastable(FGDATA *ff, FUNC *ftp)
{
  /*
    p3 = size
    p4 = GEN routine
    p5 = choose axis to use values of; 0 = x, 1 = y, 2 = z
    p6 = min output; if p6 and p7 == 0 don't scale
    p7 = max output; if p6 and p7 == 0 don't scale
    p8 = norm; 0 == normalize, -1 == don't normalize
    p9 = x start
    p10 = y start
    p11 = z start
    p12 = a constant value; 0 == default to 3.0
	p13 = b constant value; 0 == default to 2.7
	p14 = c constant value; 0 == default to 1.7
	p15 = d constant value; 0 == default to 2.0
	p16 = e constant value; 0 == default to 9.0
    p17 = time delta; 0 == default -> 0.15
    p18 = stepsize; 0 == default -> 1 
  */
  
  /* function table */
  MYFLT *fp = ftp->ftable;
  MYFLT a, b, c, d, e, time, step_size;
  
  /* default arguments */
  if (ff->e.p[12] == 0){
    a = 3.0;
  } else {
    a = ff->e.p[12];
  }
  if (ff->e.p[13] == 0){
    b = 2.7;
  } else {
    b = ff->e.p[13];
  }
  if (ff->e.p[14] == 0){
    c = 1.7;
  } else {
    c = ff->e.p[14];
  }
  if (ff->e.p[15] == 0){
    d = 2.0;
  } else {
    d = ff->e.p[15];
  }
  if (ff->e.p[16] == 0){
    e = 9.0;
  } else {
    e = ff->e.p[16];
  }  
  if (ff->e.p[17] == 0){
    time = 0.15;
  } else {
    time = ff->e.p[17];
  }  
  if (ff->e.p[18] == 0){
    step_size = 1;
  } else {
    step_size = ff->e.p[18];
  }
    
  /* p-field arguments */
  int32_t ft_size = ftp->flen;
  MYFLT min_out = ff->e.p[6];
  MYFLT max_out = ff->e.p[7];
  MYFLT norm = ff->e.p[8];
  MYFLT x = ff->e.p[9];
  MYFLT y = ff->e.p[10];
  MYFLT z = ff->e.p[11];
  MYFLT xx, yy;
  
  /* write raw data to ft from the dadras attractor
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
		/* fuctions for the dadras attractor		   
		 */
		xx = x + time * (y - (a * x) + (b * y * z));
		yy = y + time * ((c * y) - (x * y) + z);
		z =  z + time * ((d * x * y) - (e * z));
		x = xx;
		y = yy;
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
  { "dadras", dadrastable},
  { NULL, NULL}
};

FLINKAGE_BUILTIN(localfgens)
