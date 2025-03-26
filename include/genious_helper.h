#ifndef GENIOUS_HELPER_H
#define GENIOUS_HELPER_H

#include "./genious.h"

void scale_array(MYFLT *fp, MYFLT min_out, MYFLT max_out, int32_t size)
{
  MYFLT min_in = fp[0];
  MYFLT max_in = fp[0];
  
  for(int i = 1; i < size; i++){
    if(fp[i] < min_in)
      min_in = fp[i];
    if(fp[i] > max_in)
      max_in = fp[i];
  }

  MYFLT old_value;
  for(int j = 0; j < size; j++){
    old_value = fp[j];
    fp[j] =  min_out + ((old_value - min_in) / (max_in - min_in))
      * (max_out - min_out);
  }
}

#endif 
