__device__ float psf(float pxr, float pyr, float sigma2){
  float r2=pxr*pxr+pyr*pyr;
  float val=exp(-r2/2.0/sigma2)/(2.0*PI*sigma2);
  return val;
}
