__global__ void pixlight(float *pixlc, float *interpix, float *intrapix, int ntime, float* thetaX, float* thetaY, float sigma2){
  /* subpixel (thread) positinos */
  float spy = float(threadIdx.x)/float(blockDim.x);
  float spx = float(threadIdx.y)/float(blockDim.y); 
  int thInd = threadIdx.y + threadIdx.x*blockDim.y;
  
  /* pixel (grid) positions */
  float py = float(blockIdx.x);
  float px = float(blockIdx.y); 
  int grInd = blockIdx.y + blockIdx.x*gridDim.y;
  
  /* map positions */
  /* int y = threadIdx.x + blockIdx.x*blockDim.x; */
  /* int x = threadIdx.y + blockIdx.y*blockDim.y; */
  /* int offset = y + x * blockDim.y * gridDim.y; */
  
  /* total detector sensitivity */
  float sensitivity = interpix[grInd]*intrapix[thInd];
  
  /* pixel position vector from the PSF center */
  float pxr=px+spx-thetaY[0];
  float pyr=py+spy-thetaX[0];

  /* PSF Gaussian sigma2 */
  /* float sigma2=2.0; */
  int k=0;
  for (int i=0; i<ntime; i++){
    pxr=px+spx-thetaY[i];
    pyr=py+spy-thetaX[i];
    
    cache[thInd]=sensitivity*psf(pxr,pyr,sigma2);

    /*cache[thInd]=sensitivity;*/
    __syncthreads();
    
    /* thread adding */
    k = blockDim.x*blockDim.y/2;
    while (k !=0){
      if(thInd < k){
	cache[thInd] += cache[thInd+k];
      }
      __syncthreads();
      k /= 2;
    }
    if (thInd == 0){
      pixlc[i+grInd*ntime]=cache[0];
    }
    
    __syncthreads();
  }
}
