__global__ void pixlight_analytic(float *pixlc, float *interpix, float *intrapix, int ntime, float* thetaX, float* thetaY){
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
  float pxr=px+spx-thetaY[0]+0.5/float(blockDim.y);
  float pyr=py+spy-thetaX[0]+0.5/float(blockDim.x);
  int k=0;
  for (int i=0; i<ntime; i++){
 
    /* pixel position vector from the PSF center */
    pxr=px+spx-thetaY[i]+0.5/float(blockDim.y);
    pyr=py+spy-thetaX[i]+0.5/float(blockDim.x);
    
    cache[thInd]=sensitivity*psf(pxr,pyr);

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
