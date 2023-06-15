/* 

cache usage:
0 -- NINTRA(blockDim.x*blockDim.y) - 1: subpixel value
blockDim.x*blockDim -- blockDim.x*blockDim + nsubtilex*nsubtiley - 1: PSF subtile

*/



__global__ void pixlight_custom(float *pixlc, float *interpix, float *intrapix, float *psfarr, int *subtilex, int *subtiley, int ntime, float* thetaX, float* thetaY){

  /* number of threads */
  /*  unsigned int nthread = blockDim.x*blockDim.y;*/
  float rnthread=float(NINTRA);

  /* subpixel (thread) positions */
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

  /* bilinear interporated value*/
  float bilin;

  /* checking subtile indices */
  int jy = subtilex[grInd];
  int jx = subtiley[grInd];

  /* thread cooperation read of subtile */
  unsigned int rat = int(NNSUBTILE/rnthread);
  int ipickx, ipicky, ipick;
  int i;  
  /* subtile=psfarr[jx:jx+Nsubtilex,jy:jy+Nsubtiley] */
  for (unsigned int m=0; m<rat+1; m++){
    i = m*NINTRA+thInd;
    ipicky = jy+i%NSUBTILEY;
    ipickx = jx+int(float(i)/float(NSUBTILEY));
    ipick =  ipicky + ipickx*PSFDIMY;
    
    if (i < NNSUBTILE){ 
      cache[i+NINTRA]=psfarr[ipick]; 
    }
  }
  __syncthreads();

  int k=0;
  float psfposx, psfposy;
  float x,x1,x2,y,y1,y2;
  float Q11,Q12,Q21,Q22;
  float F1,F2;
  int ic,icc;
  int nout;
  nout=0;
  for (int i=0; i<ntime; i++){

    /* pixel position vector from the PSF center */
    pxr=px+spx-thetaY[i]+0.5/float(blockDim.y);
    pyr=py+spy-thetaX[i]+0.5/float(blockDim.x);

    psfposx=PSFCENTERX + pxr/PSFSCALE;
    psfposy=PSFCENTERY + pyr/PSFSCALE;
    x = psfposx - float(jx);
    y = psfposy - float(jy);

    x1 = int(x);
    x2 = x1 + 1;
    y1 = int(y);
    y2 = y1 + 1;

    /* Use shared memory */
    icc=(y1+x1*NSUBTILEY+NINTRA);
    Q11=cache[icc];
    icc=(y2+x1*NSUBTILEY+NINTRA);
    Q12=cache[icc];
    icc=(y1+x2*NSUBTILEY+NINTRA);
    Q21=cache[icc];
    icc=(y2+x2*NSUBTILEY+NINTRA);
    Q22=cache[icc];          

    /* Use global memory */
    /*       
    Q11=psfarr[int(psfposy)+int(psfposx)*PSFDIMY];
    Q12=psfarr[int(psfposy+1)+int(psfposx)*PSFDIMY];
    Q21=psfarr[int(psfposy)+int(psfposx+1)*PSFDIMY];
    Q22=psfarr[int(psfposy+1)+int(psfposx+1)*PSFDIMY];
    */
    
    F1=(float(x2)-x)*Q11 + (x-float(x1))*Q21;
    F2=(float(x2)-x)*Q12 + (x-float(x1))*Q22;
    bilin=(float(y2)-y)*F1 + (y-float(y1))*F2;

    cache[thInd]=sensitivity*bilin;

    /* reduce the tail effect (read Issue #28)*/
    /*
    if(pxr*pxr+pyr*pyr>100.0){
      cache[thInd]=0.0;
    }
    */
    
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
      pixlc[i+grInd*ntime]=cache[0]/float(ntime)/float(blockDim.x*blockDim.y);
    }
    
    __syncthreads();
  }
}
