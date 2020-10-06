#define alpha_ejas 1.82 
#define beta_ejas 0.59  

#define PI 3.141592653589793
#define NORM 8.908387538030743

__device__ float BesselJ1(float x){

  double ax,z;
  double xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
         +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
      ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
         +y*(99447.43394+y*(376.9991397+y*1.0))));
      ans=ans1/ans2;
  } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
         +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }

  return ans;
  
}


__device__ float psf(float pxr, float pyr, float sigma2){
  float r=sqrt(pxr*pxr+pyr*pyr);
  double A=0.0;
  if(r < 10.0&&r>0){
  A=BesselJ1(PI*r/alpha_ejas)/(r/alpha_ejas) - BesselJ1(PI*r/beta_ejas)/(r/beta_ejas);
  }
  
  return float(A*A)/NORM;
}


/* Reference: https://www.atnf.csiro.au/computing/software/gipsx2/sub/bessel.c */

