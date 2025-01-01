/*!****************************************************************************                                                                         **/
/**   FILE         = CHANNELFLOW.F90                                        **/
/**   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 **/
/**   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        **/
/*****************************************************************************/
/**  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.           **/
/*****************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define nx 401
#define ny 101
int main(void){
     double u[nx][ny], un[nx][ny];
     double v[nx][ny], vn[nx][ny];
     double p[nx][ny], pn[nx][ny], F[nx][ny];
     double x[nx], y[ny];
     double xmax, ymax, dx, dy, rho, vis, dt, rms;
     int i, j, k, nt;
     xmax = 1.0;
     ymax = 0.25;
     dx = xmax/(nx-1.0);
     dy = ymax/(ny-1.0);
     nt = 10000, rho = 1.0, vis = 0.1, dt  =0.00001;
 
 // Grid generation
    x[0] = 0.0;
    for(i = 1; i < nx; i++)
      {
        x[i] = x[i-1]+dx;
      }
   
    y[0] = 0.0;
    for(i = 1; i < ny; i++)
      {
        y[i] = y[i-1]+dy;
      }

 // intial condition
    for(i = 0; i < nx; i++){
       for(j = 0; j < ny; j++){
         u[i][j]=0.0;
         v[i][j]=0.0;
         p[i][j]=0.0;
         F[i][j]=1.0;
         
        }
      }

 // computation
  
   for(k = 1; k<=nt; k++){
  
    
   // navier stokes
      
          for(j=1; j<(ny-1); j++){
             for(i=1; i<(nx-1); i++){
              un[i][j]= u[i][j]-u[i][j]*dt*(u[i][j]-u[i-1][j])/dx-v[i][j]*dt*(u[i][j]-u[i][j-1])/dy+vis*dt*(u[i+1][j]-2.0*u[i][j]+u[i-1][j])/(dx*dx)+ vis*dt*(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/(dy*dy);
              
              
              vn[i][j] = v[i][j]- u[i][j]*dt*(v[i][j]-v[i-1][j])/dx-v[i][j]*dt*(v[i][j]-v[i][j-1])/dy+vis*dt*(v[i+1][j]-2.0*v[i][j]+v[i-1][j])/(dx*dx)+vis*dt*(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/(dy*dy);
              }
              }


   // boundary condition
      for (i=0; i<nx; i++){
         un[i][0]=0.0;
         un[i][ny-1]=0.0;
    
         vn[i][0]=0.0;
         vn[i][ny-1]=0.0;
      }
     for(j=0; j<ny; j++){
          
          un[nx-1][j]=1.0;
          un[0][j]=1.0;
        
          
          vn[nx-1][j]=0.0;
          vn[0][j]=0.0;
       }


//pressure poisssion
        for(j=1; j<(ny-1); j++){
          for(i=1; i<(nx-1); i++){
            pn[i][j]=((p[i+1][j]+p[i-1][j])*dy*dy + (p[i][j+1]+p[i][j-1])*dx*dx
                     -((un[i+1][j]-un[i-1][j])/(2.0*dx)+(vn[i][j+1]-vn[i][j-1])/(2.0*dy))/dt*rho*dx*dx*dy*dy)/(2.0*(dx*dx + dy*dy));
        }
       }
 
           for (j = 0; j<ny; j++){
              pn[nx-1][j]=0.0;
              pn[0][j]=0.0;
        }
        
        for( i=0; i<nx; i++){
             pn[i][0]=pn[i][1];
             pn[i][ny-1]=pn[i][ny-2];
         }
//correcting velocity

     for(j=1; j<(ny-1); j++){
             for(i=1; i<(nx-1); i++){
              u[i][j]= un[i][j]-(pn[i+1][j]-pn[i-1][j])*dt/(2.0*dx*rho);
              v[i][j] = vn[i][j]-(pn[i][j+1]-pn[i][j-1])*dt/(2.0*dy*rho);
             }
           }
    
 // boundary condition
      for (i=0; i<nx; i++){
         u[i][0]=0.0;
         u[i][ny-1]=0.0;
    
         v[i][0]=0.0;
         v[i][ny-1]=0.0;
      }
     for(j=0; j<ny; j++){
          
          u[nx-1][j]=1.0;
          u[0][j]=1.0;
        
          
          v[nx-1][j]=0.0;
          v[0][j]=0.0;
       }
 //substituing back
          rms = 0.0;
         for(i = 0; i < nx; i++){
            for(j = 0; j < ny; j++){
                rms = rms + (p[i][j]-pn[i][j])*(p[i][j]-pn[i][j]);
              p[i][j]=pn[i][j];
             }
            }

         rms = sqrt(rms/(nx*ny));
         printf("Iteration  " "%d\t" " error  " "%E\n", k, rms);

       if (rms < 1.0e-06) {
          break;
     }
}
  
// output data
 FILE *f1;
 f1 = fopen("PUV.dat","w+t");


  
fprintf(f1," TITLE = \"FLOW-FIELD\"\n");
fprintf(f1," VARIABLES = X, Y, P, U, V\n");
fprintf(f1," ZONE T=  \"N= 0\", I= " "%d\t " ", J= " "%d\t" ",F=POINT\n", nx, ny);
for(j = 0; j < ny; j++){
  for(i = 0; i < nx; i++){           
            fprintf(f1,"%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", x[i], y[j], p[i][j], u[i][j], v[i][j]);
            }
          }
    fclose(f1);
}



