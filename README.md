# Maxima-PDE-finite-differences
Some .mac files for calculating partial differential equations through the method of finite differences  

In 2013 an academic paper was published by the Turkish professor Tufan Sirin, which demonstrated how to do PDEs in Maxima using finite differences. The title of the paper is "Solutions of Partial Differential Equations with Computer Algebra." The math exposition is on a basic level; the whole paper is only four pages in length. Though I know jack about PDEs, it seems that the degree of accuracy in the described method could be increased as desired simply by finer discretization. 

I took the .mac files out of this paper and style-edited them somewhat, so they would run on M46. They have been found to run successfully in copy-paste mode from this repository.  

An attraction of this method is that it is free of dependencies; everything needed is available in the core of Maxima itself.
  
    

### The first is for the parabolic class of PDEs, such as the diffusion equation. The .mac file reads as follows:  
  
declare([r,iend,jend,fini],constant);  
 r:0.5;iend:20.0;jend:20;fini:100.0;  
 for j:1.0 thru jend do for i:1.0 thru iend do f[i,j]:fini;  
 for j:1.0 thru jend do f[1.0,j]:0.0;  
 for i:iend thru iend do for j:1.0 thru jend do  
f[iend,j]:0.0;  
 f[1.0,1.0]:fini/2.0;f[iend,1.0]:fini/2.0;  
 f[i,j]:for j:1.0 thru jend-1.0 do for i:2.0 thru iend-1.0 do  
(f[i,j+1.0]:r*(f[i+1.0,j]+f[i-1.0,j])+(1-2*r)*f[i,j]);  
 load(draw)$  
M:apply(matrix,makelist(makelist(f[i,j],j,1.0,19.0),i,1.0,19.0 ));$
draw3d (contour_levels={10, 20,30,40,50,60,70,80,90,100},  
contour = both,color = blue,elevation_grid(M,0,0,1,1),title = 
    "paraPDE.mac", xlabel  
= "x",ylabel = "y",zlabel ="Function value",surface_hide = true);  
im: apply(matrix,makelist(makelist (f[i,j],j,1.0,19.0),  
i,1.0,19.0)) $ 
/*draw2d(palette=gray,image(im,0,0,100,100))$ */
  

![paraPDE mac](https://user-images.githubusercontent.com/29483443/187046975-fef2f58b-3b2e-4da6-b41b-e02b6dc3419d.svg)






### The second is for the hyperbolic class of PDEs, such as the wave equation. The .mac file reads as follows:  

 declare([iend,jend,dx,dt],constant); iend:20; jend:20;  
 dx:.15; dt:.5;  
 for i from 1 thru 20 do for j from 1 thru 20 do (f[i,j]:0);  
 for i:1 thru iend do (xx[i]:dx*i);  
 for i:1 thru iend do (f[i,1]:sin(3.14*xx[i]));  
 for i:2 thru iend-1 do (f[i,2]:0.5*(f[i-1,1]+f[i+1,1]));  
 f[i,j]:for j:2 thru jend-1 do for i:2 thru iend-1 do  
(f[i,j+1]:(f[i+1,j]+f[i-1,j]-f[i,j-1])); 
load(draw);
M:apply(matrix,makelist(makelist(f[i,j+1],j,2,19),i,2,19))$
 for i from 1 thru iend do (y[i]:dx*i);   
 for j from 1 thru jend do (x[j]:dt*j);
load(draw); draw3d (contour_levels=  
{x[1],x[2],x[3], x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],  
x[12], x[13],x[14],x[15], x[16], x[17],x[18],x[19],x[20]},  
contour = both, color = blue, elevation_grid(M,0,0,8,24),  
title = "HypPDE.mac", xlabel = "Time", ylabel = "x", zlabel = "Function  value"); 



![hypPDE mac](https://user-images.githubusercontent.com/29483443/187046337-0094b575-94bd-4148-933b-4d6096b5e7aa.svg)

  

  
    
 ### The third and fourth are for the elliptical class of PDEs, such as Laplace's and Poisson's equations. The .mac files read as follows:  
   
     
 nx:10;ny:10;a:1.0;b:1.0;v1:10.0;v2:60.0;v3:0.0;v4:30.0;  
for i:1 thru nx-1 do for j:1 thru ny-1 
do(v[i,j]:(v1+v2+v3+v4)/4); 
fpprintprec:4$  
 for i:1 thru nx-1 do (v[i,1]:v1); 
for i:1 thru nx-1 
do(v[i,ny]:v3); for j:1 thru ny-1 do (v[1,j]:v4); 
for j:1 thru ny-1 do(v[nx,j]:v2);  
 for k:1 thru 23 do for i:2 thru nx-1 do for j:2 thru ny-1 do  
(v[i,j]:0.25*(v[i+1,j]+v[i-1,j]+v[i,j+1]+v[i,j-1])); 
 load(draw);  
M:apply(matrix,makelist(makelist(v[i,j],j,1,ny-1),  
i,1,nx)) $ 
draw3d (contour_levels = {10,15,20, 25,30  
,35, 40,45,50,60,70,80,90,100},contour = both,color =  
blue,elevation_grid(M,0,0,100,100),title = "lapPDE.mac",xlabel = "x",ylabel =  
"y",zlabel = "Function value",surface_hide = true); 
  
    
 ![lapPDEmac](https://user-images.githubusercontent.com/29483443/187046488-4c34fa0d-b03f-46ef-a8aa-0e351b4989e0.svg)
     
      
  
    
### And for Poisson's equation: 
  
    
      
      
   
   
 
nx:8;ny:8;a:1.0;h:a/nx;v1:10.0;v2:-30.0;v3:100.0;v4:-10.0; 
for i:1 thru nx-1 do for j:1 thru ny-1 do  
(v[i,j]:(v1+v2+v3+v4)/4);
for i:1 thru nx-1 do  
(v[i,1]:v1,v[i,ny]:v3) ; 
for j:1 thru ny-1 do  
(v[1,j]:v4,v[nx,j]:v2);
t:cos(3.14/nx)+cos(3.14/ny);  
 w: (8.0-sqrt(64.0-16.0 * t^2))/ (t^2);  
 w4:w/4.0; for k:1 thru 10 do for i:2 thru nx-1 do for j:2  
thru ny-1 do (v[i,j]:(x:h * i,y:h * j,g:-36.0 * 3.14 * x * (y-1.0),  
r:w4*(v[i+1,j]+v[i-1,j]+v[i,j+1]+v[i,j-1]-4* v[i,j]- 
g* h^2),v[i,j]:v[i,j]+(w4)* r));
load(draw);  
 M:apply (matrix,makelist(makelist (v[i,j],j,1,7),i,1,7))$
 draw3d (contour_levels = {10,20,30, 40,50, 60, 
70,80,90, 100}, contour= both,color = blue, elevation_grid 
(M,0,0,1,1), title = "poiPDE.mac", xlabel = "x",ylabel = "y",
zlabel = "Function value",surface_hide = true); 
  
    
      
      
![poiPDEmac](https://user-images.githubusercontent.com/29483443/187046619-1ca21454-094b-4cbb-ba82-5012399bfb60.svg)
   
     
 

  

