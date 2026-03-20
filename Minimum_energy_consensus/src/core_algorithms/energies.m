xc = -0.5 ;  
yc = 0;
IC = [3,-1.1; 
      -3, 1.1;  
      3.1, -3;
      -3.1 ,3; 
       6.5, 6.5 ;
       -7.5,-6.5];
tf = 8.158576;



Q = [tf^3/3, tf^2/2; tf^2/2, tf];                               
quad = @(dx) sum((Q \ dx).*dx,1);  

f1 = @(x,y) quad([x - IC(1,1) - tf*IC(1,2); y - IC(1,2)]);
f2 = @(x,y) quad([x - IC(2,1) - tf*IC(2,2); y - IC(2,2)]);
f3 = @(x,y) quad([x - IC(3,1) - tf*IC(3,2); y - IC(3,2)]);
f4 = @(x,y) quad([x - IC(4,1) - tf*IC(4,2); y - IC(4,2)]);
f5 = @(x,y) quad([x - IC(5,1) - tf*IC(5,2); y - IC(5,2)]);
f6 = @(x,y) quad([x - IC(6,1) - tf*IC(6,2); y - IC(6,2)]);

[f1(xc,yc), f2(xc,yc),f3(xc,yc),f4(xc,yc),f5(xc,yc),f6(xc,yc)]