function z = f(x1,x2)
         x01 = 1; x02=2; tf=30; beta = 50;
         Qtf = [tf^3/3 tf^2/2; tf^2/2 tf];
         z = [x1 - x01 - tf*x02;x2-x02]'*Qtf*[x1 - x01 - tf*x02;x2-x02]-beta;
end