function Interpolacion=Interpolacion_cubica(x,y,p)

x_new=[];
y_new=[];
n=1;

for m=1:length(x)
   if not (isnan(y(m)))
   
       x_new(n)=x(m);
       y_new(n)=y(m);
       
       n=n+1;       
   end
   
end

Interpolacion=spline(x_new,y_new,x); 

if p==1
    figure();
    plot(x_new,y_new, '*r');
    
    hold on
    plot (x,Interpolacion, 'o');

end 