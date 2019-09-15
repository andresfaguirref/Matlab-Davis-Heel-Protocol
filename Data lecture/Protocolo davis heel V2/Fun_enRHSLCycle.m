function porcentaje=Fun_enRHSLCycle( rhs1, rhs2, lhs1, lhs2)

long=lhs2-lhs1;

if rhs1>lhs1 && rhs1<lhs2
    
    y=rhs1-lhs1;
    
    x=y*(100)/long;
    
    valor=x;

elseif rhs2>lhs1 && rhs2<lhs2
    
    y=rhs2-lhs1;
    
    x=y*(100)/long;
    
    valor=x;
    
end 

porcentaje=valor;

end