function porcentaje=Fun_enLHSRCycle( lhs1, lhs2, rhs1, rhs2)

long=rhs2-rhs1;

if lhs1>rhs1 && lhs1<rhs2
    
    y=lhs1-rhs1;
    
    x=y*(100)/long;
    
    valor=x;

elseif lhs2>rhs1 && lhs2<rhs2
    
    y=lhs2-rhs1;
    
    x=y*(100)/long;
    
    valor=x;
    
end 

porcentaje=valor;

end