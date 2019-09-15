function Position= Find_name(celda, nombre)


n=0;
for x=1:length(celda)

    busqueda = strfind(celda{x},nombre);
    
    
    if  length(busqueda{1}) > 0    
        Position=x;
        n=1;
        break
    end
    
end

if n==0
fprintf ('The name "%s" you are searching for has not been found',nombre)
Position=NaN;
end

end 