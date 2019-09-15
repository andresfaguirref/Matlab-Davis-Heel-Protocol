function centro_masa=CoM(prox_RC, distal_RC, parametro)

%Se utiliza los tablas de antropometrica para sacar los centros de masa
vector=distal_RC-prox_RC;

for x=1:length(prox_RC(:,1))
   
centro(x,:)=Mover_punto(vector(x,:),prox_RC(x,:), parametro*norm(vector(x,:))); 
    
end

centro_masa = centro;

end