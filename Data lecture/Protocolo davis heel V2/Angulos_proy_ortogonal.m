function angulos=Angulos_proy_ortogonal(punto, vector_normal, vector_eje1, vector_eje2, punto_plano)

d=-1*(vector_normal(1)*punto_plano(1)+vector_normal(2)*punto_plano(2)+vector_normal(3)*punto_plano(3));
t=-1*((d+dot(vector_normal,punto))/(vector_normal(1)^2+vector_normal(2)^2+vector_normal(3)^2));

x=t*vector_normal(1)+punto(1);
y=t*vector_normal(2)+punto(2);
z=t*vector_normal(3)+punto(3);

punto_proy=[x, y, z];

punto_proy=punto_proy/norm(punto_proy);

angulos=acosd(dot(vector_eje1,punto_proy));


end
