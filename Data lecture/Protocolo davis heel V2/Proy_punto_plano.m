function punto_proy = Proy_punto_plano(punto, vector_normal, punto_plano)

d=-1*(vector_normal(1)*punto_plano(1)+vector_normal(2)*punto_plano(2)+vector_normal(3)*punto_plano(3));
t=-1*((d+dot(vector_normal,punto))/(vector_normal(1)^2+vector_normal(2)^2+vector_normal(3)^2));

x=t*vector_normal(1)+punto(1);
y=t*vector_normal(2)+punto(2);
z=t*vector_normal(3)+punto(3);

punto_proy=[x, y, z];

end