function punto_movido=Mover_punto(vector,punto,parametro)

norma=norm(vector);
vector_uni=vector/norma;

x=parametro*vector_uni(1)+punto(1);
y=parametro*vector_uni(2)+punto(2);
z=parametro*vector_uni(3)+punto(3);

punto_movido=[x,y,z];

end