function Eje_global_marcha=Eje_Global_marcha(names,static)

Eje_Global=[1 0 0; 0 1 0; 0 0 1];

X=Find_name(names,' r heel.X');
Y=Find_name(names,' r heel.Y');
Z=Find_name(names,' r heel.Z');

R_Heel=[static(:,X), static(:,Y), static(:,Z)];

X=Find_name(names,' l heel.X');
Y=Find_name(names,' l heel.Y');
Z=Find_name(names,' l heel.Z');

L_Heel=[static(:,X), static(:,Y), static(:,Z)];

vector=L_Heel-R_Heel;

for x=1:length(L_Heel(:,1)) 
    
    distancia=norm(vector(x,:));
    vector_unitario(x,:)=vector(x,:)/distancia;
    
end

Media=mean(vector_unitario);

vGaitJ=cross(Media, Eje_Global(2,:));

vGaitK=cross(Eje_Global(2,:), vGaitJ);

Eje_global_marcha=[Eje_Global(2,:); vGaitJ; vGaitK];


end