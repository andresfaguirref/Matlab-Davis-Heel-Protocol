function Medidas = Leer_archivo_medidas(Pac) 


Pnum = num2str(Pac);
Nombre_Mar_Sta = ['Medidas_P' Pnum '.txt'];

fid = fopen(Nombre_Mar_Sta,'rt');
 
oneline = fgets(fid);
oneline = strrep(oneline,'\n','');

i=1;
Medidas = [];

IND = 1;

while ischar(oneline)
    
    if i == 9
        Split_Char = strsplit(oneline,'\t');
        for mm = 1: length(Split_Char)-1
                Medidas(mm) = str2double(Split_Char(mm));
        end
    end
   
   
     
    oneline = fgets(fid);
    oneline = strrep(oneline,'\n','');
    i=i+1;
end

end