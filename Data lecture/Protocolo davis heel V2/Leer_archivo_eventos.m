function Eventos = Leer_archivo_eventos(Pac) 


Pnum = num2str(Pac);
Nombre_Mar_Sta = ['Eventos_P' Pnum '.txt'];

fid = fopen(Nombre_Mar_Sta,'rt');
 
oneline = fgets(fid);
oneline = strrep(oneline,'\n','');

i=1;
Eventos = [];

eRHS = [];
eRTO = [];
eLHS = [];
eLTO = [];

IND = 1;
names = {};
while ischar(oneline)
    
    if i == 8
       Split_Char = strsplit(oneline,'\t');
       for mm=1:length(Split_Char)
        names{mm} = Split_Char(mm);
       end
       
       eRHS_pos = Find_name(names,'eRHS');
       eRTO_pos = Find_name(names,'eRTO');
       eLHS_pos = Find_name(names,'eLHS');
       eLTO_pos = Find_name(names,'eLTO');
    end
    
    if i > 8 
        Split_Char = strsplit(oneline,'\t');
        
        if length(Split_Char) > 3
            
            eRHS(IND) = str2double(Split_Char(eRHS_pos));
            eRTO(IND) = str2double(Split_Char(eRTO_pos));
            eLHS(IND) = str2double(Split_Char(eLHS_pos));
            eLTO(IND) = str2double(Split_Char(eLTO_pos));
     
%             Eventos(IND,1) = str2double(Split_Char(eRHS_pos));
%             Eventos(IND,2) = str2double(Split_Char(eRTO_pos));
%             Eventos(IND,3) = str2double(Split_Char(eLHS_pos));
%             Eventos(IND,4) = str2double(Split_Char(eLTO_pos));
            

            IND = IND + 1;
        end
    end
    
     
    oneline = fgets(fid);
    oneline = strrep(oneline,'\n','');
    i=i+1;
end

    C_1 = sort(eRHS);
    C_2 = sort(eRTO);
    C_3 = sort(eLHS);
    C_4 = sort(eLTO);
    
    
    Eventos(:,1) = C_1;
    Eventos(:,2) = C_2;
    Eventos(:,3) = C_3;
    Eventos(:,4) = C_4;
   

end

