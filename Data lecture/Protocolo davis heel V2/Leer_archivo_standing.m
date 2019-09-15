function [Data_standing Names] = Leer_archivo_standing(Pac) 

Pnum = num2str(Pac);
Nombre_Mar_Sta = ['Marcadores_Standing_P' Pnum '.txt'];

fid = fopen(Nombre_Mar_Sta,'rt');
 
oneline = fgets(fid);
oneline = strrep(oneline,'\n','');

i=1;
Data_standing = [];
Names = {};
IND = 1;

while ischar(oneline)
    
    %fprintf('%s',oneline)
    %fprintf('%s',i)
    
    if i == 11
        Split_Char = strsplit(oneline,'\t');
        for mm = 1: length(Split_Char)-1
                Names{mm} = Split_Char(mm);
        end
    end
  
    if i > 11
        
        Split_Char = strsplit(oneline);
        Data_num = [];
        
        if length(Split_Char) > 5
            for mm = 2: length(Split_Char)-1
                Data_num(mm-1) = str2double(Split_Char(mm));
            end
            
            F_nan = 0;
            for mm = 3: length(Data_num)
                if not(isnan(Data_num(mm)))
                    F_nan = 1;
                    break
                end     
            end
            
            if IND > 1 & F_nan == 1 
                if length(Data_num) ~= length(Data_standing(IND-1,:));
                    diff = Data_standing(IND-1,:) - length(Data_num);
                    for x = length(Data_num) + 1:length(Data_standing(IND-1,:))
                        Data_num(x) = NaN;
                    end
                end
            end
            
            
            if F_nan == 1;
                Data_standing(IND,:) = Data_num;
                IND = IND + 1;
            end
            
        end 
    end
     
    oneline = fgets(fid);
    oneline = strrep(oneline,'\n','');
    i=i+1;
end

end