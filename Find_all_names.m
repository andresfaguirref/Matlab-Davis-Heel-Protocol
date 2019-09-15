function Name_pos = Find_all_names(names_cell,name)
    
    Name_pos = [];
    m = 1;
    for x=1:length(names_cell)
       
       find = strfind(names_cell{x},name);
       if length(find) > 0 
            Name_pos(m) = x;
            m=m+1;
       end
       
    end
    
    if length(Name_pos) == 0
        fprintf('The name " %s " was not found \n', name);
    end
    
end
