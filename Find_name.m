function Name_pos = Find_name(names_cell,name)
    
    Y=0;
    for x=1:length(names_cell)
       
       find = strfind(names_cell{x},name);
       if length(find) > 0 
            Name_pos = x;
            Y=1;
            break
       end
       
    end
    
    if Y == 0
        fprintf('The name " %s " was not found \n', name);
        Name_pos = nan;
    end
    
end
