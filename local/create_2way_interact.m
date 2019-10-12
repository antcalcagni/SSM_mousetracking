function D = create_2way_interact(D1,D2,type)

switch type
    case 'categorical'        
        D=[];
        for i=1:size(D1,2)
            D = [D cell2mat(arrayfun(@(j)min(D1(:,i),D2(:,j)),1:size(D2,2),'UniformOutput',false))]; 
        end
    case 'continuous'
        D = repmat(D1,1,size(D2,2)).*D2;
end
        
end