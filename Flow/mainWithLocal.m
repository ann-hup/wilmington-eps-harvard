function structMat = mainWithLocal(num_rows,num_cols)
%MAINWITHLOCALS Summary of this function goes here
%   Detailed explanation goes here
structMat(1 : num_rows, 1 : num_cols) = struct;
    for r = 1 : num_rows
        for c = 1 : num_cols
            if r == c
                structMat(r, c).num = 1;
            else
                structMat(r, c).num = 0;
            end
            structMat(r, c).ids = [];
        end
    end
    
  %  structMat = modifyStruct(structMat, num_rows, num_cols);

end

function structMat = modifyStruct(structM, num_rows, num_cols)
    for r = 1 : num_rows
        for c = 1 : num_cols
            if mod(c, 2) == 0
                structM(r, c).ids = [structM(r, c).ids, 0, 2, 4, 6, 8];
            else
                structM(r, c).ids = [structM(r, c).ids, 1, 3, 5, 7, 9];
            end
        end
    end
    structMat = structM;
end
