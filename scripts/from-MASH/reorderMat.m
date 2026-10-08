function mat = reorderMat(mat0,id)
% mat = reorderMat(mat0,id)
%
% Order rows and columns of a square matrix according to input order.
%
% mat0: [N-by-N] original matrix
% id: [1-by-N] rows and columns new order
% mat: [N-by-N] ordered matrix

J = size(mat0,1);
if numel(id)~=J
    disp(['reorderMat: the number of indexes is different than the matrix',...
        ' size.'])
    mat = [];
    return
end
mat = zeros(size(mat0));
for j1 = 1:J
    for j2 = 1:size(mat0,2)
        if j2>J
            mat(j1,j2) = mat0(id(j1),j2);
        else
            mat(j1,j2) = mat0(id(j1),id(j2));
        end
    end
end
