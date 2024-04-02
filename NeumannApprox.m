function [Matrix_inv]=NeumannApprox(GramMatrix,flag)
% Input: Matrix_ori: the matrix due to be inversed (ChannelPhase'*ChannelPhase)
%        ChannelPhase:The part of Gram matrix
%        flag: the index of Neuman Series Approximation
% Output: Matrix_inv

 
OriDiag=diag(diag(GramMatrix));
OriOff=triu(GramMatrix,1)+tril(GramMatrix,-1); 
 
Matrix_inv = 0;
for i = 0:flag-1
    Matrix_inv = Matrix_inv+(( -inv(OriDiag)*OriOff)^i)*inv(OriDiag);
end

end


function [Matrix_inv]=TriNeumannApprox(GramMatrix,flag)
% Input: Matrix_ori: the matrix due to be inversed (ChannelPhase'*ChannelPhase)
%        ChannelPhase:The part of Gram matrix
%        flag: the index of Neuman Series Approximation
% Output: Matrix_inv

 
OriDiag=diag(diag(GramMatrix));
OriOff=triu(GramMatrix,1)+tril(GramMatrix,-1); 
 
Matrix_inv = 0;
for i = 0:flag-1
    Matrix_inv = Matrix_inv+(( -inv(OriDiag)*OriOff)^i)*inv(OriDiag);
end

end

 
