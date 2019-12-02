% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Map a vector to a matrix
function [wholeMatrix] = mapToMatrix(AQN,vec,width)

wholeMatrix = reshape(vec,[width+1,width]);

% count = 0;
% 
% for i = 1:size(wholeMatrix,2)
%     wholeMatrix(i+1:end,i) = vec(count+1:count+size(wholeMatrix,2)+1-i);
%     count = count + size(wholeMatrix,2)+1-i;
% end

end