% Copyright (C) 2019 Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% All Rights Reserved.
%
% Authors: Albert Berahas, Frank E. Curtis, Baoyu Zhou
%
% Method definition for AggQN class

% Map a vector to a matrix
function [vec] = mapToVector(AQN,wholeMatrix)

vec = [];
for i = 1:size(wholeMatrix,2)
    vec = [vec ; wholeMatrix(i,1:i)'];
end

end