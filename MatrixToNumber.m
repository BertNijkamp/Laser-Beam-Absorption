function [number] = MatrixToNumber(matrix,decimals)
% This function converts the list particle numbers to a single number,
% in case a laser position is inside multiple particles.

Size = size(matrix,2);
number = 0;

for i = Size:-1:1
    number = number+matrix(Size-i+1)*10^(i*decimals-decimals);
end

end