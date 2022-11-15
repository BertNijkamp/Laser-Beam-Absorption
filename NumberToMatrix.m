function [matrix] = NumberToMatrix(number,decimals)
% This function converts a single number to a list of particle numbers,
% in case a laser position is inside multiple particles.

if number == 0 || number == 1
    matrix = number;
else
    Size = ceil(log10(number)/decimals);
    matrix = zeros(1,Size);

    for i = Size:-1:1
        spot = 10^(i*decimals-decimals);
        temp = floor(number/spot);
        matrix(1,Size-i+1) = temp;
        number = number-temp*spot;
    end
end

end