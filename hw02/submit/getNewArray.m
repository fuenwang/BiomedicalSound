%
% EE6265 Fu-En Wang 106061531 HW2 11/14/2017
%

function [new_data] = getNewArray(origin, M, N)
    new_data = zeros(1, N);
    for i = 1:N
    if i * M <= 1000
        index = (i-1)*M+1 : i*M;
    else
        index = (i-1)*M+1 : length(origin);
    end
    new_data(i) = sum(origin(index));
    end
end