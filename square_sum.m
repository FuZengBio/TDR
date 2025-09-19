function [out]  = square_sum(data, span)
%NB don't change orientation/dimension of sum. Used in other applications.
idx=1;
while (idx * span) <= size(data,2)
    out(:,idx) = sum(data(:,(idx-1)*span + [1:span]),2);
    idx=idx + 1;
end
    