function error_of_square=er(y)
ysum_time_average=sum(y)/length(y);
ysum_time_average_square=ysum_time_average^2;
ysum_square_time_average=sum(y.^2)/length(y);
error_of_square=ysum_square_time_average-ysum_time_average_square;