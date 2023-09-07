function [mae, mse, rmse] = calc_error(self)
% Mean absolute error
mae = mean(abs(self.pow - self.fooofed_spectrum));
% Mean squared error
mse = mean((self.pow - self.fooofed_spectrum).^2);
% Root mean squared error
rmse = sqrt(mean((self.pow - self.fooofed_spectrum).^2));

end