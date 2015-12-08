function [mu sigma2] = estimateGaussian(X)
%ESTIMATEGAUSSIAN This function estimates the parameters of a 
%Gaussian distribution using the data in X
%   [mu sigma2] = estimateGaussian(X), 
%   The input X is the dataset with each n-dimensional data point in one row
%   The output is an n-dimensional vector mu, the mean of the data set
%   and the variances sigma^2, an n x 1 vector
% 

% Useful variables
[m, n] = size(X);

% You should return these values correctly
mu = zeros(n, 1);
sigma2 = zeros(n, 1);

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the mean of the data and the variances
%               In particular, mu(i) should contain the mean of
%               the data for the i-th feature and sigma2(i)
%               should contain variance of the i-th feature.
%
for j = 1:n
    mu_temp = 0;
    sigma2_temp = 0;
    for i = 1:m
        mu_temp = mu_temp + X(i,j);
    end
    mu_temp = mu_temp / m;
    for i = 1:m
        sigma2_temp = sigma2_temp + (X(i,j) - mu_temp)^2;
    end
    sigma2_temp = sigma2_temp / (m);
    
    mu(j) = mu_temp;
    sigma2(j) = sigma2_temp;
end

%for i=1:n,
%    mu(i) = sum(X(:,i))/m;
%    sigma2(i) = sum((X(:,i) .- mu(i)).^2)/m;
%end









% =============================================================


end
