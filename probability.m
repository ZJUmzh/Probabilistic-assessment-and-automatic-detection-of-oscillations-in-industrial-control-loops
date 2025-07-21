function P=probability(L,sigma_T_star,T_bar_star)

mu = sigma_T_star / T_bar_star;
variance_term = (sigma_T_star / T_bar_star)^2 * (0.5 + (sigma_T_star / T_bar_star)^2);
sigma = sqrt(1 / (L - 1) * variance_term);
z1 = -(mu / sigma);
z2 = ((1/3) - mu) / sigma;
P = normcdf(z2) - normcdf(z1);

end