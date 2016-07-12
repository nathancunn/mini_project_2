load 'galaxy_perm.mat'

% Dirichlet process mixture models 
mu = 0;
sigmasq = 1;
a = 0.1;
M = 1;
[s1DP] = algorithm1_DP(galaxy, mu, sigmasq, a, M, 1000);
[s2DP] = algorithm2_DP(galaxy, mu, sigmasq, a, M,  200, 1000);
[s3DP] = algorithm3_DP(galaxy, mu, sigmasq, a, M, 1000);
[s4DP] = algorithm4_DP(galaxy, mu, sigmasq, a, M, 200, 1000);

% Dirichlet process mixture models with parameter updating
mu = 0;
sigmasq = 1;
M = 1;
[s3DPp, a3DPp] = algorithm3_DP_param(galaxy, mu, sigmasq, M, 1000);
[s3DPp, a3DPp] = algorithm4_DP_param(galaxy, mu, sigmasq, M, 200, 1000);

% Normalized generalized gamma process mixture models 
mu = 0;
sigmasq = 1;
a = 0.1;
gamma = 0.1;
M = 1;
[s3NGG] = algorithm3_NGG(galaxy, mu, sigmasq, a, gamma, M, 1000);
[s3NGG] = algorithm4_NGG(galaxy, mu, sigmasq, a, gamma, M, 200, 1000);

% Normalized generalized gamma process mixture models 
mu = 0;
sigmasq = 1;
gamma = 0.1;
M = 1;
[s3NGGp, a3NGGp] = algorithm3_NGG_param(galaxy, mu, sigmasq, gamma, M, 1000);
[s3NGGp, a3NGGp] = algorithm4_NGG_param(galaxy, mu, sigmasq, gamma, M, 200, 1000);
