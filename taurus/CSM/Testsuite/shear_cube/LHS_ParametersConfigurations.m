function LHS_ParametersConfigurations(filename, mu_train_Dimension)
%LHS_ParametersConfigurations generates a .vit file "filename.vit" 
%containing mu_train_Dimension parameter configurations
%
%   Usage example:
%
%   filename           = 'TestingParameters500'
%   mu_train_Dimension = 500


P = 3;

mu_min = [6*10^4 0.3 1000];

mu_max = [7*10^4 0.4 2000];

mu_cube            = lhsdesign(mu_train_Dimension, P); % normalized design
mu_train           = bsxfun(@plus,mu_min,bsxfun(@times,mu_cube,(mu_max-mu_min)));

fid = fopen([filename,'.vit'],'w+');
fprintf(fid, 'File with parameters (listed in this order per row: Young, Poisson, External Load)');
fprintf(fid, '\n%d %d', mu_train_Dimension, P);

for i = 1 : mu_train_Dimension
    
    fprintf(fid, '\n%f', mu_train(i,1));
    
    for j = 2 : P
        
        fprintf(fid, ' %f', mu_train(i,j));
        
    end
end

end


