% Code from Finding a sparse vector in a subspace:
% Linear sparsity using alternating directions
% by Qing Qu, Ju Sun, and John Wright
% code found here  https://sites.google.com/site/homeqingqu/miscellaneous
% solve the following problem 
% min_{q,x} 1/2*||Y*q - x||_2^2 + lambda * ||x||_1, s.t. ||q||_2 = 1 
% by alternating minimization method (ADM). 
% Y: input data
% q_init: initialization for q, lambda
% lambda: penalty parameter
% MaxIter: max iteration
% tol: tolerance for convergence
% q: output result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = ADM(Y,q_init,lambda,MaxIter,tol)
q = q_init;

for k = 1:MaxIter
    q_old = q;
    x = soft_thresholding(Y*q,lambda); % update y by soft thresholding
    q = Y'*x/norm(Y'*x,2); % update q by projection to the sphere
    res_q = norm(q_old-q,2);
    
    % this next section is to create a figure showing the algorithm at
    % work:
    out = Y*q; % compose sparsets vector
    
    
    if (res_q<=tol)
        
        return;
    end
    k;

    
    % if(mod(k,10)==0)
    %    res_obj = abs(obj-obj_old);
    %    fprintf('Running the %d-th iteration, diff_q=%f, diff_obj =%f \n',k, res_q,res_obj);
    % end
end
out = NaN.*Y*q

end

% soft-thresholding operator
function Y = soft_thresholding(X,d)
Y = sign(X).*max(abs(X)-d,0);
% Y = max(abs(X)-d,0);
end
