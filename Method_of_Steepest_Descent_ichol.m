function [ soln, niters ] = Method_of_Steepest_Descent_ichol( A, b, x0 )
    
    L_M = sparse(ichol(sparse(A), struct('type', 'ict','droptol',1e-3 ...
        ,'michol','off'))); 
    %L_M = sparse(ichol(sparse(A))); 
    L_M_transp = transpose(L_M);  %transpose of L_M
    r_i = b - A*x0; %initialize residual
    niters=0;
    x = x0;
    
    while(norm(r_i) >= eps(1)*(norm(b)))
        z = L_M \ r_i;
        p_i = L_M_transp \ z;
        q_i = A*p_i;
        scalar_alpha_i= dot(p_i, q_i); 
        alpha_i = (dot(p_i,r_i))/scalar_alpha_i;
        x = x + alpha_i*p_i;
        r_i = r_i - alpha_i * q_i;
        niters = niters+1;
    end
   soln = x;
end
