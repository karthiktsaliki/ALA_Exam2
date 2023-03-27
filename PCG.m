function [ soln, niters ]  = PCG( A, b, x0 )
   
    L_M = sparse(ichol(sparse(A), struct('type', 'ict','droptol',1e-3,'michol','off'))); 
    %L_M = sparse(ichol(sparse(A))); 
    L_M_transp = transpose(L_M);
    curr_r = b-A*x0; %initialize residual
    niters=0;
    x = x0;
    
    while(norm(curr_r) >= eps(1)*(norm(b)))
       
        t = L_M \ curr_r;
        curr_z = L_M_transp \ t;

        dot_curr_r = dot(curr_r, curr_z);
        if niters==0
            curr_p = curr_z;
        else
            gamma_i = dot_curr_r/(dot(prev_r,prev_z));
            curr_p = curr_z + gamma_i * prev_p;
        end
    
        prev_p = curr_p;
        prev_r = curr_r;
        prev_z = curr_z;
        
        q_i = A*curr_p;
        scalar_alpha_i= dot(curr_p, q_i); 
        alpha_i = dot_curr_r/scalar_alpha_i;
        x = x + alpha_i*curr_p;
        curr_r = curr_r - alpha_i * q_i;
        niters = niters+1;


    end
   soln = x;
  
end
