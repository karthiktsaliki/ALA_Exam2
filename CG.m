
function [ soln, niters ] = CG( A, b, x0 )
 
    curr_r = b-A*x0; %initialize residual
    niters=0;
    x = x0;
    
    while(norm(curr_r) >= eps*(norm(b)))
        dot_curr_r = dot(curr_r, curr_r);
        if niters==0
            curr_p = curr_r;
        else
            gamma_i = dot_curr_r/(dot(prev_r,prev_r));
            curr_p = curr_r + gamma_i * prev_p;
        end
    
        prev_p = curr_p;
        prev_r = curr_r;
        
        q_i = A*curr_p;
        scalar_alpha_i= dot(curr_p, q_i); 
        alpha_i = dot_curr_r/scalar_alpha_i;
        x = x + alpha_i*curr_p;
        curr_r = curr_r - alpha_i * q_i;
        niters = niters+1;

    end
   soln = x;

end 
