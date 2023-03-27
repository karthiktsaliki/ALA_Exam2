function [ A ] = Create_Poisson_problem_A( N )
% Create the archtypical matrix A for an N x N Poisson problem (5-point
% stencil.
    A = zeros(N^2, N^2);
    
    for j = 1:(N*N)
        r=rem(j,N);
        % Set the diagonal
        A(j,j) = 4;
          
        % Set the entries of the first sub and super diagonals
        if j-1>=1 && (r~=1)  
            A(j,j-1) = -1;
        end 

        if j+1<=N*N && (r~=0)    
            A(j,j+1) = -1;
        end
        % Set the other off-diagonal entries
        if j-N >=1 
            A(j,j-N) = -1;
        end
        if j+N<=N*N 
            A(j,j+N) = -1;
        end 

    end

end
  
