function M = compute_inverse(lambda,k,UU,SS,VV)
%COMPUTE_INVERSE Summary of this function goes here
%   Detailed explanation goes here
if k > 0
   S = eye(size(SS,1),k)*eye(k,size(SS,2))*diag(1./diag(SS));
   if lambda == 0 % Just compute the truncated solution
    M = VV*S*UU';
   else % Use truncated as initial value
	% M = VV*((SS^2 + lambda*eye(size(SS,2)))\(SS+ lambda*S)*UU');

     
    VVn = VV(:,1:k);
     UUn  = UU(:,1:k);
     SSn = SS(1:k,1:k);
     
    
   
      M = VVn*((SSn^2 + lambda*eye(size(SSn,2)))\SSn*UUn');
   end
else % Regular Tikhonov solution
    M = VV*((SS^2 + lambda*eye(size(SS,2)))\SS*UU');
    
end
end

