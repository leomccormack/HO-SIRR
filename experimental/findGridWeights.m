function [weights, order] = findGridWeights(azi, zen, order)
%findGridWeights Find / approximate quadrature weights by pseudo-inverse.
%   Grid azi and zen in rad, order optional (searches max if not given).
% CHris Hold 2020, 2022

assert(all(size(azi) == size(zen)))
assert(size(azi, 2) == 1)

if nargin < 3
    % find max grid order by condition number, this might take a while
    for itn = 1:100
        Y_it = getSHreal(itn, [azi, zen]);
        cn = cond(Y_it' * Y_it);
        % check condition number
        if cn > 2*(itn+1)  % experimental condition.
            order = itn - 1;
            break
        end
    end
end

Y_N = getSHreal(order, [azi, zen]);
%W = pinv(Y_N') * pinv(Y_N);
%W = ((Y_N.' * Y_N) \ Y_N.').' * ((Y_N.' * Y_N) \ Y_N.');  % faster
%weights = sum(W, 2);

% Fornberg, B., & Martel, J. M. (2014). 
% On spherical harmonics based numerical quadrature over the surface of a 
% sphere. 
% Advances in Computational Mathematics, 40(5–6)
%P_leftinv = pinv(Y_N);
P_leftinv = (Y_N.' * Y_N) \ Y_N.';
weights = sqrt(4*pi) * P_leftinv(1, :).';

if abs(sum(weights) - 4*pi) > 0.01 || any(weights < 0 )
    warning('Could not calculate weights')
end

end



function Y_N = getSHreal(N, dirs)
    % getSHreal adapted from 
    % https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/getSH.m
    Ndirs = size(dirs, 1);
    Nharm = (N+1)^2;

    Y_N = zeros(Nharm, Ndirs);
    idx_Y = 0;
    for n=0:N
        
        m = (0:n)';
            Lnm_real = legendre(n, cos(dirs(:,2)'));
            if n~=0
                condon = (-1).^[m(end:-1:2);m] * ones(1,Ndirs);
                Lnm_real = condon .* [Lnm_real(end:-1:2, :); Lnm_real];
            end
            norm_real = sqrt( (2*n+1)*factorial(n-m) ./ (4*pi*factorial(n+m)) );
            
            Nnm_real = norm_real * ones(1,Ndirs);
            if n~=0
                Nnm_real = [Nnm_real(end:-1:2, :); Nnm_real];
            end            
            
            CosSin = zeros(2*n+1,Ndirs);
            CosSin(n+1,:) = ones(1,size(dirs,1));
            if n~=0
                CosSin(m(2:end)+n+1,:) = sqrt(2)*cos(m(2:end)*dirs(:,1)');
                CosSin(-m(end:-1:2)+n+1,:) = sqrt(2)*sin(m(end:-1:2)*dirs(:,1)');
            end
            Ynm = Nnm_real .* Lnm_real .* CosSin;
            
        
        Y_N(idx_Y+1:idx_Y+(2*n+1), :) = Ynm;
        idx_Y = idx_Y + 2*n+1;
    end
    Y_N = Y_N.';
end
