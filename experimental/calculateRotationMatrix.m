function R = calculateRotationMatrix(fromVector, toVector)
% This implementation follows Eq. (4.1.) from 
%   Cid, Jose Ángel; Tojo, F. Adrián F. "A Lipschitz condition along a 
%   transversal foliation implies local uniqueness for ODEs". Electronic 
%   Journal of Qualitative Theory of Differential Equations. 13: 1-14.
%   doi:    https://doi.org/10.14232/ejqtde.2018.1.13
%
%
% See also: https://math.stackexchange.com/a/476311
%           https://en.wikipedia.org/wiki/Rotation_matrix#Vector_to_vector_formulation
%
% (c) Georg Götz, Aalto University, 2020

if norm(fromVector - toVector) < 1e-10
    % when both vectors are identical, no rotation is required
    R = eye(3);
elseif fromVector == - toVector
    if all(fromVector == [1, 0, 0])
        % rotate 180 degrees around z-axis
        R = [cos(pi), -sin(pi), 0; ...
             sin(pi), cos(pi), 0; ...
             0, 0, 1];
    else
        error(['A rotation matrix from a vector to its opposing vector '...
            'cannot be calculated. Check around which axes the vector '...
            'has to be rotated and do the rotation manually.']);
    end
else
    % arbitrary rotation, see above references
    v = cross(fromVector, toVector);
    c = dot(fromVector, toVector);
    vx = [0, -v(3), v(2);...
          v(3), 0, -v(1);...
          -v(2), v(1), 0];

    R = eye(3) + vx + vx*vx*(1/(1+c));
    R = R';
end

end