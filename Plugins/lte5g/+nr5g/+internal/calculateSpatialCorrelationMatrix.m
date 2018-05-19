%calculateSpatialCorrelationMatrix calculate spatial correlation matrix

%   Copyright 2017 The MathWorks, Inc.

function Rspat = calculateSpatialCorrelationMatrix(P,Rt,Rr,Gamma,a)

    % Compute overall spatial correlation matrix. Rspat is of size
    % (Nt*Nr)-by-(Nt*Nr)(-by-Np)
    Nt = size(P,1);
    Nr = size(P,2);
    Np = max([size(Rt,3) size(Rr,3) size(Gamma,3)]);
    if (Np>1)
        % If any of (Rt, Rr, Gamma) is 3-D i.e. defines a
        % matrix per path, repeat the other matrices for
        % uniformity
        if (size(Rt,3)==1)
            Rt = repmat(Rt,[1 1 Np]);
        end
        if (size(Rr,3)==1)
            Rr = repmat(Rr,[1 1 Np]);
        end
        if (size(Gamma,3)==1)
            Gamma = repmat(Gamma,[1 1 Np]);
        end
    end
    Rspat = zeros([Nt Nr Np]);
    for p = 1:Np
        Rspat(:,:,p) = P*kron(kron(Rt(:,:,p),Gamma(:,:,p)),Rr(:,:,p))*P.';
    end

    % Adjustment to make Rspat positive semi-definite (if necessary) as
    % defined in TS 36.101 / TS 36.104 Annex B
    if (a~=0)
        for p = 1:Np
            Rspat(:,:,p) = (Rspat(:,:,p) + a*eye(size(Rspat)))/(1+a);
        end
    end
    
end
