%calculateGammaMatrix calculate gamma matrix from polarization angles and XPR(s)

%   Copyright 2017 The MathWorks, Inc.

% IEEE 802.16m-08/004r5 Appendix B
function Gamma = calculateGammaMatrix(txPolAngles,rxPolAngles,XPRdB)
                    
    % Transmit polarization matrix
    Ptx = [cosd(txPolAngles); sind(txPolAngles)];

    % Receive polarization matrix
    Prx = [cosd(rxPolAngles); sind(rxPolAngles)];

    Npt = numel(txPolAngles);
    Npr = numel(rxPolAngles);
    Np = numel(XPRdB);
    Gamma = zeros(Npt*Npr,Npt*Npr,Np);
    for p = 1:Np
        Gamma(:,:,p) = calculateGammaMatrixOnePath(Ptx,Prx,XPRdB(p));
    end
    
end

% Note that the implementation here is formulated for an uplink case
% derived from IEEE 802.16m-08/004r5 Appendix B i.e. P_BS = Prx and P_MS =
% Ptx. This yields a different polarization matrix structure to IEEE
% 802.16m-08/004r5 Appendix B as written (the document is written for a
% downlink case) but results in the correct polarization matrix structure
% for 3GPP TS 36.101 / TS 36.104
function Gamma = calculateGammaMatrixOnePath(Ptx,Prx,XPRdB)

    % Linear cross-polar ratio (XPR) value
    XPR_linear = 10^(XPRdB/10);

    % Horizontal and vertical cross-polar ratios
    XPR_V = XPR_linear;
    XPR_H = XPR_linear;

    % Mean power of polarization components defined in terms of
    % XPRs:
    % XPR_V = p_vh / p_vv
    % XPR_H = p_hv / p_hh
    p_vv = 1;
    p_hh = 1;
    p_hv = XPR_H*p_hh;
    p_vh = XPR_V*p_vv;

    % Polarization matrix components. These matrix components
    % use 4-element vectors in conjunction with cell array
    % multiplication operation 'mult' (implemented below) in
    % order to satisfy the property of uncorrelated fading
    % between different elements in the channel polarization
    % matrix S:
    % E{ s_ij * conj(s_kl) } = 0, i~=k, j~=l
    s_vv = sqrt([p_vv 0 0 0]);
    s_vh = sqrt([0 p_vh 0 0]);
    s_hv = sqrt([0 0 p_hv 0]);
    s_hh = sqrt([0 0 0 p_hh]);

    % Channel polarization correlation matrix S
    S = {s_vv s_vh; s_hv s_hh};

    % Total polarization channel consisting of receive
    % polarization, channel polarization and transmit
    % polarization. This uses the cell array multiplication
    % function 'mult' to keep uncorrelated components separate
    Q = mult(mult(Prx.',S),Ptx);

    % Compute the polarization covariance matrix. This uses the
    % cell array multiplication function 'mult' to keep
    % uncorrelated components separate and finally the 'sum'
    % operation combines any remaining correlated terms
    Gamma = arrayfun(@(x)cellfun(@sum,x),mult({Q{:}}.',{Q{:}}));

    % Normalise the polarization covariance matrix for
    % unit power per antenna
    d = diag(Gamma);
    Gamma = Gamma*numel(d)/sum(d);

    function z = mult(x,y)

        [xr,xc] = size(x);
        if (~iscell(x))
            x = num2cell(x);
        end

        [~,yc] = size(y);
        if (~iscell(y))
            y = num2cell(y);
        end

        z = cell(xr,yc);

        for i = 1:xr
             for j = 1:yc
                 z{i,j} = 0;
                 for k = 1:xc
                     z{i,j} = z{i,j} + x{i,k}.*y{k,j};
                 end
             end
        end

    end

end
