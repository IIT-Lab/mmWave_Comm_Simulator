%h5gCRCDecode Cyclic redundancy check decoding and removal
%   [BLK,ERR] = h5gCRCDecode(...) checks the input data vector for a CRC
%   error assuming the vector comprises a block of data with the associated
%   CRC bits attached. The data part of the input is returned in vector
%   BLK. The logical difference (xor) between the attached CRC and the CRC
%   recalculated across the data part of the input is returned in uint32
%   scalar ERR. If ERR ~= 0 then either an error has occurred or the input
%   CRC has been masked. A logical mask can also be applied directly to
%   ERR.
%
%   [BLK,ERR] = h5gCRCDecode(BLKCRC,POLY) returns BLK, the data only part
%   of the combined data and CRC input vector BLKCRC, and uint32 ERR, the
%   logical (xor) CRC difference. The CRC polynomial is defined by a value
%   from the set ('6','11','16','24A','24B','24C'). See TS 38.212 Section
%   5.1 for the associated polynomials.
%
%   [BLK,ERR] = h5gCRCDecode(BLKCRC,POLY,MASK) behaves as above except the
%   CRC difference is also XORed with the scalar MASK parameter before it
%   is returned in ERR.
%
%   Example:
%   % This example shows the effect of CRC decoding a block of data with 
%   % and without a mask. When decoding without a mask, ERR is equal to the
%   % RNTI because the CRC had been masked during coding. Hence, the 
%   % logical difference between the original CRC and the recalculated CRC
%   % is the CRC mask. CRC decoding using the RNTI as a mask results in 
%   % zero ERR.
%   
%   % CRC encode an all-ones vector and mask with RNTI
%   RNTI = 8;
%   blkCrc = h5gCRCEncode(ones(100,1),'24A',RNTI);
%   % CRC decode the result without using a mask, ERR should be the RNTI
%   [blk,err1] = h5gCRCDecode(blkCrc,'24A');
%   err1
%   % CRC decoding using the RNTI as a mask, ERR should be zero
%   [blk,err2] = h5gCRCDecode(blkCrc,'24A',RNTI);
%   err2
%   
%   See also h5gCRCEncode, h5gCodeBlockDesegment, lteCRCDecode.

%   Copyright 2018 The MathWorks, Inc.

function [blk,err] = h5gCRCDecode(blkcrc, poly, mask)

    narginchk(2, 3);
    if nargin == 2
        mask = 0;
    end
    
    % Validate CRC type and obtain CRC parity bit length
    polylist = {'6', '11', '16', '24A', '24B', '24C'};
    crclenlist = {6, 11, 16, 24, 24, 24};
    encoderInd = strcmpi(poly,polylist);
    if ~any(encoderInd)
        error('lte:error',['The input CRC polynormial type is not one of the set (' strjoin(polylist,', ') ')']);
    end
    crcBitsLen = crclenlist{encoderInd};
    
    % Validate blkcrc
    if ~isvector(blkcrc) && ~isempty(blkcrc)
        error('lte:error','The input must be a vector and not a matrix.');
    end
    
    % Validate inputs and prepare the outputs
    reEncodedBlkcrc = h5gCRCEncode(blkcrc(1:end-crcBitsLen), poly, mask);
    
    % For XOR operation later
    blkcrc = blkcrc(:);
    
    % Calculate CRC error and prepare final output
    if isempty(blkcrc)
        blk = zeros(0,1,class(blkcrc));
        err = zeros(0,1,'uint32');
    else
        blk = reEncodedBlkcrc(1:end-crcBitsLen);
        % Calculate CRC error
        if length(blkcrc)<=crcBitsLen % For input length less than parity bit length
            if isempty(mask)
                mask = 0;
            end
            % uint32 conversion is to handle inf or nan inputs,
            % double conversion is for bitshift opertion later
            mask = double(uint32(real(mask(1))));
            
            blkcrc = [zeros(crcBitsLen-length(blkcrc),1); blkcrc];
            maskBits = mod(bitshift(mask,(1-crcBitsLen:0)'),2);
            errBits = xor(maskBits,blkcrc);
        else
            errBits = xor(reEncodedBlkcrc(end-crcBitsLen+1:end),blkcrc(end-crcBitsLen+1:end));
        end
        err = uint32(sum(errBits.*(2.^(crcBitsLen-1:-1:0)')));
    end

end