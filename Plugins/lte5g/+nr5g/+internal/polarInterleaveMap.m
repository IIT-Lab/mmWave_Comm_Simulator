function pi = polarInterleaveMap(K)
% Interleaver mapping pattern for Polar coding
%
%   pi = nr5g.internal.polarInterleaveMap(K) returns the interleaving
%   pattern for a length K.
% 
%   See also h5gPolarEncoder, h5gPolarDecoder.

%   Copyright 2018 The MathWorks, Inc.

%#codegen

%   Reference:
%   [1] 3GPP TS 38.212, "3rd Generation Partnership Project; Technical 
%   Specification Group Radio Access Network; NR; Multiplexing and channel 
%   coding (Release 15), v15.0.0, 2017-12. Section 5.3.1.1.

Kilmax = 164;                         
pat = getPattern();
pi = zeros(K,1);
k = 0;
for m = 0:Kilmax-1
    if pat(m+1) >= Kilmax-K
        pi(k+1) = pat(m+1)-(Kilmax-K);
        k = k+1;
    end
end

end

%-------------------------------------------------------------------------
function pat = getPattern()

% Table 5.3.1.1-1: Interleaving pattern for Kilmax = 164
pat = [ 0; 2; 4; 7; 9; 14; 19; 20; 24; 25; 26; 28; 31; 34; 42; 45; 49; ...
        50; 51; 53; 54; 56; 58; 59; 61; 62; 65; 66; 67; 69; 70; 71; 72; ...
        76; 77; 81; 82; 83; 87; 88; 89; 91; 93; 95; 98; 101; 104; 106; ...
        108; 110; 111; 113; 115; 118; 119; 120; 122; 123; 126; 127; 129; ...
        132; 134; 138; 139; 140; 1; 3; 5; 8; 10; 15; 21; 27; 29; 32; 35; ...
        43; 46; 52; 55; 57; 60; 63; 68; 73; 78; 84; 90; 92; 94; 96; 99; ...
        102; 105; 107; 109; 112; 114; 116; 121; 124; 128; 130; 133; 135; ...
        141; 6; 11; 16; 22; 30; 33; 36; 44; 47; 64; 74; 79; 85; 97; 100; ...
        103; 117; 125; 131; 136; 142; 12; 17; 23; 37; 48; 75; 80; 86; ...
        137; 143; 13; 18; 38; 144; 39; 145; 40; 146; 41; 147; 148; 149; ...
        150; 151; 152; 153; 154; 155; 156; 157; 158; 159; 160; 161; 162; ...
        163 ];
    
end
