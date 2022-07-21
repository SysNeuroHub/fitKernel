function P = dirDotProdTest(direction_rad, respAmp, baseAmp)
%P = dirDotProdTest(direction_rad, respAmp)
% tests if the respAmp is deviated from (0,) in 2D space
%
%P = dirDotProdTest(direction_rad, respAmp, baseAmp)
% tests if the respAmp and baseAmp is from the same distribution in 2D space
%
%Mazurek, Kager, Van Hooser 2014 frontiers fig7,9

if nargin < 3
    baseAmp = [];
end

dot_p = respAmp .* exp(1i*direction_rad);
X_p = [real(dot_p), imag(dot_p)];

if ~isempty(baseAmp)
    dot_b = baseAmp .* exp(1i*direction_rad);
    X_b = [real(dot_b), imag(dot_b)];
    [~, P] = T2Hot2d([X_p; X_b]);
else
    [~, P] = T2Hot1(X_p);
end
