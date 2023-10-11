function [D, S, Error] = compute_metrics(P, Q, n)

% MA Tschopp
% Purpose: Implementation of Cha (2007) similarity/distance metrics
% For test purposes:
% P = rand(1, 10); P = P / sum(P);
% Q = rand(1, 10); Q = Q / sum(Q);
% P(1) = 0;
% Q(2) = 0;
% Q(3) = P(3);
% P = P / sum(P);
% Q = Q / sum(Q);
% Q(3) = P(3);
% n = 3;
% [D1, S1, ~] = compute_metrics(P, Q, n)
% P(P == 0) = 1e-20;
% Q(Q == 0) = 1e-20;
% P(P == Q) = P(P == Q) + 1e-20;
% [D, S, ~] = compute_metrics(P, Q, n)

% Correct for divide by zero

P1 = P;
Q1 = Q;
P1(P == 0) = 1e-20;
Q1(Q == 0) = 1e-20;
P1(P == Q) = P1(P == Q) + 1e-20;

D = [];
S = [];
Error = [];

% Lp Minkowski family
D.euclidean = sqrt(sum((P - Q).^2)); % sqrt(dot(P-Q,P-Q));
D.cityblock = sum(abs(P - Q));
D.minkowski = sum(abs(P - Q).^n)^(1/n);
D.chebyshev = max(abs(P - Q));

% L1 family
D.sorensen = sum(abs(P - Q)) / sum(P + Q);
D.gower = sum(abs(P - Q)) / length(P);
D.soergel = sum(abs(P - Q)) / sum(max(P, Q));

if sum(min(P, Q)) ~= 0
    D.kulczynski_d = sum(abs(P - Q)) / sum(min(P, Q));
else
    D.kulczynski_d = sum(abs(P1 - Q1)) / sum(min(P1, Q1));
    Error.kulczynski_d = 'Divide by Zero Error';
end

if min(P + Q) ~= 0
    D.canberra = sum(abs(P - Q) ./ (P + Q));
else
    D.canberra = sum(abs(P1 - Q1) ./ (P1 + Q1));
    Error.canberra = 'Divide by Zero Error';
end

D.lorentzian = sum(log(1 + abs(P - Q)));

% Intersection family
S.intersection = sum(min(P, Q));
D.intersection = 1 - sum(min(P, Q));

if min(max(P, Q)) ~= 0
    D.wavehedges = sum(1 - min(P, Q) ./ max(P, Q));
else
    D.wavehedges = sum(1 - min(P1, Q1) ./ max(P1, Q1));
    Error.wavehedges = 'Divide by Zero Error';
end

S.czekanowski = 2 * sum(min(P, Q)) / sum(P + Q);
D.czekanowski = 1 - 2 * sum(min(P, Q)) / sum(P + Q);
S.motyka = sum(min(P, Q)) / sum(P + Q);
D.motyka = 1 - sum(min(P, Q)) / sum(P + Q);

if sum(abs(P - Q)) ~= 0
    S.kulczynski_s = sum(min(P, Q)) / sum(abs(P - Q));
else
    S.kulczynski_s = sum(min(P1, Q1)) / sum(abs(P1 - Q1));
    Error.kulczynski_s = 'Divide by Zero Error';
end

if sum(min(P, Q)) ~= 0
    D.kulczynski_s = 1 / S.kulczynski_s;
else
    D.kulczynski_s = 1 / (sum(min(P1, Q1)) / sum(abs(P1 - Q1)));
    Error.kulczynski_s = 'Divide by Zero Error';
end

S.ruzicka = sum(min(P, Q)) / sum(max(P, Q));
D.ruzicka = 1 - S.ruzicka;

D.tanimoto = (sum(P) + sum(Q) - 2 * sum(min(P, Q))) / (sum(P) + sum(Q) - sum(min(P, Q));

% Inner product family
S.inner = dot(P, Q);
D.inner = 1 - S.inner;

if min(P + Q) ~= 0
    S.harmonic = 2 * sum((P .* Q) ./ (P + Q));
    D.harmonic = 1 - S.harmonic;
else
    S.harmonic = 2 * sum((P1 .* Q1) ./ (P1 + Q1));
    D.harmonic = 1 - S.harmonic;
    Error.harmonic = 'Divide by Zero Error';
end

S.cosine = dot(P, Q) / sqrt(dot(P, P) * dot(Q, Q));
D.cosine = 1 - S.cosine;

S.kumarh = dot(P, Q) / (dot(P, P) + dot(Q, Q) - dot(P, Q));
D.kumarh = 1 - S.kumarh;

S.jaccard = dot(P, Q) / (dot(P, P) + dot(Q, Q) - dot(P, Q));
D.jaccard = dot(P - Q, P - Q) / (dot(P, P) + dot(Q, Q) - dot(P, Q));

S.dice = 2 * dot(P, Q) / (dot(P, P) + dot(Q, Q));
D.dice = 1 - 2 * dot(P, Q) / (dot(P, P) + dot(Q, Q));

% Fidelity (square chord) family
S.fidelity = sum(sqrt(P .* Q));
D.fidelity = 1 - S.fidelity;

if sum(sqrt(P .* Q)) ~= 0
    D.bhattacharrya = -log(sum(sqrt(P .* Q));
else
    D.bhattacharrya = -log(sum(sqrt(P1 .* Q1));
    Error.bhattacharrya = 'log(0) Error';
end

D.hellinger = sqrt(2 * sum((sqrt(P) - sqrt(Q)).^2));
D.squaredchord = dot(sqrt(P) - sqrt(Q), sqrt(P) - sqrt(Q));
S.squaredchord = 2 * sum(sqrt(P .* Q)) - 1;
D.matusita = sqrt(dot(sqrt(P) - sqrt(Q), sqrt(P) - sqrt(Q));

% Squared L2 (chi-squared) family
D.squaredeuclidean = dot(P - Q, P - Q);

if min(Q) ~= 0
    D.pearsonchi = sum((P - Q).^2 ./ Q);
else
    D.pearsonchi = sum((P1 - Q1).^2 ./ Q1);
    Error.pearsonchi = 'Divide by Zero Error';
end

if min(P) ~= 0
    D.neymanchi = sum((P - Q).^2 ./ P);
else
    D.neymanchi = sum((P1 - Q1).^2 ./ P1);
    Error.neymanchi = 'Divide by Zero Error';
end

if min(P + Q) ~= 0
    D.squaredchi = sum((P - Q).^2 ./ (P + Q));
else
    D.squaredchi = sum((P1 - Q1).^2 ./ (P1 + Q1));
    Error.squaredchi = 'Divide by Zero Error';
end

D.probsymm = 2 * D.squaredchi;

if min(P + Q) ~= 0
    D.divergence = 2 * sum((P - Q).^2 ./ (P + Q).^2);
    D.clark = sqrt(sum((abs(P - Q) ./ (P + Q)).^2);
else
    D.divergence = 2 * sum((P1 - Q1).^2 ./ (P1 + Q1).^2);
    D.clark = sqrt(sum((abs(P1 - Q1) ./ (P1 + Q1)).^2);
    Error.divergence = 'Divide by Zero Error';
    Error.clark = 'Divide by Zero Error';
end

if min(P .* Q) ~= 0
    D.additivesymm = sum((P - Q).^2 .* (P + Q) ./ (P .* Q));
else
    D.additivesymm = sum((P1 - Q1).^2 .* (P1 + Q1) ./ (P1 .* Q1));
    Error.additivesymm = 'Divide by Zero Error';
end

% Shannon's entropy family
if min(Q) ~= 0 && min(P ./ Q) ~= 0
    D.kullback_PQ = sum(P .* log(P ./ Q));
else
    D.kullback_PQ = sum(P1 .* log(P1 ./ Q1));
    Error.kullback_PQ = 'Divide by Zero or log(0) Error';
end

if min(P) ~= 0 && min(Q ./ P) ~= 0
    D.kullback_QP = sum(Q .* log(Q ./ P));
else
    D.kullback_QP = sum(Q1 .* log(Q1 ./ P1));
    Error.kullback_QP = 'Divide by Zero or log(0) Error';
end

if min(Q) ~= 0 && min(P ./ Q) ~= 0
    D.jeffreys = sum((P - Q) .* log(P ./ Q));
else
    D.jeffreys = sum((P1 - Q1) .* log(P1 ./ Q1));
    Error.jeffreys = 'Divide by Zero or log(0) Error';
end

if min(P + Q) ~= 0 && min(P ./ (P + Q)) ~= 0
    D.kdivergence = sum(P .* log(2 * P ./ (P + Q)));
else
    D.kdivergence = sum(P1 .* log(2 * P1 ./ (P1 + Q1));
    Error.kdivergence = 'Divide by Zero or log(0) Error';
end

if min(P + Q) ~= 0 && min(P ./ (P + Q)) ~= 0 && min(Q ./ (P + Q)) ~= 0
    D.topsoe = sum(P .* log(2 * P ./ (P + Q)) + Q .* log(2 * Q ./ (P + Q)));
else
    D.topsoe = sum(P1 .* log(2 * P1 ./ (P1 + Q1)) + Q1 .* log(2 * Q1 ./ (P1 + Q1)));
    Error.topsoe = 'Divide by Zero or log(0) Error';
end

if min(Q) ~= 0 && min(P ./ Q) ~= 0 && min(P) ~= 0 && min(Q ./ P) ~= 0
    D.jensen_s = 0.5 * (D.kullback_PQ + D.kullback_QP);
else
    D.jensen_s = 0.5 * (D.kullback_PQ + D.kullback_QP);
    Error.jensen_s = 'Divide by Zero or log(0) Error';
end

if min(P + Q) ~= 0 && min(P) ~= 0 && min(Q) ~= 0
    D.jensen_d = sum(0.5 * (P .* log(P) + Q .* log(Q) - (P + Q) .* log(0.5 * (P + Q))));
else
    D.jensen_d = sum(0.5 * (P1 .* log(P1) + Q1 .* log(Q1) - (P1 + Q1) .* log(0.5 * (P1 + Q1))));
    Error.jensen_d = 'log(0) Error';
end

% Combinations
if min(dot(P, Q)) ~= 0 && min((P + Q) ./ sqrt(dot(P, Q)) ~= 0
    D.taneja = sum(0.5 * (P + Q) .* log(0.5 * (P + Q) ./ sqrt(dot(P, Q))));
else
    D.taneja = sum(0.5 * (P1 + Q1) .* log(0.5 * (P1 + Q1) ./ sqrt(dot(P1, Q1)));
    Error.taneja = 'Divide by Zero or log(0) Error';
end

if min(P .* Q) ~= 0
    D.kumarj = 0.5 * sum((P.^2 - Q.^2).^2 ./ (P .* Q).^(3/2));
else
    D.kumarj = 0.5 * sum((P1.^2 - Q1.^2).^2 ./ (P1 .* Q1).^(3/2));
end

D.avgL = 0.5 * (sum(abs(P - Q)) + max(abs(P - Q)));

% Vicissitude
if min(min(P, Q)) ~= 0
    D.viciswave = sum(abs(P - Q) ./ min(P, Q));
    D.vicissymm1 = sum((P - Q).^2 ./ min(P, Q).^2);
    D.vicissymm2 = sum((P - Q).^2 ./ min(P, Q));
else
    D.viciswave = sum(abs(P1 - Q1) ./ min(P1, Q1));
    D.vicissymm1 = sum((P1 - Q1).^2 ./ min(P1, Q1).^2);
    D.vicissymm2 = sum((P1 - Q1).^2 ./ min(P1, Q1));
    Error.viciswave = 'Divide by Zero Error';
    Error.vicissymm1 = 'Divide by Zero Error';
    Error.vicissymm2 = 'Divide by Zero Error';
end

if min(max(P, Q)) ~= 0
    D.vicissymm3 = sum((P - Q).^2 ./ max(P, Q));
else
    D.vicissymm3 = sum((P1 - Q1).^2 ./ max(P1, Q1));
    Error.vicissymm3 = 'Divide by Zero Error';
end

if min(Q) ~= 0 && min(P) ~= 0
    D.maxsymm = max(D.pearsonchi, D.neymanchi);
    D.minsymm = min(D.pearsonchi, D.neymanchi);
else
    D.maxsymm = max(D.pearsonchi, D.neymanchi);
    D.minsymm = min(D.pearsonchi, D.neymanchi);
    Error.maxsymm = 'Divide by Zero Error';
    Error.minsymm = 'Divide by Zero Error';
end

