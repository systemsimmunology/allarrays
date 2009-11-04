
% Load values 
load TFAbooleVals.tsv

% Load gene ids and symbols
fp = fopen('TFAbooleGeneIDs.tsv')
C = textscan( fp, '%s\t%s')
gsymbols = C{:,2}


% Hamming distance object
distobj = pdist(TFAbooleVals, 'hamming' )

% Stress curve over possible choices of dimensionality

d=10;
S=zeros(1,d);

for i = 2:10
    [yy,stress] = mdscale(distobj,i);
    S(i)=stress;
end



% Disstances for 3 dimensions
yy=mdscale(distobj,3)


% 3D plot, with labels

plot3(yy(:,1),yy(:,2),yy(:,3),'.')
text(yy(:,1),yy(:,2),yy(:,3),gsymbols )

