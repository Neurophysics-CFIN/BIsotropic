function Q = canon_antipodes(Q)
% Pick canonical representative of q ≡ -q (lexicographic tie-break)
if size(Q,2) == 4
    neg = Q(:,1) < 0 | (Q(:,1)==0 & (Q(:,2) < 0 | (Q(:,2)==0 & (Q(:,3) < 0 | (Q(:,3)==0 & Q(:,4) < 0)))));
    Q(neg,:) = -Q(neg,:);
elseif size(Q,2) ==3
    neg = Q(:,3) < 0 | (Q(:,3)==0 & (Q(:,2) < 0 | (Q(:,2)==0 & (Q(:,1) < 0 ))));
    Q(neg,:) = -Q(neg,:);
else
    error('Input matrix Q must have 3 or 4 columns.');
end
