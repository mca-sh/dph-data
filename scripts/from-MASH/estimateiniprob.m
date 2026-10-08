function ipv = estimateiniprob(tp,ip,id)
% ip = estimateiniprob(tp,ip,id)
%
% Estimate initial probabilities of state aggregate in absorbing Markov 
% model from full Markov model transition probability matrix.
%
% tp: [J-by-J] full Markov model transition probability matrix
% ip: [1-by-J] initial state probabilities in full Markov model.
% id: [1-by-D] indexes in tp that correspond to aggregated states or,
%     [1-by-J] logicals indicating if state is in aggregate
% ipv: [1-by-D] initial state probabilities of aggregated states in
%     absorbing Markov model.

% initialize output
ipv = [];

% transform state indexes to logicals
if ~islogical(id)
    logid = false(1,size(tp,1));
    logid(id) = true;
    id = logid;
end

% % determine steady-state probabilities
% if sum(ip(~id))>0
%     wght = ip(~id)'/sum(ip(~id));
% else
%     wght = ones(sum(~id),1)/sum(~id);
% end
% tpv = [tp(id,id),sum(tp(id,~id),2);sum(wght.*tp(~id,id),1),0];
% tpv = tpv./repmat(sum(tpv,2),1,(sum(id)+1));
% tpv(isnan(tpv)) = 0;
% [~,Dev,Wev] = eig(tpv);
% ev = rnd2tol(diag(Dev),eps("double"));
% isev1 = find(ev==1);
% for j = 1:numel(isev1)
%     if all(Wev(:,isev1(j))>=0)
%         ipv = Wev(id,isev1(j))';
%         break
%     end
% end

% estimate initial state probabilities from summed transition probabilities
if sum(ipv)==0
    W = tp;
    W(~~eye(size(tp,1))) = 0;
    W = W./repmat(sum(W,2),1,size(W,2));
    ipv = ip(id)+sum(repmat(ip(~id)',[1,sum(id)]).*W(~id,id),1);
end

% normalize initial state probabilities
ipv = ipv/sum(ipv);
