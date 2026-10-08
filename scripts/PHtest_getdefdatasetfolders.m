function fldr0 = PHtest_getdefdatasetfolders
% fldr0 = PHtest_getdefdatasetfolders
%
% Returns cell vector containing dataset folder names.

% defaults
sets = [100:102,200:202,300:302,400:402,500,600,210:212,310:312,410:412,...
    320:322]; % indexes of datasets

% build cell array containing dataset fodler names
fldr0 = cellstr([repmat('dataset',size(sets,2),1),num2str(sets')])';
fldr0 = [fldr0, 'EBS-IBS', 'D135'];