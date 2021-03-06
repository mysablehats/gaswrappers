function [n1, n2, ni1, ni2, distvector] = findnearest(varargin)
%% findnearest
% use as: findnearest(p, data)
% or findnearest() without arguments to run test
%
% the original function had a sort, which is awfully slow!!!!
%so I am rewriting to see if I can make it better
%
% I think pdist2 'smallest' is much faster than this, but I haven't
% benchmarked it...

if isempty(varargin)||strcmp(varargin{1},'test')
    findnearest_test()
else
    p = varargin{1};
    data =varargin{2};
    
%     if length(varargin)>=3&&strcmp(varargin{3},'slow') %% this is the fastest
%         
%         [n1, n2, ni1, ni2, distvector] = findnearest_slow(p,data);
%         
%     elseif length(varargin)>=3&&strcmp(varargin{3},'notsofast')
%         
%         [n1, n2, ni1, ni2, distvector] = findnearest_faster(p,data);
%     elseif length(varargin)>=3&&strcmp(varargin{3},'bad')
%         
%         [n1, n2, ni1, ni2, distvector] = findnearest_fastest(p,data);
%     else
         [n1, n2, ni1, ni2, distvector] = findnearest_fast(p,data);
%    end
end
end
function findnearest_test()
%do some test
p = rand(1000,10000);
data = rand(1000,30000);
disp('how much time it takes to run it slow')
tic()
for i = 1:10000
    [n1f, n2f, ni1f, ni2f, distvectorf] = findnearest(p(:,i),data,'fast');
end
toc()
% disp('how much time it takes to run it fast')
% tic()
% for i = 1:10000
%     [n1, n2, ni1, ni2, distvector] = findnearest(p(:,i),data,'fast');
% end
% toc()
% disp('how much time it takes to run it fastest')
% tic()
% for i = 1:10000
%     [n1, n2, ni1, ni2, distvector] = findnearest(p(:,i),data,'notsofast');
% end
% toc()
% 
% 
% if all(n1f==n1)&&all(n2f==n2)&&all(ni1f==ni1)&&all(ni2f==ni2)&&all(distvectorf==distvector)
%     disp('ERFOLG: ALL Same results!')
% else
%     disp('Disagreeing results, check functions!')
% end
% 
% if all(n1f==n1)
%     disp('ERFOLG: N1 IS THE SAME !')
% else
%     n1f
%     n1
%     disp('Disagreeing results, check functions! N1 IS NOT THE SAME ')
%     
% end
% 
% if all(n2f==n2)
%     disp('ERFOLG:N2 IS THE SAME ')
% else
%     n2f
%     n2
%     disp('Disagreeing results, check functions! N2 IS NOT THE SAME')
%     
% end
% 
% if all(ni1f==ni1)
%     disp('ERFOLG: N1 index IS THE SAME!')
% else
%     ni1f
%     ni1
%     disp('Disagreeing results, check functions! N1 index IS NOT THE SAME ')
% end
% 
% if all(ni2f==ni2)
%     disp('ERFOLG: N2 index IS THE SAME!!')
% else
%     ni2f
%     ni2
%     disp('Disagreeing results, check functions! N2 index IS NOT THE SAME')
% end
% 
% if all(distvectorf==distvector)
%     disp('ERFOLG: DISTVECTOR, - no one uses this... !')
% else
%     disp('Disagreeing results, check functions! DISTVECTOR, - no one uses this...')
% end



end
function [n1, n2, ni1, ni2, distvector] = findnearest_slow(p,data)
maxindex = size(data,2);
distvector = zeros(1,maxindex);
for i = 1:maxindex
    distvector(i) = norm(data(:,i)- p);
end
[~, index] = sort(distvector);
ni1 = index(1);
ni2 = index(2);
n1 = data(:,ni1);
n2 = data(:,ni2);
end
function [n1, n2, ni1, ni2, distvector] = findnearest_fast(p,data)
maxindex = size(data,2);
distvector = zeros(1,maxindex);
for i = 1:maxindex
    distvector(i) = norm(data(:,i)- p);
end
[~,ni1] = min(distvector);
n1 = data(:,ni1);
pushdist = distvector(ni1);
distvector(ni1) = NaN; % I use some cleverness. Hopefully this is fast.
[~,ni2] = min(distvector);
n2 = data(:,ni2);
distvector(ni1) = pushdist;
end
function [n1, n2, ni1, ni2, distvector] = findnearest_faster(p,data)
maxindex = size(data,2);
distvector = inf(1,maxindex);
ni1 = 0;
ni2 = 0;
a = 0;
for i = 1:maxindex
    a = norm(data(:,i)- p);
    if a < distvector(1)
        distvector(2) = distvector(1);
        ni2 = ni1;
        distvector(1) = a;
        ni1 = i;
    elseif a < distvector(2)
        distvector(2) = a;
        ni2 = i;
    end
end
n1 = data(:,ni1);
n2 = data(:,ni2);
end
function [n1, n2, ni1, ni2, distvector] = findnearest_fastest(p,data)
maxindex = size(data,2);
distvector = inf(1,maxindex);
ni1 = 0;
ni2 = 0;
a = 0;
ppp = repmat(p,1,maxindex);
dada = data - ppp;
for i = 1:maxindex
    a = norm(dada(:,i));
    if a < distvector(1)
        distvector(2) = distvector(1);
        ni2 = ni1;
        distvector(1) = a;
        ni1 = i;
    elseif a < distvector(2)
        distvector(2) = a;
        ni2 = i;
    end
end
n1 = data(:,ni1);
n2 = data(:,ni2);
end