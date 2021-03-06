function [gas]= gwr_core(eta,gas)

A = gas.A;
C = gas.C;
C_age = gas.C_age;
h = gas.h;
r = gas.r;
hizero = gas.hizero;
hszero = gas.hszero;
awk = gas.awk;
params = gas.params;
en = params.en;
eb = params.eb;

%%%%%%%%%%%%%%%%%%% ATTENTION STILL MISSING FIRING RATE! will have problems
%%%%%%%%%%%%%%%%%%% when algorithm not static!!!!
%%%%%%%%%%%%%%%%%%%

%eta = data(:,k); % this the k-th data sample
[ws, ~, s, t, ~] = findnearest(eta, A); %step 2 and 3
if C(s,t)==0 %step 4
    C = spdi_bind(C,s,t);
else
    C_age = spdi_del(C_age,s,t);
end
a = exp(-norm((eta-ws).*awk)); %step 5

%algorithm has some issues, so here I will calculate the neighbours of
%s
[neighbours] = findneighbours(s, C);
num_of_neighbours = size(neighbours,2);

if a < params.at && r <= params.nodes %step 6
    wr = 0.5*(ws+eta); %too low activity, needs to create new node r
    A(:,r) = wr;
    C = spdi_bind(C,t,r);
    C = spdi_bind(C,s,r);
    C = spdi_del(C,s,t);
    r = r+1;
else %step 7
    for j = 1:num_of_neighbours % check this for possible indexing errors
        i = neighbours(j);
        %size(A)
        wi = A(:,i);
        A(:,i) = wi + en*h(i)*(eta-wi);
    end
    A(:,s) = ws + eb*h(s)*(eta-ws); %adjusts nearest point MORE;;; also, I need to adjust this after the for loop or the for loop would reset this!!!
end
%step 8 : age edges with end at s
%first we need to find if the edges connect to s

for j = 1:num_of_neighbours % check this for possible indexing errors
    i = neighbours(j);
    C_age = spdi_add(C_age,s,i);
end

%step 9: again we do it inverted, for loop first
%%%% this strange check is a speedup for the case when the algorithm is static
if params.STATIC % skips this if algorithm is static
    h = hizero;
    h(s) = hszero;
else
    for i = 1:r %%% since this value is the same for all I can compute it once and then make all the array have the same value...
        h(i) = gas.hi(gas.time,params); %ok, is this sloppy or what? t for the second nearest point and t for time
    end
    h(s) = gas.hs(gas.time,params);
    gas.time = (cputime - gas.t0)*1;
end

%step 10: check if a node has no edges and delete them
%[C, A, C_age, h, r ] = removenode(C, A, C_age, h, r);
%check for old edges

% makes the algorithm slightly faster when the matrices are not full
if r > params.nodes
    R = params.nodes;
else
    R = r;
end

if r>2 % don't remove everything
    
    [C(1:R,1:R), C_age(1:R,1:R) ] = removeedge(C(1:R,1:R), C_age(1:R,1:R), params.amax);
    [C(1:R,1:R), A(:,1:R), C_age(1:R,1:R), h, r ] = removenode(C(1:R,1:R), A(:,1:R), C_age(1:R,1:R), h, r);  %inverted order as it says on the algorithm to remove points faster
end

gas.A = A;
gas.C = C;
gas.C_age = C_age;
gas.h = h;
gas.r = r;
gas.a = a;
end

function sparsemat = spdi_add(sparsemat, a, b) %increases the number so that I don't have to type this all the time and forget it...
sparsemat(a,b) = sparsemat(a,b) + 1;
sparsemat(b,a) = sparsemat(a,b) + 1;
end

function sparsemat = spdi_bind(sparsemat, a, b) % adds a 2 way connection, so that I don't have to type this all the time and forget it...
sparsemat(a,b) = 1;
sparsemat(b,a) = 1;
end

function sparsemat = spdi_del(sparsemat, a, b) % removes a 2 way connection, so that I don't have to type this all the time and forget it...
sparsemat(a,b) = 0;
sparsemat(b,a) = 0;
end

function [C, C_age ] = removeedge(C, C_age, amax) 
[row, col] = find(C_age > amax);
a = size(row,2);
if ~isempty(row)
    for i = 1:a
        C_age(row(i),col(i)) = 0;
        C_age(col(i),row(i)) = 0;
        C(row(i),col(i)) = 0;
        C(col(i),row(i)) = 0;
    end
end
end

function [C, A, C_age, h,r ] = removenode(C, A, C_age, h,r) %depends only on C operates on everything

[row,~] = find(C);
%a = [row;col];
maxa = max(row);
for i = 1:maxa% r
    if max(row)<maxa %ok, lets try this, if the old maximum is not valid anymore stop the for loop.
        break % I am assuming that this also means that all of the remaining rows and columns are zeroed
    end
    if isempty(find(row == i, 1))
        
        %has to do this to every matrix and vector
        C = clipsimmat(C,i);
        if i>size(A,2)
            disp('wrong stuff going on')
        end
        A = clipA(A,i); 
        C_age = clipsimmat(C_age,i);
        h = clipvect(h,i);
        r = r-1;
        if r<1||r~=fix(r)
            error('something fishy happening. r is either zero or fractionary!')
        end
        [row,~] = find(C);
    end
   
end
end

function C = clipsimmat(C,i)
C(i,:) = [];
C(:,i) = [];
C = padarray(C,[1 1],'post');
end

function V = clipvect(V, i)
V(i) = [];
V = [V 0];
end

function A = clipA(A, i)
A(:,i) = [];
ZERO = zeros(size(A,1),1);
A = [A ZERO];
end

function neighbours = findneighbours(s,C)
[row, col] = find(C); % there will likely be infinite mistakes in indexing here...
neighbours = [];
for i = 1:length(row)
    if row(i) == s
        %s
        %i
        neighbours = [neighbours col(i)];
        %C_age = spdi_add(C_age,s,col(i)); %omg, the horror.... but s = row(i) already is this it? I have no idea...
        %cummax(cummax(C_age))
    end
end
%neighbours
end

