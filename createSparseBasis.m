function [B,hashmap,invmap,X,Y] = createSparseBasis(N)
%Create sparse DG basis \hat{V}_N^0 (constant functions)

%Intput parameters -----
%   N::int -- level.  Full grid would by 2^N basis functions in each
%       dimension

%Output variables ------
%   B::cell -- Cell containing the sparse basis vectors.  Written as a
%       matrix structrue as if on the domain [0,1]\times[0,1]
%   hashmap::hashtable -- map from [l1,j1,l2,j2] -> basis function
%       for j1=0:2^(l1-1)-1, j2=0:2^(l2-1)-1
%   invmap::cell -- inverse map of hashmap.  
%   [X,Y] - meshgrid of the finiest mesh on [0,1]\times[0,1]

%Create [X,Y] with meshgrid
[X,Y] = meshgrid((1/2^N:1/2^N:1)-1/2^(N+1),(1/2^N:1/2^N:1)-1/2^(N+1));


%Create V_N^0 in 1D.  W is the container.
W = zeros(2^N);
W(1,:) = ones(1,2^N);
count = 2;
for n=1:N
    chunk = 2^(N-n);
    %2^(n-1) basis vectors
    for i=1:2^(n-1)
        temp = zeros(1,2^N);
        temp(1+2*(i-1)*chunk:2*(i-1)*chunk+chunk) = -1*ones(1,chunk);
        temp(2*(i-1)*chunk+chunk+1:2*i*chunk) = 1*ones(1,chunk);
        temp = temp*2^((n-1)/2);
        W(count,:) = temp;
        count = count + 1;
    end
end    


l2gst = @(n) 1*(n == 0) + (2^(n-1)+1).*(n > 0); %local-2-global start
l2ged = @(n) 1*(n == 0) + (2^(n)    ).*(n > 0); %local-2-global end


%Construct sparse grid space
B = cell(2^N,1);
keys = cell(2^N,1);
values = cell(2^N,1);
invmap = cell(2^N,1);
count = 1;
for l1=0:N
    for j1=l2gst(l1):l2ged(l1)
        %Since d=2, multi-indices l=[l1,l2] with |l|_1 \leq N
        %require l2\leq N-1
        for l2=0:N-l1
            for j2=l2gst(l2):l2ged(l2)
                keys{count} = sprintf('%02di%02di%02di%02d',l1,j1,l2,j2);
                values{count} = count;
                invmap{count} = [l1,j1-1,l2,j2-1];
                B{count} = flipud(W(j1,:)')*W(j2,:); %%flipud so functions
                    %match mesh; that is, entry (1,1) corresponds to
                    %x=0,y=1.  
                count = count + 1;
            end
        end
    end
end

%Collapse unneeded values
B = B(1:count-1);
keys = keys(1:count-1);
values = values(1:count-1);
invmap = invmap(1:count-1);
hashmap = containers.Map(keys,values);

end

