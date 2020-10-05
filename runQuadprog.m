%Projecting a sparse DG function onto the set of positive sparse
%DG functions while maintaing conservation

% minimize F(a) := (1/2)\| a - \bar{a} \|_2^2
% subject to:
%                    -Bmat*a  \leq  0
%                   e_1^T * a   =   0

%Here \bar{a} are the coordinates for sthandard L^2 projeciton which is not
%guaranteed to be positive even if the projecting function is positive

%%% Data ------------------

%Finest level
N = 4;

%function to project

        %f(x,y) = x*y
f = @(x,y) x.*y;
        % Indicator function on [1/3,2/3]^2
%f = @(x,y) (x >= 1/3).*(x <= 2/3).*(y >= 1/3).*(y <= 2/3);

%%%------------------------

%Create sparse DG basis
[B,hashmap,invmap,X,Y] = createSparseBasis(N);

%Create Bmat, the matrix that maps coordinates to piecewise constant
%functions
M = numel(B);
Bmat = zeros(numel(B{1}),M);
for i=1:M
    Bmat(:,i) = B{i}(:)/sqrt(numel(B{1})); %Scaling to make inner product
                            % on \W=[0,1]^2 equivalent to the dot product
end

%Get L2 projection of f(x,y) = xy onto sparse DG space.
proj = calcProj(B,f);

%L^2 projection coordinates
abar = Bmat'*proj(:);

%Intial point (not feasible but whatever).  Not used by interior-point
%algorithm
x0 = abar;

%F(a) can be written as F(a) = a'*I*a + f_quad'*a where f_quad = -\bar(a)
%Note that minimizing against \bar{a} is the same as minimizing against f
f_quad = -abar;

%Inequality constaints Ax \leq b can be enforeced by A=-Bmat and b being
%the zero vector
b = zeros(numel(B{1}),1);

%Equality constants for global conservation.  Since the integral of every
%basis function besides the first one vanishes on \W, only need to enforce
%the first one
Aeq = zeros(1,M);
Aeq(1) = 1;
beq = abar(1);

%Set quadprog options
options = optimoptions('quadprog','Display','iter-detailed',...
    'OptimalityTolerance',1e-12);

%Run quadprog
[x,fval,exitflag,output] = quadprog(speye(M),f_quad,-Bmat,b,Aeq,beq,[],[],x0,options);

%Produce solution from coordinates
uQuad = reshape(Bmat*x,sqrt(numel(B{1})),sqrt(numel(B{1})));

%Verify positivity
fprintf("Number of elements that are negative is\n");
fprintf("\t%d for Standard L2 projection\n",sum(proj(:) < -(1e-12)));
fprintf("\t%d for Positive L2 projection\n",sum(uQuad(:) < -(1e-12)));
fprintf('\n');

%Verify global conservation |\int_W v-v_pos \dx|
fprintf("Global conservation error is %e\n",abs(sum(proj(:)-uQuad(:)))/numel(B{1}));
fprintf('\n');

%Calculate errors
fprintf('L2 error is\n');
fprintf('\t%e for Standard L2 projection\n',errFuncSparse(proj,f));
fprintf('\t%e for Positive L2 projection\n',errFuncSparse(uQuad,f));

