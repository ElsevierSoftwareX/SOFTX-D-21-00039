%% Execute tests in the console by running: 
% result = runtests('unitTests')

%% Main function
function tests = unitTests
    tests = functiontests(localfunctions);
end


%% Unit tests

% cdvinearray
function test_cdvinearray(testCase)

% assert errors
assertError(testCase,@()cdvinearray('t',1),'cdvinearray:WrongVineType')

% invalid dimensions
actual = cdvinearray('c',-3);
expected = [];

assertEqual(testCase,actual,expected)

% test 1: 1-dimensional C-vine
actual = cdvinearray('c',1);
expected = 1;

assertEqual(testCase,actual,expected)

% test 2: 1-dimensional D-vine
actual = cdvinearray('d',1);

assertEqual(testCase,actual,expected)

% test 3: 3-dimensional C-vine
actual = cdvinearray('c',3);
expected = [1 1 1; 0 2 2; 0 0 3];

assertEqual(testCase,actual,expected)

% test 4: 3-dimensional D-vine
actual = cdvinearray('d',3);
expected = [1 1 2; 0 2 1; 0 0 3];

assertEqual(testCase,actual,expected)    
    
end


% check_cube
function test_check_cube(testCase)

% test 1.1: line is ok
line = {0.2 0.7 1; 2 3 4};
actual = check_cube(line);
expected = 1;

assertEqual(testCase,actual,expected)

% test 1.2: line with 0 at the beginning 
line = {0 0.2 0.7 1; 1 2 3 4};
actual = check_cube(line);
expected = 0;

assertEqual(testCase,actual,expected)

% test 1.3: line without 1 at the end 
line = {0.2 0.7; 2 3};
actual = check_cube(line);
expected = 0;

assertEqual(testCase,actual,expected)

% test 1.4: line not ordered
line = {0.7 0.2 1; 2 3 4};
actual = check_cube(line);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2.1: square is ok
square = {[0 0.5;0 1],[0.5 1; 0 1]; -0.5, 0.5};
actual = check_cube(square);
expected = 1;

assertEqual(testCase,actual,expected)

% test 2.2: square not ordered properly, case 1
square = {[0.5 0;1 0],[0.5 1; 0 1]; -0.5, 0.5};
actual = check_cube(square);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2.3: square not ordered properly, case 2
square = {[0 0;0.5 1],[0.5 1; 0 1]; -0.5, 0.5};
actual = check_cube(square);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2.4: square doesn't contain (0,0)
square = {[0.1 0.6;0 1],[0.5 1; 0 1]; -0.5, 0.5};
actual = check_cube(square);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2.5: square doesn't contain (1,1)
square = {[0 0.5;0 1],[0.6 1; 0 1.1]; -0.5, 0.5};
actual = check_cube(square);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2.6: area less than 1
square = {[0 0.5;0 1],[0.5 1; 0 0.8]; -0.5, 0.5};
actual = check_cube(square);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2.7: area larger than 1
square = {[0 0.5;0 1],[0.4 1; 0 1]; -0.5, 0.5};
actual = check_cube(square);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2.8: square incomplete, area is 1, though
square = {[0 0.5;0 1],[0 0.5;0 1]; -0.5, 0.5};
actual = check_cube(square);
expected = 0;

assertEqual(testCase,actual,expected)
    
end


% copulapdfadv
function test_copulapdfadv(testCase)

u = [0.5 0.5];

% assert errors
assertError(testCase,@()copulapdfadv('gumbel',u,0),'copulapdfadv:InvalidParameter')

% test 1: independence
actual = copulapdfadv('ind',u,0);
expected = 1;

assertEqual(testCase,actual,expected)

% test 2: all families
families = {'gumbel','clayton','frank','t','gauss','ind',...
            'amhaq','tawn','fgm','plackett','joe','surclayton',...
            'surgumbel','surjoe'};
theta = {2,2,2,[0.3,3],0.5,0,...
         0.5,0.5,0.5,2,2,2,...
         2,2};

for ii = 1:1:length(families)
    
    actual = copulapdfadv(families{ii},u,theta{ii});
    
    assertSize(testCase,actual,[1 1])
    
end % ii

end


% copulaselect
function test_copulaselect(testCase)

n = 100;
u = copulasim('clayton',5,n);

% assert errors
assertError(testCase,@()copulaselect(u),'copulaselect:InvalidNumberOfInputs')

% test 1: copula family pre-specified
[actualfam,actualthetahat,actualloglik,actualmcrit] = copulaselect(u,'clayton','frank');

assertTrue(testCase,any(strcmp(actualfam,{'clayton','frank'})))
assertSize(testCase,actualthetahat,[1 1])
assertSize(testCase,actualloglik,[n 1])
assertSize(testCase,actualmcrit,[1 1])

% test 2: copula family not pre-specified
families = {'gumbel','clayton','frank','t','gauss','ind',...
            'amhaq','tawn','fgm','plackett','joe','surclayton',...
            'surgumbel','surjoe'};
[actualfam,actualthetahat,actualloglik,actualmcrit] = copulaselect(u,'sll');

assertTrue(testCase,any(strcmp(actualfam,families)))
assertSize(testCase,actualthetahat,[1 1])
assertSize(testCase,actualloglik,[n 1])
assertSize(testCase,actualmcrit,[1 1])
assertEqual(testCase,actualmcrit,sum(actualloglik))

end


% copulasim
function test_copulasim(testCase)

families = {'gumbel','clayton','frank','t','gauss','ind',...
            'amhaq','tawn','fgm','plackett','joe','surclayton',...
            'surgumbel','surjoe'};
theta = {2,2,2,[0.3,3],0.5,0,...
         0.5,0.5,0.5,2,2,2,...
         2,2};
n = 2;

% assert errors
assertError(testCase,@()copulasim('gumbel',0,n),'copulasim:InvalidParameter')

% test 1: simulate n points from each supported copula
for ii = 1:1:length(families)
    
    actual = copulasim(families{ii},theta{ii},n);
    
    assertSize(testCase,actual,[n 2])
    
end % ii

end


% cpcheck
function test_cpcheck(testCase)

families = {'gumbel','clayton','t','gauss',...
            'amhaq','tawn','fgm','plackett','joe','surclayton',...
            'surgumbel','surjoe'};

% test 1: everything ok
theta = {2,2,[0.3,3],0.5,...
         0.5,0.5,0.5,2,2,2,...
         2,2};
expected = 1;
     
for ii = 1:1:length(families)
    
    actual = cpcheck(families{ii},theta{ii});
    
    assertEqual(testCase,actual,expected)
    
end % ii    

% test 2: nothing ok
theta = {-1,-1,[5,-1],4,...
         -2,-1,-2,-1,-1,-1,...
         -1,-1};
expected = 0;

for ii = 1:1:length(families)
    
    actual = cpcheck(families{ii},theta{ii});
    
    assertEqual(testCase,actual,expected)
    
end % ii    

end


% ecopula
function test_ecopula(testCase)

% test 1: only one data point
u = [0.5 0.5];
actual = ecopula(u);
expected = 0.5;

assertEqual(testCase,actual,expected)

% test 2: only one data point, based on one data point
v = [1 1];
actual = ecopula(u,v);
expected = 0;

assertEqual(testCase,actual,expected)

end


% goftest_a2
function test_goftest_a2(testCase)

A = [1 1 2; 0 2 1; 0 0 3];
family = {'gumbel', 'gumbel';'gumbel',0};
theta = {2,3;2,0};
n = 10;
u = simrvine(n,A,family,theta,0);

% assert errors
theta_err = {0,3;2,0};
assertError(testCase,@()goftest_a2(u,A,family,theta_err,0),'goftest_a2:InvalidParameter')

A_err = [1 1 2; 0 2 1];
assertError(testCase,@()goftest_a2(u,A_err,family,theta,0),'goftest_a2:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()goftest_a2(u,A_err,family,theta,0),'goftest_a2:InvalidVineArray')

assertError(testCase,@()goftest_a2(u,A,family,theta,2),'goftest_a2:InvalidParallelizationMode')

end


% goftest_a4
function test_goftest_a4(testCase)

A = [1 1 2; 0 2 1; 0 0 3];
family = {'gumbel', 'gumbel';'gumbel',0};
theta = {2,3;2,0};
n = 10;
u = simrvine(n,A,family,theta,0);

% assert errors
theta_err = {0,3;2,0};
assertError(testCase,@()goftest_a4(u,A,family,theta_err,0),'goftest_a4:InvalidParameter')

A_err = [1 1 2; 0 2 1];
assertError(testCase,@()goftest_a4(u,A_err,family,theta,0),'goftest_a4:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()goftest_a4(u,A_err,family,theta,0),'goftest_a4:InvalidVineArray')

assertError(testCase,@()goftest_a4(u,A,family,theta,2),'goftest_a4:InvalidParallelizationMode')

end


% hfunc
function test_hfunc(testCase)

u = 0.5;
v = 0.5;

% assert errors
assertError(testCase,@()hfunc(u,v,'gumbel',0),'hfunc:InvalidParameter')

% test 1: test all families
families = {'gumbel','clayton','frank','t','gauss','ind',...
            'amhaq','tawn','fgm','plackett','joe','surclayton',...
            'surgumbel','surjoe'};
theta = {2,2,2,[0.3,3],0.5,0,...
         0.5,0.5,0.5,2,2,2,...
         2,2};

for ii = 1:1:length(families)

    actual = hfunc(u,v,families{ii},theta{ii});
    
    assertSize(testCase,actual,[1 1])
    assertTrue(testCase,actual<=0.99999 && actual>=0.00001)
    
end % ii

end


% hinv
function test_hinv(testCase)

u = 0.5;
v = 0.5;

% assert errors
assertError(testCase,@()hinv(u,v,'gumbel',0),'hinv:InvalidParameter')

% test 1: test all families
families = {'gumbel','clayton','frank','t','gauss','ind',...
            'amhaq','tawn','fgm','plackett','joe','surclayton',...
            'surgumbel','surjoe'};
theta = {2,2,2,[0.3,3],0.5,0,...
         0.5,0.5,0.5,2,2,2,...
         2,2};

for ii = 1:1:length(families)

    actual = hinv(u,v,families{ii},theta{ii});
    
    assertSize(testCase,actual,[1 1])
    assertTrue(testCase,actual<=0.99999 && actual>=0.00001)
    
end % ii

end


% iarray_rvine
function test_iarray_rvine(testCase)

% assert errors
A_err = [1 1 2; 0 2 1];
assertError(testCase,@()iarray_rvine(A_err),'iarray_rvine:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()iarray_rvine(A_err),'iarray_rvine:InvalidVineArray')

% test 1: simple vine array
A = 1;
actual = iarray_rvine(A);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2: 3-dimensional vine array
A = [1 1 2; 0 2 1; 0 0 3];
actual = iarray_rvine(A);
expected = [0 1 0;0 0 0;0 0 0];

assertEqual(testCase,actual,expected)

end


% llrvine
function test_llrvine(testCase)

A = [1 1 2; 0 2 1; 0 0 3];
family = {'gumbel', 'gumbel';'gumbel',0};
theta = {2,3;2,0};
n = 10;
u = simrvine(n,A,family,theta,0);

% assert errors
theta_err = {0,3;2,0};
assertError(testCase,@()llrvine(u,A,family,theta_err),'llrvine:InvalidParameter')

A_err = [1 1 2; 0 2 1];
assertError(testCase,@()llrvine(u,A_err,family,theta),'llrvine:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()llrvine(u,A_err,family,theta),'llrvine:InvalidVineArray')

% test 1: compute loglikelihood of a few points
[actualloglik,actualsll] = llrvine(u,A,family,theta);

assertSize(testCase,actualloglik,[n 1])
assertSize(testCase,actualsll,[1 1])
assertEqual(testCase,sum(actualloglik),actualsll);

end


% nktest
function test_nktest(testCase)

A = [1 1 2; 0 2 1; 0 0 3];
family = {'gumbel', 'gumbel';'gumbel',0};
theta = {2,3;2,0};
n = 10;
u = simrvine(n,A,family,theta,0);

% assert errors
assertError(testCase,@()nktest(u,0,family,theta,A,family,theta,A,family),'nktest:InvalidNumberOfArguments')

theta_err = {0,3;2,0};
assertError(testCase,@()nktest(u,0,family,theta_err,A,family,theta,A),'nktest:InvalidParameter')
assertError(testCase,@()nktest(u,0,family,theta,A,family,theta_err,A),'nktest:InvalidParameter')
assertError(testCase,@()nktest(u,0,family,theta,A,family,theta,A,family,theta_err,A),'nktest:InvalidParameter')

A_err = [1 1 2; 0 2 1];
assertError(testCase,@()nktest(u,0,family,theta,A_err,family,theta,A),'nktest:InvalidVineArray')
assertError(testCase,@()nktest(u,0,family,theta,A,family,theta,A_err),'nktest:InvalidVineArray')
assertError(testCase,@()nktest(u,0,family,theta,A,family,theta,A,family,theta,A_err),'nktest:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()nktest(u,0,family,theta,A_err,family,theta,A),'nktest:InvalidVineArray')
assertError(testCase,@()nktest(u,0,family,theta,A,family,theta,A_err),'nktest:InvalidVineArray')
assertError(testCase,@()nktest(u,0,family,theta,A,family,theta,A,family,theta,A_err),'nktest:InvalidVineArray')

assertError(testCase,@()nktest(u,2,family,theta,A,family,theta,A),'nktest:InvalidParallelizationMode')

end


% nopcount
function test_nopcount(testCase)

% test 1: empty family 
family = {};
actual = nopcount(family);
expected = 0;

assertEqual(testCase,actual,expected)

% test 2: one parameter copula
family = {'gumbel'};
actual = nopcount(family);
expected = 1;

assertEqual(testCase,actual,expected)

% test 3: two parameter copula
family = {'t'};
actual = nopcount(family);
expected = 2;

assertEqual(testCase,actual,expected)

% test 4: more than one copula
family = {'t', 'gumbel';'t',0};
actual = nopcount(family);
expected = 5;

assertEqual(testCase,actual,expected)

end


% pobs
function test_pobs(testCase)

x = [(1:1:9)' (1:1:9)'];

% test 1: simple check of correctness
actual = pobs(x);
expected = [(0.1:0.1:0.9)' (0.1:0.1:0.9)'];

assertEqual(testCase,actual,expected,'AbsTol',0.00000001)

end


% primmaxst
function test_primmaxst(testCase)

% assert errors
A_err = [1 1; 2 2];
assertError(testCase,@()primmaxst(A_err),'primmaxst:MatrixAsymmetric')

A_err = [-1 -1; -1 -1];
assertError(testCase,@()primmaxst(A_err),'primmaxst:InvalidMatrix')

% test 1: two node graph
A = [0 2; 2 0];
actual = primmaxst(A);
expected = [1 2];

assertEqual(testCase,actual,expected)

% test 2: three node graph
A = [0 2 1; 2 0 3; 1 3 0];
actual = primmaxst(A);
expected = [2 3; 2 1];

assertEqual(testCase,actual,expected)

end


% rvineselect
function test_rvineselect(testCase)

n = 10;
d = 3;

% assert errors
u = [0.5 0.5];
assertError(testCase,@()rvineselect(u),'rvineselect:InvalidNumberOfDimensions')

% test 1: estimation from a few points
u = rand(n,d);
[actualA,actualfamily,actualtheta] = rvineselect(u);

assertSize(testCase,actualA,[3 3])

assertClass(testCase,actualfamily,'cell')
assertSize(testCase,actualfamily,[2 2])

assertClass(testCase,actualtheta,'cell')
assertSize(testCase,actualtheta,[2 2])

end


% simrvine
function test_simrvine(testCase)

A = [1 1 2; 0 2 1; 0 0 3];
family = {'gumbel', 'gumbel';'gumbel',0};
theta = {2,3;2,0};
n = 10;

% assert errors
theta_err = {0,3;2,0};
assertError(testCase,@()simrvine(n,A,family,theta_err,0),'simrvine:InvalidParameter')

A_err = [1 1 2; 0 2 1];
assertError(testCase,@()simrvine(n,A_err,family,theta,0),'simrvine:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()simrvine(n,A_err,family,theta,0),'simrvine:InvalidVineArray')

assertError(testCase,@()simrvine(n,A,family,theta,2),'simrvine:InvalidParallelizationMode')

% test 1: simulate a few points
actual = simrvine(n,A,family,theta,0);

assertSize(testCase,actual,[n 3])
assertTrue(testCase,all(all(actual<=1 & actual>=0)))

end


% simrvinens
function test_simrvinens(testCase)

A = [1 1 2; 0 2 1; 0 0 3];
family = {'gumbel', 'gumbel';'gumbel',0};
theta = {2,3;'@(u)(thetalink((4*u(3) - 2).^2,''gumbel''))',0};
n = 10;

% assert errors
theta_err = {0,3;'@(u)(thetalink((4*u(3) - 2).^2,''gumbel''))',0};
assertError(testCase,@()simrvinens(n,A,family,theta_err,0),'simrvinens:InvalidParameter')

A_err = [1 1 2; 0 2 1];
assertError(testCase,@()simrvinens(n,A_err,family,theta,0),'simrvinens:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()simrvinens(n,A_err,family,theta,0),'simrvinens:InvalidVineArray')

assertError(testCase,@()simrvinens(n,A,family,theta,2),'simrvinens:InvalidParallelizationMode')

% test 1: simulate a few points
actual = simrvinens(n,A,family,theta,0);

assertSize(testCase,actual,[n 3])
assertTrue(testCase,all(all(actual<=1 & actual>=0)))

end


% simrvinetess
function test_simrvinetess(testCase)

A = [1 1 2; 0 2 1; 0 0 3];
family = {'gumbel', 'gumbel';{'gumbel' 'clayton' 'gumbel'},0};
theta = {2,3;{0.2 0.7 1; 2 3 4},0};
n = 10;

% assert errors
theta_err = {0,3;{0.2 0.7 1; 2 3 4},0};
assertError(testCase,@()simrvinetess(n,A,family,theta_err),'simrvinetess:InvalidParameter')

theta_err = {2,3;{0.2 0.7 1; 0 3 4},0};
assertError(testCase,@()simrvinetess(n,A,family,theta_err),'simrvinetess:InvalidParameter')

theta_err = {2,3;{0.2 0.7 0.8; 2 3 4},0};
assertError(testCase,@()simrvinetess(n,A,family,theta_err),'simrvinetess:InvalidTesselation')

A_err = [1 1 2; 0 2 1];
assertError(testCase,@()simrvinetess(n,A_err,family,theta),'simrvinetess:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()simrvinetess(n,A_err,family,theta),'simrvinetess:InvalidVineArray')

% test 1: simulate a few points
actual = simrvinetess(n,A,family,theta);

assertSize(testCase,actual,[n 3])
assertTrue(testCase,all(all(actual<=1 & actual>=0)))

end


% ssp
function test_ssp(testCase)

A = [1 1 2; 0 2 1; 0 0 3];
family = {'gumbel', 'gumbel';'gumbel',0};
theta = {4,5;6,0};
n = 10;
u = simrvine(n,A,family,theta,0);

% assert errors
A_err = [1 1 2; 0 2 1];
assertError(testCase,@()ssp(u,A_err,family),'ssp:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()ssp(u,A_err,family),'ssp:InvalidVineArray')

% test 1: estimation from a few points 
[actualtheta,actualloglik,actualsll,actualfamily] = ssp(u,A,family);

assertClass(testCase,actualtheta,'cell')
assertSize(testCase,actualtheta,[2 2])

assertClass(testCase,actualfamily,'cell')
assertSize(testCase,actualfamily,[2 2])

assertSize(testCase,actualloglik,[n 1])
assertSize(testCase,actualsll,[1 1])
assertEqual(testCase,sum(actualloglik),actualsll)

end


% transforma
function test_transforma(testCase)

% assert errors
A_err = [1 1 2; 0 2 1];
assertError(testCase,@()transforma(A_err),'transforma:InvalidVineArray')

A_err = [1 1 2; 0 1 1; 0 0 3];
assertError(testCase,@()transforma(A_err),'transforma:InvalidVineArray')

% test 1: 3-dimensional vine array
A = [3 3 2; 0 2 3; 0 0 1];
actual = transforma(A);
expected = [1 1 2; 0 2 1; 0 0 3];

assertEqual(testCase,actual,expected);

end



%% File fixtures  
function setupOnce(~)  
% unused        
end

function teardownOnce(~)  
% unused
end


%% Fresh fixtures  
function setup(~) 

rng(123456789);

end

function teardown(~)  
% unused
end