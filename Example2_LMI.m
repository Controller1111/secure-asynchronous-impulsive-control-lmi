clc
clear

n=3; N=4; m=6;
A=[-1 1 1;-100/7 0 0;9 0 -18/7];
D=[0;0;27/7];
G=[1 0 0];
B=eye(n);

c=0.3;
mu=0.9;
theta0=40;
delta1=0.005;
delta2=0.01;
barw=0.25;

% Generate one impulse period.
rng(2026,'twister');

Sequence=zeros(1,m);
randposition=randperm(m,2);
randagent=randi(N,1,2);
Sequence(randposition(1))=randagent(1);
Sequence(randposition(2))=randagent(2);

randseq=randperm(N);
q=1;

for z=1:m
    if Sequence(z)==0
        Sequence(z)=randseq(q);
        q=q+1;
    end
end

K=c*eye(n);
kronK=kron(eye(N),K);

L=[3 -1 -1 -1
   -1 3 -1 -1
   -1 -1 3 -1
   -1 -1 -1 3];

L1=[3 -1 -1 -1
    3 -1 -1 -1
    3 -1 -1 -1
    3 -1 -1 -1];

L2=[-1 3 -1 -1
    -1 3 -1 -1
    -1 3 -1 -1
    -1 3 -1 -1];

L3=[-1 -1 3 -1
    -1 -1 3 -1
    -1 -1 3 -1
    -1 -1 3 -1];

L4=[-1 -1 -1 3
    -1 -1 -1 3
    -1 -1 -1 3
    -1 -1 -1 3];

w=[1/4 1/4 1/4 1/4];

T1=diag([3/4,-1/4,-1/4,-1/4]);
T2=diag([-1/4,3/4,-1/4,-1/4]);
T3=diag([-1/4,-1/4,3/4,-1/4]);
T4=diag([-1/4,-1/4,-1/4,3/4]);

KT1=kron(T1,eye(n)); KT2=kron(T2,eye(n));
KT3=kron(T3,eye(n)); KT4=kron(T4,eye(n));
KL1=kron(L1,eye(n)); KL2=kron(L2,eye(n));
KL3=kron(L3,eye(n)); KL4=kron(L4,eye(n));

KT={KT1,KT2,KT3,KT4};
KL={KL1,KL2,KL3,KL4};

KD=kron(eye(N),D);
KbarG=kron(eye(N),G'*G);

U=[1 0 0
   0 1 0
   0 0 1
   -1 -1 -1];

KU=kron(U,eye(n));
KA=kron(eye(N),A);
KB=kron(eye(N),B);

disp('Impulse period:')
disp(Sequence)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the LMI variables for this period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setlmis([]);

Q=cell(1,m);
P=cell(1,m);
epsilon1=cell(1,m);
epsilon2=cell(1,m);

for j=1:m
    Q{j}=lmivar(1,[n*N,1]);
end

for j=1:m
    P{j}=lmivar(1,[n*N,1]);
end

for j=1:m
    epsilon1{j}=lmivar(1,[1,0]);
    epsilon2{j}=lmivar(1,[1,0]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous-time LMIs at delta1 and delta2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta=[delta1,delta2];

for j=1:m
    for r=1:2
        deltai=delta(r);

        CONT1(j,r)=newlmi;
        lmiterm([CONT1(j,r) 1 1 P{j}],(log(mu))/deltai+theta0,1);
        lmiterm([CONT1(j,r) 1 1 P{j}],1,KA,'s');
        lmiterm([CONT1(j,r) 1 1 Q{j}],1,-1/deltai);
        lmiterm([CONT1(j,r) 1 1 P{j}],1,1/deltai);
        lmiterm([CONT1(j,r) 1 1 epsilon1{j}],1,KbarG);
        lmiterm([CONT1(j,r) 1 1 epsilon2{j}],N*barw,KbarG);
        lmiterm([CONT1(j,r) 1 2 P{j}],1,KD);
        lmiterm([CONT1(j,r) 1 3 P{j}],1,KD);
        lmiterm([CONT1(j,r) 2 2 epsilon1{j}],-1,1);
        lmiterm([CONT1(j,r) 3 3 epsilon2{j}],-1,1);

        CONT2(j,r)=newlmi;
        lmiterm([CONT2(j,r) 1 1 Q{j}],(log(mu))/deltai+theta0,1);
        lmiterm([CONT2(j,r) 1 1 Q{j}],1,KA,'s');
        lmiterm([CONT2(j,r) 1 1 Q{j}],1,-1/deltai);
        lmiterm([CONT2(j,r) 1 1 P{j}],1,1/deltai);
        lmiterm([CONT2(j,r) 1 1 epsilon1{j}],1,KbarG);
        lmiterm([CONT2(j,r) 1 1 epsilon2{j}],N*barw,KbarG);
        lmiterm([CONT2(j,r) 1 2 Q{j}],1,KD);
        lmiterm([CONT2(j,r) 1 3 Q{j}],1,KD);
        lmiterm([CONT2(j,r) 2 2 epsilon1{j}],-1,1);
        lmiterm([CONT2(j,r) 3 3 epsilon2{j}],-1,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impulsive LMIs for the selected period.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LMIIMP=zeros(1,m);
CurrentSequence=Sequence;

for j=1:m
    Pin=P{j};

    if j<m
        Qout=Q{j+1};
    else
        Qout=Q{1};
    end

    agent=CurrentSequence(j);
    JumpMatrix=KU-c*KT{agent}*KL{agent}*KU;

    LMIIMP(j)=newlmi;
    lmiterm([LMIIMP(j) 1 1 Pin],-mu*KU',KU);
    lmiterm([LMIIMP(j) 1 2 Qout],JumpMatrix',1);
    lmiterm([LMIIMP(j) 2 2 Qout],-1,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Positive-definiteness constraints for all P_j and Q_j.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
positiveMargin=1;

for j=1:m
    POSP(j)=newlmi;
    lmiterm([POSP(j) 1 1 P{j}],-1,1);
    lmiterm([POSP(j) 1 1 0],positiveMargin*eye(n*N));

    POSQ(j)=newlmi;
    lmiterm([POSQ(j) 1 1 Q{j}],-1,1);
    lmiterm([POSQ(j) 1 1 0],positiveMargin*eye(n*N));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve one combined feasibility problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmisys=getlmis;
[tmin,xfeas]=feasp(lmisys,[1e-8,300,1e9,0,1]);

fprintf('Feasibility tmin = %.12g\n',tmin);

if tmin<0
    disp('The impulse period is feasible.')

    format short g

    vP=cell(1,m);
    vQ=cell(1,m);
    minEigP=zeros(1,m);
    minEigQ=zeros(1,m);

    for j=1:m
        vP{j}=dec2mat(lmisys,xfeas,P{j});
        vQ{j}=dec2mat(lmisys,xfeas,Q{j});
        minEigP(j)=min(eig(vP{j}));
        minEigQ(j)=min(eig(vQ{j}));

        fprintf('\nP%d =\n',j-1);
        disp(vP{j});
        fprintf('Q%d =\n',j-1);
        disp(vQ{j});
    end

else
    disp('The impulse period is infeasible.')
end
