clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common settings of Example 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2;
N=6;
m=6;

A=[0 -0.5;0.5 0];
D=[-0.5 0;0 0.5];
G=1.5*eye(n);

delta1=0.005;
delta2=0.01;
p=0.1;

% One complete asynchronous period. A cyclic shift gives an equivalent
% description because Q_m=Q_0.
sigma=[1 2 3 4 5 6];

% Parameters of Example 1, Case 1
L=[0 0 0 0 0 0
   -1 1 0 0 0 0
   0 -1 1 0 0 0
   0 0 -1 1 0 0
   0 0 0 -1 1 0
   0 0 0 0 -1 1];

w=[1 0 0 0 0 0];

U=[0 0 0 0 0
   1 0 0 0 0
   0 1 0 0 0
   0 0 1 0 0
   0 0 0 1 0
   0 0 0 0 1];

c=0.990;
theta0=16.28;
mu=0.9;

    barw=sum(w.^2);
    KU=kron(U,eye(n));
    KA=kron(eye(N),A);
    KD=kron(eye(N),D);
    KbarG=kron(eye(N),G'*G);

    KT=cell(1,N);
    KL=cell(1,N);

    for agent=1:N
        t=-w(agent)*ones(N,1);
        t(agent)=1-w(agent);
        Tagent=diag(t);

        Lagent=ones(N,1)*L(agent,:);
        KT{agent}=kron(Tagent,eye(n));
        KL{agent}=kron(Lagent,eye(n));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define the LMI variables for the current case.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    setlmis([]);

    P=cell(1,m);
    Q=cell(1,m);
    epsilon1=cell(1,m);
    epsilon2=cell(1,m);

    for j=1:m
        P{j}=lmivar(1,[n*N,1]);
        Q{j}=lmivar(1,[n*N,1]);
        epsilon1{j}=lmivar(1,[1,0]);
        epsilon2{j}=lmivar(1,[1,0]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Continuous-time LMIs in (17), for delta1 and delta2 and for both
    % H=P_j and H=Q_j.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta=[delta1,delta2];

    for j=1:m
        for r=1:2
            deltai=delta(r);

            CONT1(j,r)=newlmi;
            lmiterm([CONT1(j,r) 1 1 P{j}], ...
                (log(mu))/deltai+theta0,1);
            lmiterm([CONT1(j,r) 1 1 P{j}],1,KA,'s');
            lmiterm([CONT1(j,r) 1 1 Q{j}],1,-1/deltai);
            lmiterm([CONT1(j,r) 1 1 P{j}],1,1/deltai);
            lmiterm([CONT1(j,r) 1 1 epsilon1{j}],1,KbarG);
            lmiterm([CONT1(j,r) 1 1 epsilon2{j}],N*barw,KbarG);
            lmiterm([CONT1(j,r) 1 2 P{j}],1,KD);
            lmiterm([CONT1(j,r) 1 3 P{j}],1,KD);
            lmiterm([CONT1(j,r) 2 2 epsilon1{j}],-1,eye(N*n));
            lmiterm([CONT1(j,r) 3 3 epsilon2{j}],-1,eye(N*n));

            CONT2(j,r)=newlmi;
            lmiterm([CONT2(j,r) 1 1 Q{j}], ...
                (log(mu))/deltai+theta0,1);
            lmiterm([CONT2(j,r) 1 1 Q{j}],1,KA,'s');
            lmiterm([CONT2(j,r) 1 1 Q{j}],1,-1/deltai);
            lmiterm([CONT2(j,r) 1 1 P{j}],1,1/deltai);
            lmiterm([CONT2(j,r) 1 1 epsilon1{j}],1,KbarG);
            lmiterm([CONT2(j,r) 1 1 epsilon2{j}],N*barw,KbarG);
            lmiterm([CONT2(j,r) 1 2 Q{j}],1,KD);
            lmiterm([CONT2(j,r) 1 3 Q{j}],1,KD);
            lmiterm([CONT2(j,r) 2 2 epsilon1{j}],-1,eye(N*n));
            lmiterm([CONT2(j,r) 3 3 epsilon2{j}],-1,eye(N*n));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Six impulsive LMIs in Corollary 1, one for each agent in sigma.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:m
        if j<m
            Qout=Q{j+1};
        else
            Qout=Q{1};
        end

        agent=sigma(j);
        Pi=(1-p)*(KU-c*KT{agent}*KL{agent}*KU);

        LMIIMP(j)=newlmi;
        lmiterm([LMIIMP(j) 1 1 P{j}],-mu*KU',KU);
        lmiterm([LMIIMP(j) 1 2 Qout],p*KU',1);
        lmiterm([LMIIMP(j) 1 3 Qout],Pi',1);
        lmiterm([LMIIMP(j) 2 2 Qout],-p,1);
        lmiterm([LMIIMP(j) 3 3 Qout],-(1-p),1);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Strict positive-definiteness conditions for all P_j and Q_j.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    positiveMargin=1;

    for j=1:m
        POSP(j)=newlmi;
        lmiterm([POSP(j) 1 1 P{j}],-1,1);
        lmiterm([POSP(j) 1 1 0],positiveMargin*eye(N*n));

        POSQ(j)=newlmi;
        lmiterm([POSQ(j) 1 1 Q{j}],-1,1);
        lmiterm([POSQ(j) 1 1 0],positiveMargin*eye(N*n));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve and extract the complete matrices.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lmisys=getlmis;
    [tmin,xfeas]=feasp(lmisys,[1e-8,300,1e9,0,1]);

    fprintf('tmin = %.12g\n',tmin);

    if tmin<0
        disp('The LMIs are feasible.')

        for j=1:m
            vP{j}=dec2mat(lmisys,xfeas,P{j});
            vQ{j}=dec2mat(lmisys,xfeas,Q{j});

            fprintf('\nP%d =\n',j-1);
            disp(vP{j});

            fprintf('Q%d =\n',j-1);
            disp(vQ{j});
        end
    else
        disp('The LMIs are infeasible.')
    end
