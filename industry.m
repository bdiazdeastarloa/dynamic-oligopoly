%
% MPE.
% Continuous-time industry equilibrium.
% Run MEX/EXE file.
%

clear all;
delete('industry.lst');
diary('industry.lst');

% Start timer.
tic;

% Setup.
run setup;

% Program control.
plot3d = 1;

% Starting values: Value and policy functions.
V0 = pi./rho;
x0 = zeros(N,S);
y0 = zeros(N,S);

% Convert to unsigned 32-bit integer (ulong).
binom = uint32(binom);

% Convert to unsigned 8-bit integer (uchar).
state = uint8(state);

% Save variables and clear them from workspace.
filestr = sprintf('indP%dD%dN%dM%d.mat',[pro,D,N,M]);
eval(['save ',filestr,vars]);
eval(['clear ',vars]);

% Run MEX file.
compMPE(filestr);

% Run EXE file.
% eval(['!compMPE ',filestr]);

% Load variables.
load(filestr);

% Convert to double.
binom = double(binom);
state = double(state);
% eval(['save ',filestr,' binom state -append']);

% List of states for plots.
for d=1:D
    for Nstar=1:N
        dil = [repmat(d,1,M);[1:M];repmat(1,Nstar-1,M);repmat(M+1,N-Nstar,M)];
        dim = [repmat(d,1,M);[1:M];repmat(round(M./2),Nstar-1,M);repmat(M+1,N-Nstar,M)];
        dih = [repmat(d,1,M);[1:M];repmat(M,Nstar-1,M);repmat(M+1,N-Nstar,M)];
        diel = [repmat(d,1,M);repmat(M+1,1,M);[1:M];repmat(1,Nstar-1,M);repmat(M+1,N-Nstar-1,M)];
        diem = [repmat(d,1,M);repmat(M+1,1,M);[1:M];repmat(round(M./2),Nstar-1,M);repmat(M+1,N-Nstar-1,M)];
        dieh = [repmat(d,1,M);repmat(M+1,1,M);[1:M];repmat(M,Nstar-1,M);repmat(M+1,N-Nstar-1,M)];
        [code,pos] = allencode(dil); indl(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(dim); indm(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(dih); indh(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(diel); indel(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(diem); indem(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(dieh); indeh(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
    end
end
for d=1:D
    di = [repmat(d,1,M.^2);cartcols(repmat([1:M]',1,2))';repmat(M+1,N-2,M.^2)];
    die = [repmat(d,1,M.^2);repmat(M+1,1,M.^2);cartcols(repmat([1:M]',1,2))';repmat(M+1,N-3,M.^2)];
    [code,pos] = allencode(di); ind3(d,:) = sub2ind([N,S],pos(1,:),code);
    [code,pos] = allencode(die); inde3(d,:) = sub2ind([N,S],pos(1,:),code);
end

% Prepare plots.
for d=1:D
    for Nstar=1:N
        f((d-1).*N+Nstar) = figure;
    end
end

% Plots.
for d=1:D
    for Nstar=1:N
        figure(f((d-1).*N+Nstar));
        if plot3d & (Nstar==2)
            subplot(2,3,1);
            mesh([1:M],[1:M],reshape(V1(ind3(d,:)),M,M)');
            v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
            xlabel('i_1');
            ylabel('i_2');
            zlabel('V(d,i)');
            subplot(2,3,2);
            mesh([1:M],[1:M],reshape(x1(ind3(d,:)),M,M)');
            v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
            xlabel('i_1');
            ylabel('i_2');
            zlabel('x(d,i)');
            subplot(2,3,3);
            mesh([1:M],[1:M],reshape(y1(ind3(d,:)),M,M)');
            v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
            xlabel('i_1');
            ylabel('i_2');
            zlabel('y(d,i)');
            if Nstar<N
                subplot(2,3,4);
                mesh([1:M],[1:M],reshape(V1(inde3(d,:)),M,M)');
                v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
                xlabel('i_2');
                ylabel('i_3');
                zlabel('V^e(d,i)');
                subplot(2,3,5);
                mesh([1:M],[1:M],reshape(x1(inde3(d,:)),M,M)');
                v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
                xlabel('i_2');
                ylabel('i_3');
                zlabel('x^e(d,i)');
                subplot(2,3,6);
                mesh([1:M],[1:M],reshape(y1(inde3(d,:)),M,M)');
                v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
                xlabel('i_2');
                ylabel('i_3');
                zlabel('y^e(d,i)');
            end
        else
            subplot(2,3,1);
            plot([1:M],squeeze(V1(indl(d,Nstar,:))),'r-.',[1:M],squeeze(V1(indm(d,Nstar,:))),'g--',[1:M],squeeze(V1(indh(d,Nstar,:))),'b-');
            v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
            xlabel('i_1');
            ylabel('V(d,i)');
            legend(sprintf('i_n=%d',1),sprintf('i_n=%d',round(M./2)),sprintf('i_n=%d',M),0);
            subplot(2,3,2);
            plot([1:M],squeeze(x1(indl(d,Nstar,:))),'r-.',[1:M],squeeze(x1(indm(d,Nstar,:))),'g--',[1:M],squeeze(x1(indh(d,Nstar,:))),'b-');
            v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
            xlabel('i_1');
            ylabel('x(d,i)');
            legend(sprintf('i_n=%d',1),sprintf('i_n=%d',round(M./2)),sprintf('i_n=%d',M),0);
            subplot(2,3,3);
            plot([1:M],squeeze(y1(indl(d,Nstar,:))),'r-.',[1:M],squeeze(y1(indm(d,Nstar,:))),'g--',[1:M],squeeze(y1(indh(d,Nstar,:))),'b-');
            v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
            xlabel('i_1');
            ylabel('y(d,i)');
            legend(sprintf('i_n=%d',1),sprintf('i_n=%d',round(M./2)),sprintf('i_n=%d',M),0);
            if Nstar<N
                subplot(2,3,4);
                plot([1:M],squeeze(V1(indel(d,Nstar,:))),'r-.',[1:M],squeeze(V1(indem(d,Nstar,:))),'g--',[1:M],squeeze(V1(indeh(d,Nstar,:))),'b-');
                v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
                xlabel('i_2');
                ylabel('V^{e}(d,i)');
                legend(sprintf('i_1=%d, i_n=%d',M+1,1),sprintf('i_1=%d, i_n=%d',M+1,round(M./2)),sprintf('i_1=%d, i_n=%d',M+1,M),0);
                subplot(2,3,5);
                plot([1:M],squeeze(x1(indel(d,Nstar,:))),'r-.',[1:M],squeeze(x1(indem(d,Nstar,:))),'g--',[1:M],squeeze(x1(indeh(d,Nstar,:))),'b-');
                v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
                xlabel('i_2');
                ylabel('x^{e}(d,i)');
                legend(sprintf('i_1=%d, i_n=%d',M+1,1),sprintf('i_1=%d, i_n=%d',M+1,round(M./2)),sprintf('i_1=%d, i_n=%d',M+1,M),0);
                subplot(2,3,6);
                plot([1:M],squeeze(y1(indel(d,Nstar,:))),'r-.',[1:M],squeeze(y1(indem(d,Nstar,:))),'g--',[1:M],squeeze(y1(indeh(d,Nstar,:))),'b-');
                v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
                xlabel('i_2');
                ylabel('y^{e}(d,i)');
                legend(sprintf('i_1=%d, i_n=%d',M+1,1),sprintf('i_1=%d, i_n=%d',M+1,round(M./2)),sprintf('i_1=%d, i_n=%d',M+1,M),0);
            end
        end
        subtitle(sprintf('Continuous-Time Industry Equilibrium (d=%d, D=%d, N^{*}=%d, N=%d, M=%d)\n(done=%d, iter=%d, tol(V, x, y)=(%g, %g, %g))',[d,D,Nstar,N,M,info]));
    end
end

% Print plots.
for d=1:D
    for Nstar=1:N
        figure(f((d-1).*N+Nstar));
        filestr = sprintf('indP%dD%dN%dM%d-%d-%d.ps',[pro,D,N,M,d,Nstar]);
        print('-dpsc',filestr);
    end
end

% Stop timer.
elapsed = toc./60;

% Print results.
disp(' ');
disp('-------');
disp('Summary');
disp('-------');
disp(' ');
disp(sprintf('Elapsed time (min.)  :  %12.4f',elapsed));

pause;
close all;
diary off;
