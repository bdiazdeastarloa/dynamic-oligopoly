function plotpol(par,sol)

% =========================================================================
% Dynamic oligopoly.
% Continuous time version.
%
% plotpol
%
% Plot policy functions for different industry configurations. 
% 'par' is the file containing parameter values.
% 'sol' is the file containing the solution to be plotted.
%
% Written by Bernardo Diaz de Astarloa @ PSU 2015 based on code by Uli
% Doraszleski. (Optimized for speed.)
% =========================================================================

load(sol);
load(par);
close all

% Assign parameters.
binom = double(par.binom); 
N     = par.N;
M     = par.M;
D     = par.D;
S     = par.S;
mgc   = par.mgc;

% Set 3D plots on(1)/off(0);
plot3d = 1;

% Round values to zero if extremely low.
err = 1e-25;
y1(y1<err) = 0.0;      %#ok
V1(V1<err) = 0.0;      %#ok

%% Plot value functions and policy for different industry configurations.

% List of states for plots.
% Policies of firm in row 1: for every common state, for every number of 
% competitors, for all competitors at different productivity levels.
indl = zeros(D,N,M);
indm = indl;
indh = indl;
indel = indl;
indem = indl;
indeh = indl;
for d=1:D
    for Nstar=1:N
        % Low: competitors at lowest productivity level.
        dil = [d*ones(1,M);(1:M);ones(Nstar-1,M);M+1*ones(N-Nstar,M)];
        % Medium: competitors at medium productivity level.
        dim = [d*ones(1,M);(1:M);round(M/2)*ones(Nstar-1,M);(M+1)*ones(N-Nstar,M)];
        % High: competitors at highest productivity level.
        dih = [d*ones(1,M);(1:M);M*ones(Nstar-1,M);(M+1)*ones(N-Nstar,M)];
        % Entry at low.
        diel = [d*ones(1,M);(M+1)*ones(1,M);(1:M);ones(Nstar-1,M);(M+1)*ones(N-Nstar-1,M)];
        % Entry at medium.
        diem = [d*ones(1,M);(M+1)*ones(1,M);(1:M);round(M/2)*ones(Nstar-1,M);(M+1)*ones(N-Nstar-1,M)];
        % Entry at high.
        dieh = [d*ones(1,M);(M+1)*ones(1,M);(1:M);M*ones(Nstar-1,M);(M+1)*ones(N-Nstar-1,M)];
        % Encode all industry configurations and save linear indices to recover policies.
        [code,pos] = allencode(dil,binom,N,D); indl(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(dim,binom,N,D); indm(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(dih,binom,N,D); indh(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(diel,binom,N,D); indel(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(diem,binom,N,D); indem(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
        [code,pos] = allencode(dieh,binom,N,D); indeh(d,Nstar,:) = sub2ind([N,S],pos(1,:),code);
    end
end

% 3D graphs for pairs of firms.
ind3 = zeros(D,M^2);
inde3 = ind3;
for d=1:D
    di = [d*ones(1,M.^2);cartx((1:M),(1:M))';(M+1)*ones(N-2,M.^2)];
    die = [d*ones(1,M.^2);(M+1)*ones(1,M.^2);cartx((1:M),(1:M))';(M+1)*ones(N-3,M.^2)];
    [code,pos] = allencode(di,binom,N,D); ind3(d,:) = sub2ind([N,S],pos(1,:),code);
    [code,pos] = allencode(die,binom,N,D); inde3(d,:) = sub2ind([N,S],pos(1,:),code);
end


%% Do plots.
%
for d=1:D
    for Nstar=1:N
        figure('visible','off');        
        if plot3d && (Nstar==2)
            % 3D graphs (only if there are 2 firms).
            % Value function (incumbent).
            subplot(2,4,1);
            mesh((1:M),(1:M),reshape(V1(ind3(d,:)),M,M)');
            v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
            xlabel('w_1');
            ylabel('w_2');
            zlabel('V(w,c)');
            % Price (incumbent).
            subplot(2,4,2);
            mesh((1:M),(1:M),reshape(p1(ind3(d,:)),M,M)');
            v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; axis(v);
            xlabel('w_1');
            ylabel('w_2');
            zlabel('p(w,c)');
            % Investment (incumbent).
            subplot(2,4,3);
            mesh((1:M),(1:M),reshape(x1(ind3(d,:)),M,M)');
            v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
            xlabel('w_1');
            ylabel('w_2');
            zlabel('x(w,c)');
            % Exit (incumbent).
            subplot(2,4,4);
            mesh((1:M),(1:M),reshape(y1(ind3(d,:)),M,M)');
            v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
            xlabel('w_1');
            ylabel('w_2');
            zlabel('y(w,c)');
            if Nstar<N
                % Entrants.
                % Axes are for states of two incumbents.
                % Value of entry.
                subplot(2,4,5);
                mesh((1:M),(1:M),reshape(V1(inde3(d,:)),M,M)');
                v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
                xlabel('w_2');
                ylabel('w_3');
                zlabel('V^e(w,c)');
                % Price of potential entrant (zeros).
                subplot(2,4,6);
                mesh((1:M),(1:M),reshape(p1(inde3(d,:)),M,M)');
                v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
                xlabel('w_2');
                ylabel('w_3');
                zlabel('p^e(w,c)');
                % Investment of potential entrant (zeros).
                subplot(2,4,7);
                mesh((1:M),(1:M),reshape(x1(inde3(d,:)),M,M)');
                v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
                xlabel('w_2');
                ylabel('w_3');
                zlabel('x^e(w,c)');
                % Entry policy.
                subplot(2,4,8);
                mesh((1:M),(1:M),reshape(y1(inde3(d,:)),M,M)');
                v = axis; v(1) = 1; v(2) = M; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
                xlabel('w_2');
                ylabel('w_3');
                zlabel('y^e(w,c)');
            end
        else
            % Single firm plots.
            % Value function.
            subplot(2,4,1);
            plot((1:M),squeeze(V1(indl(d,Nstar,:))),'r-.',(1:M),squeeze(V1(indm(d,Nstar,:))),'g--',(1:M),squeeze(V1(indh(d,Nstar,:))),'b-');
            v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
            xlabel('w_1');
            ylabel('V(w,c)');
            legend(sprintf('w_n=%d',1),sprintf('w_n=%d',round(M./2)),sprintf('w_n=%d',M),'Location','southeast');
            % Price.
            subplot(2,4,2);
            plot((1:M),squeeze(p1(indl(d,Nstar,:))),'r-.',(1:M),squeeze(p1(indm(d,Nstar,:))),'g--',(1:M),squeeze(p1(indh(d,Nstar,:))),'b-');
            v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
            xlabel('w_1');
            ylabel('p(w,c)');
            legend(sprintf('w_n=%d',1),sprintf('w_n=%d',round(M./2)),sprintf('w_n=%d',M),'Location','northeast');
            % Investment.
            subplot(2,4,3);
            plot((1:M),squeeze(x1(indl(d,Nstar,:))),'r-.',(1:M),squeeze(x1(indm(d,Nstar,:))),'g--',(1:M),squeeze(x1(indh(d,Nstar,:))),'b-');
            v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
            xlabel('w_1');
            ylabel('x(w,c)');
            legend(sprintf('w_n=%d',1),sprintf('w_n=%d',round(M./2)),sprintf('w_n=%d',M),'Location','northeast'); 
            % Exit.
            subplot(2,4,4);
            plot((1:M),squeeze(y1(indl(d,Nstar,:))),'r-.',(1:M),squeeze(y1(indm(d,Nstar,:))),'g--',(1:M),squeeze(y1(indh(d,Nstar,:))),'b-');
            v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
            xlabel('w_1');
            ylabel('y(w,c)');
            legend(sprintf('w_n=%d',1),sprintf('w_n=%d',round(M./2)),sprintf('w_n=%d',M),0);
            if Nstar<N
                % Entrants.
                % Value of entry.
                subplot(2,4,5);
                plot((1:M),squeeze(V1(indel(d,Nstar,:))),'r-.',(1:M),squeeze(V1(indem(d,Nstar,:))),'g--',(1:M),squeeze(V1(indeh(d,Nstar,:))),'b-');
                v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
                xlabel('w_2');
                ylabel('V^{e}(w,c)');
                legend(sprintf('w_1=%d, w_n=%d',M+1,1),sprintf('w_1=%d, w_n=%d',M+1,round(M./2)),sprintf('w_1=%d, w_n=%d',M+1,M),'Location','southeast');
                % Price.
                subplot(2,4,6);
                plot((1:M),squeeze(p1(indel(d,Nstar,:))),'r-.',(1:M),squeeze(p1(indem(d,Nstar,:))),'g--',(1:M),squeeze(p1(indeh(d,Nstar,:))),'b-');
                v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
                xlabel('w_2');
                ylabel('x^{e}(w,c)');
                legend(sprintf('w_1=%d, w_n=%d',M+1,1),sprintf('w_1=%d, w_n=%d',M+1,round(M./2)),sprintf('w_1=%d, w_n=%d',M+1,M),'Location','northeast');
                % Investment.
                subplot(2,4,7);
                plot((1:M),squeeze(x1(indel(d,Nstar,:))),'r-.',(1:M),squeeze(x1(indem(d,Nstar,:))),'g--',(1:M),squeeze(x1(indeh(d,Nstar,:))),'b-');
                v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
                xlabel('w_2');
                ylabel('x^{e}(w,c)');
                legend(sprintf('w_1=%d, w_n=%d',M+1,1),sprintf('w_1=%d, w_n=%d',M+1,round(M./2)),sprintf('w_1=%d, w_n=%d',M+1,M),'Location','northeast');
                % Entry.
                subplot(2,4,8);
                plot((1:M),squeeze(y1(indel(d,Nstar,:))),'r-.',(1:M),squeeze(y1(indem(d,Nstar,:))),'g--',(1:M),squeeze(y1(indeh(d,Nstar,:))),'b-');
                v = axis; v(1) = 1; v(2) = M; v(3) = 0; axis(v);
                xlabel('w_2');
                ylabel('y^{e}(w,c)');
                legend(sprintf('w_1=%d, w_n=%d',M+1,1),sprintf('w_1=%d, w_n=%d',M+1,round(M./2)),sprintf('w_1=%d, w_n=%d',M+1,M),'Location','southeast');
            end
        end
        set(gcf,'PaperOrientation','landscape');
        set(gcf,'PaperUnits','normalized');
        set(gcf,'PaperPosition', [0 0 1 1]);
        figname = sprintf('figures/ind_N%dM%dD%d-%d-%d.pdf',[N,M,D,d,Nstar]);
        print(gcf, '-dpdf', figname)
    end
end
%}

clear indl indm indh indel indem indeh ind3 inde3

%% Plot specific indicators.

% States: one firm at wn, given Nstar competitors all at the
% same productivity level, for all productivities, for all common states.
indo = zeros(M,D*M,N-1);
for wn=1:M
    for Nstar=1:N-1
        cart = cartx((1:D),(1:M))';
        dipc = [cart(1,:); wn*ones(1,D*M); repmat(cart(2,:),Nstar,1);(M+1)*ones(N-Nstar-1,D*M)];
        [code,pos] = allencode(dipc,binom,N,D); indo(wn,:,Nstar) = sub2ind([N,S],pos(1,:),code);
    end
end

% Price/Cost margins.
for Nstar=1:N-1
    % Low prod.
    w = 1;
    pc = reshape(p1(indo(w,:,Nstar)),M,D)./mgc(w);
    figure('visible','off');
    mesh((1:D),(1:M),pc);
    v = axis; v(1) = 1; v(2) = D; v(3) = 1; v(4) = M; v(5) = 0.2; axis(v);
    xlabel('Common state');
    ylabel('Competitors'' productivity');
    zlabel('Price/Cost');   
    view([138 30]);
    set(get(gca,'xlabel'),'rotation',23,'Position',get(get(gca,'xlabel'),'Position') - [0 0.5 0]);
    set(get(gca,'ylabel'),'rotation',-28,'Position',get(get(gca,'ylabel'),'Position') - [0.5 0 0]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition', [0.5 0.5 13.5 13.5]);
    set(gcf, 'PaperSize', [14 14]); 
    figname = sprintf('figures/pcmgn_N%dM%dD%d-%d-%d.pdf',[N,M,D,w,Nstar]);
    saveas(gcf,figname)
    
    % Medium productivity.
    w = round(M/2);
    pc = reshape(p1(indo(w,:,Nstar)),M,D)./mgc(w);
    figure('visible','off');
    mesh((1:D),(1:M),pc);
    v = axis; v(1) = 1; v(2) = D; v(3) = 1; v(4) = M; v(5) = 0.2; axis(v);
    xlabel('Common state');
    ylabel('Competitors'' productivity');
    zlabel('Price/Cost');    
    view([138 30]);
    set(get(gca,'xlabel'),'rotation',23,'Position',get(get(gca,'xlabel'),'Position') - [0 0.5 0]);
    set(get(gca,'ylabel'),'rotation',-28,'Position',get(get(gca,'ylabel'),'Position') - [0.5 0 0]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition', [0.5 0.5 13.5 13.5]);
    set(gcf, 'PaperSize', [14 14]); 
    figname = sprintf('figures/pcmgn_N%dM%dD%d-%d-%d.pdf',[N,M,D,w,Nstar]);
    saveas(gcf,figname)
    
    % High productivity.
    w = M;
    pc = reshape(p1(indo(w,:,Nstar)),M,D)./mgc(w);
    figure('visible','off');
    mesh((1:D),(1:M),pc);
    v = axis; v(1) = 1; v(2) = D; v(3) = 1; v(4) = M; v(5) = 0.2; axis(v);
    xlabel('Common state');
    ylabel('Competitors'' productivity');
    zlabel('Price/Cost');     
    view([138 30]);
    set(get(gca,'xlabel'),'rotation',23,'Position',get(get(gca,'xlabel'),'Position') - [0 0.5 0]);
    set(get(gca,'ylabel'),'rotation',-28,'Position',get(get(gca,'ylabel'),'Position') - [0.5 0 0]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition', [0.5 0.5 13.5 13.5]);
    set(gcf, 'PaperSize', [14 14]); 
    figname = sprintf('figures/pcmgn_N%dM%dD%d-%d-%d.pdf',[N,M,D,w,Nstar]);
    saveas(gcf,figname)
end



% Investment.
for Nstar=1:N-1
    % Low prod.
    w = 1;
    pc = reshape(x1(indo(w,:,Nstar)),M,D);
    figure('visible','off');
    mesh((1:D),(1:M),pc);
    v = axis; v(1) = 1; v(2) = D; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
    xlabel('Common state');
    ylabel('Competitors'' productivity');
    zlabel('Investment');   
    view([138 30]);
    set(get(gca,'xlabel'),'rotation',23,'Position',get(get(gca,'xlabel'),'Position') - [0 0.5 0]);
    set(get(gca,'ylabel'),'rotation',-28,'Position',get(get(gca,'ylabel'),'Position') - [0.5 0 0]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition', [0.5 0.5 13.5 13.5]);
    set(gcf, 'PaperSize', [14 14]); 
    figname = sprintf('figures/rdx_N%dM%dD%d-%d-%d.pdf',[N,M,D,w,Nstar]);
    saveas(gcf,figname)
    
    % Medium productivity.
    w = round(M/2);
    pc = reshape(x1(indo(w,:,Nstar)),M,D);
    figure('visible','off');
    mesh((1:D),(1:M),pc);
    v = axis; v(1) = 1; v(2) = D; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
    xlabel('Common state');
    ylabel('Competitors'' productivity');
    zlabel('Investment');    
    view([138 30]);
    set(get(gca,'xlabel'),'rotation',23,'Position',get(get(gca,'xlabel'),'Position') - [0 0.5 0]);
    set(get(gca,'ylabel'),'rotation',-28,'Position',get(get(gca,'ylabel'),'Position') - [0.5 0 0]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition', [0.5 0.5 13.5 13.5]);
    set(gcf, 'PaperSize', [14 14]); 
    figname = sprintf('figures/rdx_N%dM%dD%d-%d-%d.pdf',[N,M,D,w,Nstar]);
    saveas(gcf,figname)
    
    % High productivity.
    w = M;
    pc = reshape(x1(indo(w,:,Nstar)),M,D);
    figure('visible','off');
    mesh((1:D),(1:M),pc);
    v = axis; v(1) = 1; v(2) = D; v(3) = 1; v(4) = M; v(5) = 0; axis(v);
    xlabel('Common state');
    ylabel('Competitors'' productivity');
    zlabel('Investment');     
    view([138 30]);
    set(get(gca,'xlabel'),'rotation',23,'Position',get(get(gca,'xlabel'),'Position') - [0 0.5 0]);
    set(get(gca,'ylabel'),'rotation',-28,'Position',get(get(gca,'ylabel'),'Position') - [0.5 0 0]);
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition', [0.5 0.5 13.5 13.5]);
    set(gcf, 'PaperSize', [14 14]); 
    figname = sprintf('figures/rdx_N%dM%dD%d-%d-%d.pdf',[N,M,D,w,Nstar]);
    saveas(gcf,figname)
end

close all;