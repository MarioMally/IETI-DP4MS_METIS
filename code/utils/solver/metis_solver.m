function [sol,lambda,condEst,iter] = ieti_dp_solver(mp_strct,tol)
    
    %% Prepare splitting of local DOFs
    dirCell = {mp_strct.patch_arr.dir_dofs};
    priCell = {mp_strct.patch_arr.pri_dofs};

    depCell = {mp_strct.patch_arr.dep_dofs};
    idepCell = {mp_strct.patch_arr.idep_dofs};
    intfCell = cellfun(@(dep,idep) union(dep,idep), depCell, idepCell, 'UniformOutput', false);

    remCell = cellfun(@(sp,t) setdiff(1:sp.ndof,t),...
        {mp_strct.patch_arr.space}, {mp_strct.patch_arr.tree}, 'UniformOutput', false);
    % Remove dirichlet from cotree
    remCell = cellfun(@(rem,dir) setdiff(rem,dir), remCell, dirCell, 'UniformOutput', false);
    % Remove primal dofs from cotree
    remCell = cellfun(@(rem,pri) setdiff(rem,pri), remCell, priCell, 'UniformOutput', false);

    [~,remIntfCell,~] = cellfun(@(rem,int) intersect(rem,int), remCell, intfCell, 'UniformOutput', false);
    [~,remVolCell] = cellfun(@(rem,int) setdiff(rem,int), remCell, intfCell, 'UniformOutput', false);

    
    %% Prepare blocks for stiffness matrices
    ArrCell = cellfun(@(A,rem) A(rem,rem), {mp_strct.patch_arr.Ai}, remCell, 'UniformOutput', false);
    ArdCell = cellfun(@(A,rem,dir) A(rem,dir), {mp_strct.patch_arr.Ai}, remCell, dirCell, 'UniformOutput', false);
    ArpCell = cellfun(@(A,rem,pri) A(rem,pri), {mp_strct.patch_arr.Ai}, remCell, priCell, 'UniformOutput', false);
    ApdCell = cellfun(@(A,pri,dir) A(pri,dir), {mp_strct.patch_arr.Ai}, priCell, dirCell, 'UniformOutput', false);
    AppCell = cellfun(@(A,pri) A(pri,pri), {mp_strct.patch_arr.Ai}, priCell, 'UniformOutput', false);

    %% Prepare blocks for constraint matrices
    BrCell = cellfun(@(B,rem) B(:,rem), {mp_strct.patch_arr.Bi}, remCell, 'UniformOutput', false);
    Br = [BrCell{:}];
    mult2keep = sum(abs(Br),2)~=0;
    BrCell = cellfun(@(B,rem) B(mult2keep,rem), {mp_strct.patch_arr.Bi}, remCell, 'UniformOutput', false);
    
    %% Prepare Cp blocks
    CpCell = {mp_strct.patch_arr.Cp};

    %% Compute dirichlet contribution
    adCell = {mp_strct.patch_arr.gDi};

    %% Prepare blocks for rhs
    frCell = cellfun(@(f,rem,Ard,ad) sum([f(rem);-(Ard*ad)'],1,"omitmissing"),...
        {mp_strct.patch_arr.fi}, remCell, ArdCell, adCell, 'UniformOutput', false);
    fpCell = cellfun(@(f,pri,Apd,ad) sum([f(pri);-(Apd*ad)'],1,"omitmissing"),...
        {mp_strct.patch_arr.fi}, priCell, ApdCell, adCell, 'UniformOutput', false);

    %% Solve system after first Schur complement
    tic;
    decArrCell = cellfun(@(A) decomposition(A,"chol","upper"), ArrCell, 'UniformOutput', false);
    t1 = toc;

    fprintf('\t\tLocal factorization of all Arr^(i) (seq.) took %fs\n',t1);
    
    tic;
    F = sparse(size(CpCell{1},2),size(CpCell{1},2));
    G = sparse(size(BrCell{1},1),size(CpCell{1},2));
    W = sparse(size(BrCell{1},1),size(BrCell{1},1));
    d = sparse(size(CpCell{1},2),1);
    e = sparse(size(BrCell{1},1),1);
    bp = sparse(size(CpCell{1},2),1);
    for i=1:numel(mp_strct.patch_arr)
        F = F + CpCell{i}'*ArpCell{i}'*(decArrCell{i}\(ArpCell{i}*CpCell{i})) - CpCell{i}'*AppCell{i}*CpCell{i};
        G = G + BrCell{i}*(decArrCell{i}\(ArpCell{i}*CpCell{i}));
        W = W + BrCell{i}*(decArrCell{i}\(BrCell{i}'));
        d = d + CpCell{i}'*ArpCell{i}'*(decArrCell{i}\frCell{i}');
        e = e + BrCell{i}*(decArrCell{i}\frCell{i}');
        bp = bp + CpCell{i}'*fpCell{i}';
    end
    d = d - bp;

    decF = decomposition(-F,"chol","upper");
    t2 = toc;
    fprintf('\t\tAssembling of coarse problem and factorization of F took %fs\n',t2);

    Sfun = @(x) G*(decF\(G'*x))+W*x;


    %% Solver after second Schur complement
%     fprintf('\t\tStart Solving\n');
    tic;
    [lambda,~,~,iter,~,eigEst] = pcg_w_eigest(Sfun,(G*(decF\d)+e),tol,1e3,@(x) dirPrec(BrCell,ArrCell,remIntfCell,remVolCell,x));
    t3 = toc();
    fprintf('\t\tSolving interface probl. took %fs for %i iterations\n',t3,iter);

    condEst = eigEst(2)/eigEst(1);
    
    % Compute p
    p = decF\(d-G'*lambda);
    % Compute ap
    apCell = cellfun(@(Cp) (-Cp*p)',...
        CpCell, 'UniformOutput', false);
    % Compute ar
    arCell = cellfun(@(decArr,fr,Arp,ap,Br) (decArr\(fr'-Arp*ap'-Br'*lambda))',...
        decArrCell, frCell, ArpCell, apCell, BrCell, 'UniformOutput', false);

  

    %% Change back the ordering
    rem = []; pri = []; dir = [];
    ad = []; 
    for i=1:numel(mp_strct.patch_arr)
        rem = [rem;remCell{i}(:)+mp_strct.cumu_dofs(i)];
        pri = [pri;priCell{i}(:)+mp_strct.cumu_dofs(i)];
        dir = [dir;dirCell{i}(:)+mp_strct.cumu_dofs(i)];
        ad = [ad;adCell{i}];
    end
    sol = zeros(mp_strct.cumu_dofs(end),1);
    sol(rem) = [arCell{:}]';
    sol(pri) = [apCell{:}]';
    sol(dir) = ad;
end

function res = lumpedPrec(BrCell,ArrCell,remIntfCell,x)
    res = zeros(size(x));
    for i=1:numel(BrCell)
        intf = remIntfCell{i};
        % TODO: Add material scaling if discontinuous
        res = res + BrCell{i}(:,intf)*(ArrCell{i}(intf,intf)*(BrCell{i}(:,intf)'*x));
    end
end

function res = dirPrec(BrCell,ArrCell,remIntfCell,remVolCell,x)
    res = zeros(size(x));
    for i=1:numel(BrCell)
        intf = remIntfCell{i};
        vol = remVolCell{i};
        % TODO: Add material scaling if discontinuous
        vec = (BrCell{i}(:,intf)'*x);
        res = res + BrCell{i}(:,intf)*...
            (ArrCell{i}(intf,intf)*vec ...
            - ArrCell{i}(intf,vol)*(ArrCell{i}(vol,vol)\(ArrCell{i}(vol,intf)*vec)));
    end
end