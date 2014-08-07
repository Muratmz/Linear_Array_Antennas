function MVMOSH_demo
% ================ MVMO-SH Search Progress Visualization ================ %
% Application to N-Dimensional Continuous Minimization Problems           %
% By:                                                                     %
% Dr.-Ing. Jose L. Rueda (jose.rueda@uni-duisburg-essen.de)               %
% Prof. Dr.-Ing. István Erlich (istvan.erlich@uni-due.de)                 %
% Sebastian Wildenhues, M.Sc. (sebastian.wildenhues@uni-due.de)           %
%                                                                         %
%                                                                         %
% Date: 26/06/2013                                                        %                                        
% ======================================================================= %                        

%% Hints:
% (1)
% While testing and checking out, use the MATLAB-inherent 'Zoom'- and 
% 'Rotate'-features that can be found in the corresponding toolbars!
% The 'Freeze'-button can, thereby, help you to observe the charcteristics 
% of MVMO-SH-Search-Progress more carefully.

% (2)
% For better readability and clearness, collapse all of the code by placing
% your cursor anywhere within the file, then right-click, and select
% 'Code Folding > Collapse All (or) Fold All' from the context menu.
% Thereafter, you can stepwise expand the folded code again.

    MVMO_SH_run
    
    %% Job control (for nomenclature see function 'MVMOSH')
    function MVMO_SH_run
        close all
        close all hidden
        clear all
        clear variables
        clear classes
        clc 
    
        global sn              
        global iopt
        global nnsave
        global hfig_i
        global act_con_val
        global testcase
        global parameter
        global max_run
        global printff fe_count1
        global ratio_gute_max ratio_gute_min local_search0_percentage
            
        %% Default settings
        testcase=1; %Test case
        printff=25; %Frequency of output result update
        parameter.MaxEval=1000; %Max number of function evaluations
        parameter.n_var=2; %Problem dimension
        
        % Regular MVMO-SH control
        parameter.n_par=1; %Number of particles     
        ratio_gute_max=0.7; %Max percentage of good particles
        ratio_gute_min=0.2; %Min percentage of good particlest
        parameter.n_tosave=4; %Size of the solution archive
        parameter.n_random_ini=2; %Initial number of variables selected for mutation
        parameter.n_random_last=1; %Final number of variables selected for mutation                      
        parameter.Indep_run=2; %Independent evaluation of every particle
        
        % Advanced MVMO-SH control
        %Selection strategy for offspring creation [1,2,3,4,5]: 
        %1-Random, 2-Neighbor group block stepping, 3-Neighbor group single stepping
        %4-Sequential random selection, 5-Roulette tournament
        parameter.mode=4; 
        parameter.fs_factor_start=1; %Initial shape scaling factor
        parameter.fs_factor_end=1; %Final shape scaling factor
        parameter.dddd = 1d0; %Initial value of alternative shape
        parameter.delta_dddd_start=0.05; %Initial scale factor for alternative shape 
        parameter.delta_dddd_end=0.05; %Final scale factor for alternative shape         
        parameter.nrand_method=41; %Decrement strategy of variables selected for mutation: [1,2,3]:
                                    %1-Linear, 2-Quadratic progressive, 3-Quadratic degressive
                                    %plus 40 means variable randomly changed within the range
        local_search0_percentage = 0; % Probability of local search (percentage / number of optimization variables)
        rand('state',sum(100*clock)); %Seed of the random number generator
        h_function_show([],[],[],[],[],[],[],[],[],[])
        sn = true; %Print output results
        nnsave=parameter.MaxEval;
        max_run=1;
        act_con_val=1;                              
        hfig_i=1;

        parameter.scaling=parameter.x_max - parameter.x_min; %Control range of optimization variables for normalization 
        for iopt = 1:max_run %Number of stochastically independent opimization runs
            [FConv(:,iopt),xbest]=MVMOSH(iopt);
        end
    end

    %% Core of MVMOSH (slightly modified to facilitate visualization)
    function [ofcn,best] = MVMOSH(krun)
    %By:  
    % Dr. Jose L. Rueda (jose.rueda@uni-duisburg-essen.de)
    % Prof. Istvan Erlich (istvan.erlich@uni-due.de)
    % Hybrid variant of the Mean-Variance Mapping Optimization (MVMO-SH) - Version 2013 ©
    %              -------------------*-------------------

    %              -------------------*-------------------
    %============================== Nomenclature ==============================
    % max_eval :=  Maximum number of objective function evaluations
    % n_var := Number of parameters to be optimized, i.e. the dimensionality of the problem.
    % n_par := Number of particles
    % n_to_save := Size of the solution archive
    % n_randomly_ini := Initial number of variables selected for mutation 
    % n_randomly_last := Minimum number of variables selected for mutation
    % Indep_run := Independent evaluation of every particle
    % mode := Selection strategy for offspring creation
    % fs_factor_start := Initial shape scaling factor
    % fs_factor_end := Final shape scaling factor
    % dddd := Initial value of alternative shape 
    % delta_dddd_start := Initial scale factor for alternative shape 
    % delta_dddd_end := Final scale factor for alternative shape
    % nrand_method := Decrement strategy of variables selected for mutation

    %================================ Syntax ==================================
    %  [ofcn, best] = MVMOS_SH(fitfcn,repcont,varargin) 
    %  Inputs:
    %  fitfcn   - Function handle for calculation of objective function and constraints
    %  repcont  - Actual repetition of the optimization
    %  varargin - Problem number
    %  Outputs:
    %  ofcn - Recorded objective function values in steps as defined by parameter.rec_step
    %  best - Values of optimization variables in the last iteration
    %==========================================================================

    global parameter;
    global x_normalized_best  
    global considered changed fs_factor
    global printff sn best_history
    global xi_s ipp_gbest ipp 
    global shape
    global meann izm izz l_vari
    global testcase
    global x_normalized_pic
    global table 
    global variance best_results
    global na_pg fe_count1 
    global ratio_gute_max ratio_gute_min local_search0_percentage
    
    xi_s = 0:0.01:1; %Setpoints for visualization of mapping functions

    %%----------------- Create initial random population ------------------   
    xx=zeros(parameter.n_par,parameter.n_var);
    x_norm=xx;
    for iijj=1:parameter.n_par %Start initialization: x_normalized
        for jjkk=1:parameter.n_var
            xx(iijj,jjkk)=parameter.x_min (jjkk) + rand*(parameter.x_max(jjkk)-parameter.x_min(jjkk));
        end
        x_norm(iijj,:)=(xx(iijj,:)-parameter.x_min)./parameter.scaling;
    end % End initialization
    x_normalized=x_norm;

    
    %% ------------------ Initialize control parameters --------------------
    max_eval = parameter.MaxEval;
    n_par = parameter.n_par; 
    n_var = parameter.n_var; 
    n_to_save = parameter.n_tosave; 
    n_randomly_ini = parameter.n_random_ini;
    n_randomly_last = parameter.n_random_last;
    Indep_run=parameter.Indep_run;
    mode = parameter.mode;
    fs_factor_start = parameter.fs_factor_start;
    fs_factor_end = parameter.fs_factor_end;
    dddd = ones(n_par,n_var)*parameter.dddd;
    delta_dddd_start = parameter.delta_dddd_start; 
    delta_dddd_end = parameter.delta_dddd_end; 
    nrand_method=parameter.nrand_method;
    
    
    %% --------------------- Data structure for the table ---------------------
    table.bests = zeros(parameter.n_tosave,parameter.n_var,n_par);
    table.fitness = Inf*ones(parameter.n_tosave,1,n_par);
    table.objective = Inf*ones(parameter.n_tosave,1,n_par);
    best_history=NaN*zeros(parameter.MaxEval,parameter.n_var);

    %% ----------------------------- Mapping ----------------------------------
    shape = zeros(n_par,n_var);
    x_normalized_best = x_normalized;
    meann = x_normalized;
    meann_app = x_normalized;
    gbest_sss = zeros(2,n_var);

    %% ------------------------ Variable selection ----------------------------
    izm = zeros(1,n_par);
    izz = zeros(1,n_par);
    considered = true(n_par,n_var);
    variance =ones(n_par,n_var);
    probab=ones(n_par,n_var);
    values_noneq = zeros(n_to_save,1);
    
    if (n_randomly_last<1)
        n_randomly_last=1;
    end
    
    if (n_randomly_last>n_var)
        n_randomly_last=n_var;
    end
    
    if (n_randomly_ini>n_var)
        n_randomly_ini=n_var;
    end  
    
    if (mode<1)
        mode=4d0;
    end
    
    if (mode>5)
        mode=4d0;
    end
    
    if (max_eval<=0)
        max_eval=100000d0;
    end

    if (fs_factor_start<=0)
        fs_factor_start=1d0;
    end
    
    if (fs_factor_end<=0)
        fs_factor_end=1d0;
    end

    fs_factor0=fs_factor_start;
    
    if (n_to_save<=1)
        n_to_save=2d0;
    end
    
    if (delta_dddd_start<=0)
        delta_dddd_start=0.2d0;
    end
    
    if (delta_dddd_end<=0)
        delta_dddd_end=0.2d0;
    end
    
    delta_dddd0=delta_dddd_start;
    delta_dddd1=1.d0+delta_dddd0;
    
    yes_fs_factor=true;
    if (fs_factor_start==fs_factor_end)
        yes_fs_factor=false;
    end
    
    yes_delta_dddd=true;
    if (delta_dddd_start==delta_dddd_end)
        yes_delta_dddd=false;
    end
  
    yes_n_randomly=true;
    if (n_randomly_ini==n_randomly_last)
        n_randomly=n_randomly_last;
        yes_n_randomly=false;
    end

    %% ----------------------------- Counters ---------------------------------
    fe_count1w=0;
    fe_count1=0;
    na_pg=n_par;

    no_in = zeros(1,n_par);
    no_inin = zeros(1,n_par);
    n_par_last=n_par;
    
    indpendent_run_p=zeros(1,n_par);
    local_search0=local_search0_percentage/100/n_var; % Probability of local search (percentage / number of optimization variables)
    best_results=NaN*zeros(max_eval,1);
    goodbad=zeros(n_par,1);
    make_sense=max_eval-round(n_var*100.d0/2);

    counter_fmincon=0;
    counter_fmincon_FEVALS=0;

    l_vari=2; 
    if (l_vari<2)
        l_vari=2;
    end

    %% --------------------------------- MVMO-SH ---------------------------------
    if (n_var==2)
    fprintf('===============================================================================\n');
    fprintf('                   Mean-Variance Mapping Optimization Results                  \n');
    fprintf('===============================================================================\n');
    fprintf('Nomenclature:                                                                  \n');
    fprintf('FE-No. => Objective function evaluation number                                 \n');
    fprintf('NP     => Number of active particles                                           \n');
    fprintf('AGBP   => Actual global best particle                                          \n');
    fprintf('OF     => Global best objective function                                       \n');
    fprintf('x1,x2  => Optimization variables                                               \n');
    fprintf('No_fmincon  => Calls of local search                                           \n');
    fprintf('Ev_fmincon  => Number of function evaluations performed by local search        \n');
    fprintf('\n');
    fprintf('  FE-No.      NAP        AGBP            OF                    x1                  x2          No_fmincon     Ev_fmincon \n');
    fprintf('  ------     ------     ------        ---------             --------            ---------     ------------    ---------- \n');
    else
    fprintf('===============================================================================\n');
    fprintf('                   Mean-Variance Mapping Optimization Results                  \n');
    fprintf('===============================================================================\n');
    fprintf('Nomenclature:                                                                  \n');
    fprintf('FE-No. => Objective function evaluation number                                 \n');
    fprintf('NP     => Number of active particles                                           \n');
    fprintf('AGBP   => Actual global best particle                                          \n');
    fprintf('OF     => Global best objective function                                       \n');
    fprintf('No_fmincon  => Calls of local search                                           \n');
    fprintf('Ev_fmincon  => Number of function evaluations performed by local search        \n');
    fprintf('\n');
    fprintf('  FE-No.      NAP        AGBP            OF          No_fmincon       Ev_fmincon  \n');
    fprintf('  ------     ------     ------        ---------     ------------      ----------  \n');
    end
    
   delta_nrandomly=n_randomly_ini-n_randomly_last;
   tol_flag=0;
   A=zeros(n_par,1); 
   
   while (tol_flag==0)       
        ff=real(fe_count1/max_eval);
        ff2=ff*ff;
        one_minus_ff2=1.d0-ff2*0.99;
        shift=one_minus_ff2*0.5;
        
        %Determining the relative number of particles belonging to the group of
        %good particles
        border_gute0=ratio_gute_max-ff*(ratio_gute_max-ratio_gute_min);
        border_gute=round(n_par_last*border_gute0);
         
        %Selecting the subset of variables to be mutated
        if yes_n_randomly
            switch nrand_method
                case 1
                    n_randomly=round(n_randomly_ini-ff*delta_nrandomly);
                case 2
                    n_randomly=round(n_randomly_ini-ff2*delta_nrandomly); 
                case 3
                    n_randomly=round(n_randomly_last+(1.0-ff)*(1.0-ff)*delta_nrandomly); 
                case 4
                    n_randomly=round(unifrnd(n_randomly_last,n_randomly_ini)); 
                otherwise
                    if nrand_method==42
                        n_randomly_X=round(n_randomly_ini-ff2*delta_nrandomly); 
                    elseif nrand_method==43
                        n_randomly_X=round(n_randomly_last+(1.0-ff)*(1.0-ff)*delta_nrandomly); 
                    else
                        n_randomly_X=round(n_randomly_ini-ff*delta_nrandomly);
                    end
                    n_randomly=round(n_randomly_last+rand*(n_randomly_X-n_randomly_last));
            end
        end

        %Quadratic variation of fs_factor0
        if yes_fs_factor
            fs_factor0=fs_factor_start+ff2*(fs_factor_end-fs_factor_start); 
        end
        
        %Quadratic variation of delta_dddd0
        if yes_delta_dddd
            delta_dddd0=delta_dddd_start+ff2*(delta_dddd_end-delta_dddd_start); 
            delta_dddd1=1.d0+delta_dddd0;
        end  
        
        %Evaluating the particles.....
        for ipp=1:parameter.n_par  
            if rand < local_search0 &&  fe_count1 < make_sense &&  fe_count1 > 1    
                [msgstr, msgid] = lastwarn;
                TFrcond = strcmp('MATLAB:nearlySingularMatrix',msgid); % Only informative from 'fmincon' function 
                if TFrcond~=0
                    rcond_value0=str2num(msgstr(81:end-1));
                end
                
                [xt.fitness,x_normalized(ipp,:),FEVALS] = LocalSearchMVMOSH(x_normalized(ipp,:),testcase); %Local search
                
                [msgstr1, msgid1] = lastwarn;
                TFrcond1 = strcmp('MATLAB:nearlySingularMatrix',msgid1);
                if TFrcond1~=0
                    rcond_value=str2num(msgstr1(81:end-1));
                    if (exist('rcond_value0','var')==1) && (rcond_value0 ~= rcond_value) 
                        local_search0=0;
                    end
                end
                
                ffx1=xt.fitness;
                x_normalized_pic = parameter.x_min + parameter.scaling .* x_normalized(ipp,:); 
                counter_fmincon=counter_fmincon+1; %Only for printing
                counter_fmincon_FEVALS=counter_fmincon_FEVALS+FEVALS; %Only for printing
                
                for iq=fe_count1+1:fe_count1+FEVALS-1
                   if iq <= max_eval
                       best_results(iq)=best_results(iq-1);
                       best_history(iq,:)=best_history(iq-1,:);
                   else
                       break
                   end
                end
                
                fe_count1=fe_count1+FEVALS;
                if fe_count1 <= max_eval
                   if fe_count1 >= 2
                       if xt.fitness < best_results(fe_count1-1)
                           best_results(fe_count1)=xt.fitness;
                           ipp_gbest=ipp;
                       else
                              best_results(fe_count1)=best_results(fe_count1-1);
                       end
                   else
                          best_results(fe_count1)=xt.fitness;
                          ipp_gbest=ipp;
                   end
                end
            else
                [xt.fitness,x_normalized(ipp,:)] = fitness_evaluation(x_normalized(ipp,:),testcase); %Evaluation of fitness
                fe_count1=fe_count1+1;
                ffx1=xt.fitness;
                x_normalized_pic = parameter.x_min + parameter.scaling .* x_normalized(ipp,:); 
                if fe_count1 <= max_eval
                    if fe_count1 >= 2
                        if xt.fitness < best_results(fe_count1-1)
                            best_results(fe_count1)=xt.fitness;
                            ipp_gbest=ipp;
                        else
                            best_results(fe_count1)=best_results(fe_count1-1);
                        end
                    else
                        best_results(fe_count1)=xt.fitness;
                        ipp_gbest=ipp;
                    end
                end
            end
            best_history(fe_count1,:)=parameter.x_min+(parameter.x_max-parameter.x_min).*table.bests(1,:,ipp_gbest);
                
            % Store the n-best solution corresponding to the corresponding particle's archive
            Fill_solution_archive();
            meann_app(ipp,:)= meann(ipp,:);
            
            x_normalized4commandwindow = parameter.x_min + parameter.scaling .* x_normalized_best(ipp_gbest,:);
 
            if (n_var>2)
                if (fe_count1 == 1) 
                    fprintf('%8d    %5d     %5d      %17.7E     %5d        %5d\n',...
                        fe_count1,na_pg,ipp_gbest,best_results(fe_count1),counter_fmincon,counter_fmincon_FEVALS); 
                elseif (mod(fe_count1,printff) == 0)
                    fec_prev=fe_count1;
                    fprintf('%8d    %5d     %5d      %17.7E     %5d        %5d\n',...
                        fe_count1,na_pg,ipp_gbest,best_results(fe_count1),counter_fmincon,counter_fmincon_FEVALS); 
                elseif (exist('fec_prev','var') && (fe_count1>(fec_prev+printff)))
                    fec_prev=fec_prev+printff;
                    fprintf('%8d    %5d     %5d      %17.7E     %5d        %5d\n',...
                        fec_prev,na_pg,ipp_gbest,best_results(fec_prev),counter_fmincon,counter_fmincon_FEVALS); 
                end
            elseif ((fe_count1 == 1) || (mod(fe_count1,printff) == 0)) && sn
                if (fe_count1 == 1) 
                    fprintf('%8d    %5d     %5d      %17.7E    %17.7E    %17.7E     %5d        %5d\n',...
                        fe_count1,na_pg,ipp_gbest,best_results(fe_count1),x_normalized4commandwindow(1),x_normalized4commandwindow(2),counter_fmincon,counter_fmincon_FEVALS); 
                elseif (mod(fe_count1,printff) == 0)
                    fec_prev=fe_count1;
                    fprintf('%8d    %5d     %5d      %17.7E    %17.7E    %17.7E     %5d        %5d\n',...
                        fe_count1,na_pg,ipp_gbest,best_results(fe_count1),x_normalized4commandwindow(1),x_normalized4commandwindow(2),counter_fmincon,counter_fmincon_FEVALS); 
                elseif (exist('fec_prev','var') && (fe_count1>(fec_prev+printff)))
                    fec_prev=fec_prev+printff;
                    fprintf('%8d    %5d     %5d      %17.7E    %17.7E    %17.7E     %5d        %5d\n',...
                        fec_prev,na_pg,ipp_gbest,best_results(fec_prev),x_normalized4commandwindow(1),x_normalized4commandwindow(2),counter_fmincon,counter_fmincon_FEVALS); 
                end
            end
              
            if (na_pg==1) %Only one particle
                x_normalized(ipp,:)=x_normalized_best(ipp,:);
                goodbad(1,1)=1;
            
            else %Several particles
                indpendent_run_p(ipp)=indpendent_run_p(ipp)+1;
                if (indpendent_run_p(ipp)>=Indep_run)  %Non-independent evaluation of particles  
                        
                    if fe_count1w<1 
                        fe_count1w=fe_count1;
                    end
                    
                    %Determining the proportion of good particles
                    A(1:n_par)=table.fitness(1,:,1:n_par);
                    [~,IX] = sort(A);
                                        
                    goodbad(IX(border_gute+1:n_par),1)=0;
                    goodbad(IX(1:border_gute),1)=1;
                    iec=randi(border_gute-2,1,1);                                      

                    bestp=IX(1); %First (global best) particle of the group of good particles
                    onep=IX(iec+1); %Randomly selected intermediate particle of the group of good particles
                    worstp=IX(border_gute); %Last particle of the group of good particles
                    
                    %Multi-parent strategy for bad particles
                    if ~goodbad(ipp)   
                        beta1 = 2.0*(rand - shift);
                        for jl=1:n_var
                            x_normalized(ipp,jl) =x_normalized_best(onep,jl)+beta1*(x_normalized_best(bestp,jl)-x_normalized_best(worstp,jl));
                            while x_normalized(ipp,jl) > 1.0d0 ||  x_normalized(ipp,jl)<0.0d0  
                                beta2 = 2.0*(rand - shift);  
                                x_normalized(ipp,jl) =x_normalized_best(onep,jl)+beta2*(x_normalized_best(bestp,jl)-x_normalized_best(worstp,jl));
                            end
                        end
                        meann_app(ipp,1:n_var) =x_normalized(ipp,1:n_var); %Update the mean values for all variables of the ipp particle
                    else
                        x_normalized(ipp,:)= x_normalized_best(ipp,:); %Local best-based parent assignment for good particles
                    end
                    else
                    x_normalized(ipp,:)=x_normalized_best(ipp,:); %Local best-based parent assignment during independent evaluations
                end
                
            end
            
            %Change the corresponding variable(s)
            if mode==5 %Roulette tournament
                probab(ipp,:) =1.d0./(dddd(ipp,:));
            end
            considered(ipp,1:n_var) = false;
                 
            VariableSelect1(); % Call random variable selection strategy
            
            %Generation of random input for the mapping function
            for jl=1:n_var
                if considered(ipp,jl) 
                    x_normalized(ipp,jl) = rand();
                end
            end
            
            %Applying the mapping function for seleced variables
            for ivar=1:n_var
                if  considered(ipp,ivar)
                    sss1 = shape(ipp,ivar);
                    sss2 = sss1; 
                    delta_ddd_x=delta_dddd0*(rand()-0.5d0)*2.0d0+delta_dddd1; 
                    if (shape(ipp,ivar) > dddd(ipp,ivar))
                        dddd(ipp,ivar) = dddd(ipp,ivar)*delta_ddd_x;
                    else
                        dddd(ipp,ivar) = dddd(ipp,ivar)/delta_ddd_x;
                    end
                    
                    isteur = randi(2,1)-1;
                    if isteur %Random assignment of shape
                        sss1=dddd(ipp,ivar) ;
                    else
                        sss2=dddd(ipp,ivar);  
                    end
                    
                    fs_factor=fs_factor0*(1d0+rand);
                    sss1=sss1*fs_factor;
                    sss2=sss2*fs_factor;
                    
                    x_normalized_temp = x_normalized(ipp_gbest,:);
                    
                    x_normalized(ipp,ivar)=...
                        h_function(meann_app(ipp,ivar),sss1,sss2,x_normalized(ipp,ivar)); %Mapping function
                    
                    if (ipp_gbest==ipp)
                        gbest_sss(1,ivar)=sss1;
                        gbest_sss(2,ivar)=sss2;
                    end
                    
                        h_function_show(fe_count1,ivar,meann_app(ipp_gbest,:),shape(ipp_gbest,:),gbest_sss(1,ivar),gbest_sss(2,ivar),x_normalized_temp,x_normalized(ipp_gbest,:),x_normalized4commandwindow,ffx1)           
                end
            end                                     
           
            if (fe_count1>=max_eval)
                tol_flag=1;
                break;
            end    
        end %End n_par loop
    end %End while loop
       
    ofcn=best_results;
    best=parameter.x_min+(parameter.x_max-parameter.x_min).*table.bests(1,:,ipp_gbest);
    best_history(fe_count1,:) = best;
    
    %% ----------------------- Complementary functions ------------------------
    function [FUN,xn_out,FEVALS] = LocalSearchMVMOSH(xx_yy1,testcase)   
        lb=parameter.x_min;
        ub=parameter.x_max;
        Aeq=[]; %AeqX=Beq
        Beq=[];
        AA=[]; %AX<=B
        BB=[];
        
        %Denormalize to the original [min, max] range 
        xx_yy = parameter.x_min+parameter.scaling.* xx_yy1;
        
        options=optimset('Display','off','algorithm','interior-point','UseParallel','never'); 
        [Xsqp, FUN , ~ , output]=...
            fmincon(@(xx_yy)LSearch(xx_yy,testcase),xx_yy,AA,BB,Aeq,Beq,lb,ub,[],options);
        
        FEVALS=output.funcCount;
        for nvar=1:size(xx_yy,2)
            if isnan(Xsqp(1,nvar))
                Xsqp=xx_yy;
                break;
            end
        end
        
        % Normalize again to [0, 1]
        xn_out = (Xsqp-parameter.x_min)./parameter.scaling;
    end

    function J=LSearch(xx_yy,testcase)
        xx_yy2 = (xx_yy-parameter.x_min)./parameter.scaling;
        [J,~] = fitness_evaluation(xx_yy2,testcase);
    end

        function Fill_solution_archive()
        no_in(ipp) = no_in(ipp)+1;
        changed = false;

        if no_in(ipp) ==1 % the solution coming to the table for the first time
            table.bests(1:n_to_save,:,ipp) = repmat(x_normalized(ipp,:),n_to_save,1);
            table.fitness(1:n_to_save,:,ipp) = repmat(xt.fitness,n_to_save,1);
            table.objective(1:n_to_save,:,ipp) = repmat(xt.fitness,n_to_save,1);

            x_normalized_best(ipp,:) = table.bests(1,:,ipp); % set the best solution to the one of the first rank
            no_inin(ipp)=no_inin(ipp)+1;
            
        else % not for the first time and check for the update
           i_position =0;
           % check if the new coming solution is less than any in the table
           if (xt.fitness<table.fitness(n_to_save,:,ipp))
               for i=1:n_to_save 
                   if (xt.fitness<table.fitness(i,:,ipp))                                  
                       i_position = i;
                       changed =true;
                       if (i<n_to_save)
                           no_inin(ipp) = no_inin(ipp)+1; % how many times good solutions were found   
                       end
                       break;
                   end
               end
           end
        end

        if changed   % if the new individual is better than any archived individual.
                     % Move the individuals and corresponding fitness values
                     % downward so that the individuals are sorted based on the
                     % fitness value in a descending order             
            nnnnnn = n_to_save;
            if (no_inin(ipp) < n_to_save); nnnnnn = no_inin(ipp); end
            isdx=nnnnnn:-1:i_position+1;
            table.bests(isdx,:,ipp) = table.bests(isdx-1,:,ipp);
            table.fitness(isdx,:,ipp)= table.fitness(isdx-1,:,ipp);
            table.objective(isdx,:,ipp)= table.objective(isdx-1,:,ipp);
 
            % save the new best
            table.bests(i_position,:,ipp) = x_normalized(ipp,:);
            table.fitness(i_position,:,ipp) = xt.fitness;
            table.objective(i_position,:,ipp) = xt.fitness;

            % calculation of mean and variance
            if ((no_inin(ipp)>=l_vari))
                for ivvar = 1:n_var
                    [meann(ipp,ivvar),variance(ipp,ivvar)] = mv_noneq(nnnnnn,table.bests(1:nnnnnn,ivvar,ipp));
                 end
                id_nonzero = (variance(ipp,:) > 1.1d-100);
                shape(ipp,id_nonzero)=-log(variance(ipp,id_nonzero));
            end
        end
            x_normalized_best(ipp,:) = table.bests(1,:,ipp);   
        end
        
        function VariableSelect1()   
          switch mode
            case 1
                for ii=1:n_randomly
                    isrepeat = false;
                    while ~isrepeat
                        inn=round(rand*(n_var-1))+1;
                        if (~considered(ipp,inn))
                            isrepeat = true;
                        end
                    end
                    considered(ipp,inn)=true;
                end
            case 2
                in_randomly=0;
                isrepeat = false;
                izz(ipp)=round(rand*(n_var-1))+1; 
                while ~isrepeat
                    in_randomly=in_randomly+1;
                    if (izz(ipp)< 1) 
                        izz(ipp)=n_var;
                    end
                    considered(ipp,izz(ipp))=true;
                    izz(ipp)=izz(ipp)-1;
                    if (~(in_randomly<n_randomly)) 
                        isrepeat = true;
                    end
                end
            case 3
                in_randomly=0;
                izm(ipp)=izm(ipp)-1;
                if (izm(ipp)< 1) 
                    izm(ipp)=n_var;
                end
                izz(ipp)=izm(ipp);
                isrepeat = false;
                while ~isrepeat
                    in_randomly=in_randomly+1;
                    if (izz(ipp)< 1) 
                         izz(ipp)=n_var;
                    end
                    considered(ipp,izz(ipp))=true;
                    izz(ipp)=izz(ipp)-1;
                    if (~(in_randomly<n_randomly)) 
                        isrepeat = true;
                    end
                end   
            case 4
                izm(ipp)=izm(ipp)-1;
                if (izm(ipp)< 1) 
                    izm(ipp)=n_var;
                end
                considered(ipp,izm(ipp))=true;
                if (n_randomly>1)  
                    for ii=1:n_randomly-1
                        isrepeat = false;
                        while ~isrepeat
                            inn=round(rand*(n_var-1))+1;
                            if (~considered(ipp,inn))
                                isrepeat = true;
                            end
                        end
                        considered(ipp,inn)=true;
                    end
                end
            case 5
                  summep=sum(probab(ipp,:));
                  wahr=probab(ipp,:)/summep;
                  SS0=0.d0;
                  SS=zeros(1,n_var);
                  for imm=1:(n_var-1)
                      SS0=SS0+wahr(imm);
                      SS(imm+1)=SS0;
                  end
                  for ijlr=2:n_var
                      wahr(ijlr)=wahr(ijlr)+SS(ijlr);
                  end 
                  for ltt=1:n_randomly
                       isrepeat = false;
                       while  ~isrepeat
                        rnn=rand;
                        for irw=1:n_var
                          if considered(ipp,irw)
                           continue
                          end
                          if (irw==1)
                              unten=0.d0;
                          else
                              unten=wahr(irw-1);
                          end
                         if (rnn>=unten)&&(rnn<wahr(irw))
                              isrepeat = true;
                             considered(ipp,irw)=true;
                             break;
                          end
                        end
                       end
                  end
          end
        end

        function [vmean,vvariance] = mv_noneq(n_to_save,values)
            iz =1; 
            values_noneq(iz)=values(1);
            for ii_jj=2:n_to_save
                izz = iz;
                gleich = false;
                for kk_ii=1:izz
                    if  abs(values_noneq(kk_ii) - values(ii_jj)) <  1.d-70 ; 
                        gleich = true;
                        break;
                    end
                end
                if ~gleich;
                    iz = iz+1;
                    values_noneq(iz)=values(ii_jj);
                end
            end
             vmean = values_noneq(1);
            if (iz>1)
               for kk_ii=2:iz
                    vmean = vmean+values_noneq(kk_ii);
               end
               
               vmean = vmean/iz;
           end
           vvariance = 0.d0;
           if (iz>1)
                for kk_ii=1:iz
                    vvariance =vvariance+(values_noneq(kk_ii)-vmean)*(values_noneq(kk_ii)-vmean);
                end
                  vvariance = vvariance/iz;
           else
                  vvariance=1.0d-100 ;
            end
        end 
    end

    %% Evacuated h-function
    function x = h_function(x_bar,s1,s2,x_p)
        H = x_bar .* (1.d0 - exp(-x_p.*s1)) + ...
            (1.0d0 - x_bar) .* exp(-(1.d0-x_p).*s2);              
        H0 = (1.d0 - x_bar) .* exp(-s2);
        H1 = x_bar .* exp(-s1);
        x = H + H1 .* x_p + H0 .* (x_p - 1.d0);
    end

 
    %% Core of visualization
    function h_function_show(fe_count1,ivar,meann,shape,s1,s2,x_random_normalized,xi_real,x_normalized4commandwindow,fitness2show)
        global hAxes_surfplot
        global parameter x_normalized_pic best_results
        global printff sn ipp_gbest
        global testcase best_history
        global X Y Z
        global xi_s na_pg
        global hfig_i 
        global myColormap    
        global considered
        global act_con_val
        global table
        global iterf x12plot x22plot
        global Show_Fit_Value Show_Archive_Value
        global variance 
        persistent hAxes
        persistent fnccall_cnt                                              
        persistent s14fig_vec s24fig_vec s14fig_temp s24fig_temp
        persistent title1_handle Cancel_button
        persistent ph
        persistent considered_count
        persistent considered_count_handle na_pg_handle fitness_h_axes_handle
        persistent considered_count_overall 
        persistent cb
        persistent x_random_normalized_save
        persistent xi_real_save
        persistent cf 
        persistent Best_fitness_handle Best_particle_no_handle main_axes
        persistent Freeze
        persistent archivedata
        persistent row_names col_names look4archivewindow
            if (isempty(fe_count1) && isempty(ivar) && isempty(meann) && isempty(shape) && isempty(s1) && ...
                    isempty(s2) && isempty(x_random_normalized) && isempty(xi_real) && ...
                    isempty(x_normalized4commandwindow) && isempty(fitness2show))
                        Show_Fit_Value = 1;
                        Show_Archive_Value = 1;
                        Freeze = uicontrol ('style','togglebutton', ...
                        'units','normalized', ...
                        'position', [0.02 0.0325 0.07 0.045], ...
                        'callback',@Freeze_Figure, ...
                        'string','Freeze Figure','Visible','off');

                        Quit_button = uicontrol('style','pushbutton', ...
                        'units','normalized', ...
                        'position', [0.02 0.0325 0.07 0.045], ...
                        'callback','quit', ...
                        'string','Close Demo');

                        hhhh = uicontrol ('style','pushbutton', ...
                        'units','normalized', ...
                        'position', [0.47 0.0325 0.07 0.045], ...
                        'callback','uiresume(gcbf)', ...
                        'string','Start Evaluation','FontWeight','bold');
                        
                        myColormap = colormap('autumn');
                        myColormap = flipud(myColormap);
                        
                        cc = uicontrol('Style', 'popup',...
                       'String', 'Change the display colormap||jet|hsv|hot|cool|autumn|winter',...
                       'Units','normalized',...
                       'Position', [0.02 0.72 0.125 0.045],...
                       'Callback', @setmap);

                        pp = uicontrol('Style', 'popup',...
                       'String', 'Select Problem||De Jong’s function|Axis parallel hyper-ellipsoid function|Rotated hyper-ellipsoid function|Rosenbrock’s valley|Rastrigin’s function|Schwefel’s function|Griewangk’s function|Sum of different power functions|Ackley’s function|Michalewicz’s function|Branin’s function|Easom’s function|Goldstein-Price’s function|“Drop wave” function|Shubert function|Deceptive function Type III, piecewise linear (Only available for 2D)|Deceptive function Type III, piecewise non-linear (Only available for 2D)',...
                       'Units','normalized',...
                       'Position', [0.02 0.75 0.125 0.045],...
                       'Callback', @Prob_Selection);

                        boundaries(testcase);

                        whitebg([0 .29 .57])
                        surfplot_x = 0.06;
                        surfplot_y = 0.1;
                        surfplot_xsize = 0.4;
                        surfplot_ysize = 0.7;

                        hAxes_surfplot = axes('Units','normalized','Position',...
                            [surfplot_x surfplot_y-0.005 surfplot_xsize surfplot_ysize],...
                            'XLim',[parameter.x_min(1) parameter.x_max(1)],...
                            'YLim',[parameter.x_min(2) parameter.x_max(2)],...
                            'DrawMode','fast','DataAspectRatioMode','manual',...
                            'PlotBoxAspectRatioMode','manual');

                        x_surfplot = linspace(parameter.x_min(1),parameter.x_max(1),200);
                        y_surfplot = x_surfplot;
                        [X Y] = meshgrid(x_surfplot,y_surfplot);

                        [Z,pname] = fitness_evaluation_surfplot({X Y},testcase);
           
                        sp = surf(hAxes_surfplot,X,Y,Z);
                        xlabel('X_1 \rightarrow','Units','Normalized')
                        ylabel('X_2 \rightarrow','Units','Normalized')
                        zlabel('Fitness or Objective Value \rightarrow','Units','Normalized')
                        title(['Problem: ' '"' pname],'Fontsize',12,'Fontweight','bold')
                        set(sp,'AlphaDataMapping','none')
                        set(sp,'EdgeColor','none','Facecolor','texture')

                        set(sp,'FaceLighting','phong',...
                          'FaceColor','interp',...
                          'EdgeColor','none',...
                          'BackFaceLighting','reverselit',...
                          'CDataMapping','direct')

                        light('Position',[1 3 2]);
                        light('Position',[-3 -1 3]);
                        material shiny
                        alpha(0.3)

                        cf = gcf;
                        set(gcf, 'Position', get(0,'Screensize'));
                        set(cf,'NumberTitle','off')
                        set(cf,'Name','MVMO-SH Search Progress Visualization')
                        set(cf,'ToolBar','none')

                        hhh = uicontrol ('style','pushbutton', ...
                        'units','normalized', ...
                        'position', [0.9 0.0325 0.09 0.045], ...
                        'callback',@Algorithm_Settings_Dialogue, ...
                        'string','Algorithm Settings');
                    
                        lll = uicontrol ('style','pushbutton', ...
                        'units','normalized', ...
                        'position', [0.81 0.0325 0.09 0.045], ...
                        'callback',@Procedure_Settings_Dialogue, ...
                        'string','Procedure Settings');

                        text(0.71,1.225,'MVMO-SH Search Progress Visualization',...
                                'Units','normalized','FontSize',16,'FontWeight','bold') 

                        text(0.4,1.15,' Application of n-dimensional Continuous and Discontinuous Unconstrained Minimization Problems',...
                                'Units','normalized','FontSize',12,'FontWeight','bold');
                            
                        an = annotation(gcf,'textbox',[0.4949 0.175 0.45 0.6698],'FontSize',12,...
                            'String',{'In this demo, the hybrid variant of Mean Variance Mapping Optimization algorithm (MVMO-SH) can be explored in an intuitive and accessible manner.', char(10), 'MVMO-SH is a swarm intelligence based procedure which incorporates local search and multi-parent crossover strategies to increase the search diversity while striving for a balance between exploration and exploitation. Based on recognized n-dimensional continuous and discontinuous unconstrained minimization problems, the user obtains insights into the basic ideas that underlie the method.', char(10),'In order to get a comprehensible theoretical background to understand the search process of classical MVMO-SH, it is recommended to make use of the information provided at http://www.uni-due.de/mvmo/'},...
                            'FitBoxToText','off',...
                            'EdgeColor','none',...
                            'FontUnits','normalized');

                        
                         hText = findall(gca,'type','text');
                         hUIControls = findall(gcf,'type','uicontrol');
                         set([hText;hUIControls],'units','normalized','fontunits','normalized');

                        uiwait(gcf);

                        if (parameter.n_var > 2)
                            hfig_i = 0;
                        else
                            set(hhh,'Visible','off')
                            set(lll,'Visible','off')
                            set(hhhh,'Visible','off')
                            set(cc,'Visible','off')
                            set(pp,'Visible','off')
                            set(an,'Visible','off')
                            set(Quit_button,'Visible','off')
                            set(Freeze,'Visible','on')

                            Cancel_button = uicontrol ('style','pushbutton', ...
                            'units','normalized', ...
                            'position', [0.02 0.0755 0.07 0.045], ...
                            'callback',@my_finishdlg, ...
                            'string','Cancel');

                            uicontrol ('style','pushbutton', ...
                            'units','normalized', ...
                            'position', [0.435 0.0325 0.06 0.045], ...
                            'string','Search History','Callback',@Show_Fit);

                            uicontrol ('style','pushbutton', ...
                            'units','normalized', ...
                            'position', [0.375 0.0325 0.06 0.045], ...
                            'string','Best Archive','Callback',@Show_Archive);
                         hText = findall(gca,'type','text');
                         hUIControls = findall(gcf,'type','uicontrol');
                         set([hText;hUIControls],'units','normalized','fontunits','normalized');
                        end
            else 
                iterf = 1:parameter.MaxEval;
                x12plot = NaN*zeros(parameter.MaxEval,1);
                x22plot = NaN*zeros(parameter.MaxEval,1);
                n_var = parameter.n_var;
                archivedata = NaN*zeros(parameter.n_tosave+2,n_var+1);
                n_randomly = parameter.n_random_ini; 
                if ((hfig_i == 1) & (n_var==2))
                    if isempty(fnccall_cnt)  
                        x_random_normalized_save = zeros(1,n_var);
                        x_random_normalized_in = zeros(1,n_var);
                        considered_count_overall = zeros(1,n_var);
                        xi_real_in = zeros(1,n_var);
                        s14fig_vec = zeros(1,n_var);                                  
                        s14fig_temp = zeros(1,n_randomly);                            
                        s24fig_vec = zeros(1,n_var);                                  
                        s24fig_temp = zeros(1,n_randomly);                            
                        xx = zeros(n_var,length(xi_s));                             
                        title1_handle = zeros(1,n_var);
                        considered_count_handle = cell(1,n_var);
                        na_pg_handle = cell(1,n_var);
                        fitness_h_axes_handle = cell(1,n_var);
                        ph = zeros(1,n_var);
                        considered_count = zeros(1,n_var);
                        lineh_horizontal = zeros(1,n_var);
                        lineh_vertical = zeros(1,n_var);
                        row_names = cell(1,parameter.n_tosave+2);
                        fnccall_cnt = 1;
                    end

                    if fnccall_cnt == 1
                        for hi = 1:n_var
                            if considered(ipp_gbest,hi) == true
                                considered_count_overall(hi) = considered_count_overall(hi) + 1;
                            end
                        end
                    end

                    if parameter.mode == 1 && considered(ipp_gbest,ivar) && ( mod(fe_count1,printff) == printff - 5 || mod(fe_count1,printff) == printff - 4 || mod(fe_count1,printff) == printff - 3 || mod(fe_count1,printff) == printff - 2 || mod(fe_count1,printff) == printff - 1 )
                        x_random_normalized_save(ivar) = x_random_normalized(ivar);
                        xi_real_save(ivar) = xi_real(ivar);
                        s14fig_temp(ivar) = s1;
                        s24fig_temp(ivar) = s2;
                    elseif parameter.mode ~= 1 && mod(fe_count1,printff) == printff - 1
                        x_random_normalized_save(fnccall_cnt) = x_random_normalized(ivar);
                        xi_real_save(fnccall_cnt) = xi_real(ivar);
                        s14fig_temp(fnccall_cnt) = s1;
                        s24fig_temp(fnccall_cnt) = s2;
                    end

                    if fe_count1 > 1 && mod(fe_count1,printff/10) == 0 && sn
                        if fe_count1 == parameter.MaxEval
                            set(cb,'YTick',[1 parameter.MaxEval])
                        else
                            set(cb,'YTick',[1 fe_count1 parameter.MaxEval])
                        end
                    end

                    if ((fe_count1 == 1) || (mod(fe_count1,printff) == 0)) && sn
                        num_print_h_fcnt = 0:printff:parameter.MaxEval;                                                   
                        subplot_t = 1;
                        subplot_s = 0.65/(2);
                        if fe_count1 == 1 && fnccall_cnt == 1
                            x_bar4fig_vec = meann;
                            s14fig_vec = shape;
                            s24fig_vec = shape;
                            for hi = 1:n_var
                                xx(hi,:) = h_function(x_bar4fig_vec(hi),s14fig_vec(hi),...
                                    s24fig_vec(hi),xi_s);
                            end

                            for hi = 1:n_var            
                                subplot_y = (hi-0.825)/2.35;
                                subplot_x = 0.55;
                                hAxes(subplot_t) = axes('Units','normalized','Position',...
                                    [subplot_x subplot_y-0.005 subplot_s subplot_s-0.0075]);               
                                subplot_t = subplot_t+1;
                            end

                            set(hAxes_surfplot,'Nextplot','add')
                            plot3(hAxes_surfplot,x_normalized4commandwindow(1),...
                                x_normalized4commandwindow(2),fitness2show(1,1),...
                                'LineWidth',2,'Marker','v')

                            for hi = 1:n_var             

                                set(gcf,'CurrentAxes',hAxes(hi)) 
                                set(hAxes(hi),'Box','on')
                                ph(hi) = plot(xi_s,xx(hi,:),'Parent',hAxes(hi));
                                set(get(hAxes(1),'XLabel'),'String','x_1''  /  X_1'' \in [problem domain] \rightarrow')
                                set(get(hAxes(2),'XLabel'),'String','x_2''  /  X_2'' \in [problem domain] \rightarrow')
                                set(get(hAxes(1),'YLabel'),'String','x_1  /  X_1 \in [problem domain] \rightarrow')
                                set(get(hAxes(2),'YLabel'),'String','x_2  /  X_2 \in [problem domain] \rightarrow')

                                set(hAxes(hi),'XLim',[0 1])
                                set(hAxes(hi),'YLim',[0 1])

                                set(hAxes(hi),'XTick',0:0.5:1,'YTick',0:0.5:1)

                                XTick = get(hAxes(hi),'XTick');
                                hx = get(hAxes(hi),'XLabel');
                                set(hx,'Units','normalized');
                                xpos = get(hx,'Position');
                                ypos = xpos(2);

                                text(XTick(1), ypos, {XTick(1);num2str(parameter.x_min(2))},'HorizontalAlignment','center');
                                text(XTick(2), ypos, {XTick(2);[]},'HorizontalAlignment','center');
                                text(XTick(3), ypos, {XTick(3);num2str(parameter.x_max(2))},'HorizontalAlignment','center');

                                set(gca,'XTickLabel',[])

                                YTick = get(hAxes(hi),'YTick');
                                hy = get(hAxes(hi),'YLabel');
                                set(hy,'Units','normalized');
                                xxpos = get(hy,'Position');
                                yypos = xxpos(1)+0.05;

                                text(yypos, YTick(1), [num2str(parameter.x_min(2)) '    ' num2str(YTick(1))],'HorizontalAlignment','right');
                                text(yypos, YTick(2), ['        ' num2str(YTick(2))],'HorizontalAlignment','right');
                                text(yypos, YTick(3), [num2str(parameter.x_max(2)) '    ' num2str(YTick(3))],'HorizontalAlignment','right');

                                set(gca,'YTickLabel',[])

                                colormap(myColormap)
                                caxis([1 parameter.MaxEval])
                                cb = colorbar([0.91 0.065 0.025 0.745]);
                                set(cb,'YTick',[])
                                set(cb,'YTickMode','manual')
                                set(cb,'YColor',[1 1 1])

                                set(ph(hi),'Color',myColormap(1,:))
                                title1_handle(1) = title(['s_{11}' '=' num2str(s14fig_vec(1),'%4.3f') ...
                                      ' || ' 'mean=' num2str(x_bar4fig_vec(1),'%4.3f') ...
                                      ' || ' 's_{12}' '=' num2str(s24fig_vec(1),'%4.3f')],...
                                    'FontSize',8,'Parent',hAxes(1));

                                title1_handle(2) = title(['s_{21}' '=' num2str(s14fig_vec(2),'%4.3f') ...
                                      ' || ' 'mean=' num2str(x_bar4fig_vec(2),'%4.3f') ...
                                      ' || ' 's_{22}' '=' num2str(s24fig_vec(2),'%4.3f')],...
                                    'FontSize',8,'Parent',hAxes(2));

                                if hi == 2
                                    text(1.3,-0.9,'Iteration-No. or Function-Evaluation-No. \rightarrow',...
                                        'Units','normalized','rotation',90,'FontSize',11,'FontWeight','bold')
                                end

                                hText = findall(gca,'type','text');
                                hUIControls = findall(gcf,'type','uicontrol');
                                set([hText;hUIControls],'units','normalized','fontunits','normalized');
                                if rem(hi,n_var) == 0
                                    drawnow
                                end  
                                
                                
                            end

                        main_axes = gca;
                        
                        end

                        if fe_count1 > 1                                                            
                            s14fig_vec = s14fig_temp;
                            s24fig_vec = s24fig_temp;
                            x_random_normalized_in = x_random_normalized_save;
                            xi_real_in = xi_real_save;

                            x_bar4fig_vec = meann;

                            if (parameter.n_random_ini==1)
                               for hi = 1:parameter.n_random_ini
                                color_sf = round(((find((fe_count1-considered_count(hi).*printff) == ...
                                num_print_h_fcnt))/length(num_print_h_fcnt))*length(myColormap));
                                xx(hi,:) = h_function(x_bar4fig_vec(hi),s14fig_vec(hi),...
                                    s24fig_vec(hi),xi_s);
                               end  
                            else
                            for hi = 1:n_var
                                color_sf = round(((find((fe_count1-considered_count(hi).*printff) == ...
                                num_print_h_fcnt))/length(num_print_h_fcnt))*length(myColormap));
                                xx(hi,:) = h_function(x_bar4fig_vec(hi),s14fig_vec(hi),...
                                    s24fig_vec(hi),xi_s);
                            end  
                            end

                            if fnccall_cnt == nnz(considered(ipp_gbest,:))
                                for hi = 1:n_var            
                                    subplot_y = (hi-0.825)/2.35;
                                    subplot_x = 0.55;
                                    hAxes(subplot_t) = axes('Units','normalized','Position',...
                                        [subplot_x subplot_y-0.005 subplot_s subplot_s-0.0075],'Parent',cf); 
                                    subplot_t = subplot_t+1;
                                end
                                set(hAxes,'Nextplot','add')
                                set(hAxes,'DrawMode','fast')
                                
                                if fe_count1 ~= parameter.MaxEval
                                    plot3(hAxes_surfplot,x_normalized_pic(1),...
                                    x_normalized_pic(2),fitness2show(1,1),...
                                    'LineWidth',2,'Marker','v','MarkerFaceColor',myColormap(color_sf,:),...
                                    'MarkerEdgeColor',myColormap(color_sf,:))
                                else
                                    plot3(hAxes_surfplot,x_normalized_pic(1),...
                                    x_normalized_pic(2),fitness2show(1,1),...
                                    'LineWidth',500,'Marker','v','MarkerFaceColor',myColormap(color_sf,:),...
                                    'MarkerEdgeColor',myColormap(color_sf,:))
                                end

                                for hi = 1:n_var
                                    set(hAxes(hi),'Box','on','Visible','off')
                                    if (parameter.n_random_ini==1)
                                        for ii=1:length(considered(ipp_gbest,:))
                                            if (considered(ipp_gbest,ii)==1)
                                                ph(hi) = plot(xi_s,xx(1,:),'-','Parent',hAxes(hi),'LineWidth',2);
                                                set(hAxes(hi),'XLim',[0 1])
                                                set(hAxes(hi),'YLim',[0 1])
                                            end
                                        end
                                    else
                                        for ii=1:length(considered(ipp_gbest,:))
                                            if (considered(ipp_gbest,ii)==1)
                                                ph(hi) = plot(xi_s,xx(hi,:),'-','Parent',hAxes(hi),'LineWidth',2);
                                                set(hAxes(hi),'XLim',[0 1])
                                                set(hAxes(hi),'YLim',[0 1])
                                            end
                                        end
                                    end

                                    set(ph(hi),'Color',myColormap(color_sf,:))

                                    if strcmp(get(considered_count_handle{hi},'Visible'),'on')
                                        set(considered_count_handle{hi},'Visible','off')
                                    end

                                    if strcmp(get(na_pg_handle{hi},'Visible'),'on')
                                        set(na_pg_handle{hi},'Visible','off')
                                    end

                                    if strcmp(get(fitness_h_axes_handle{hi},'Visible'),'on')
                                        set(fitness_h_axes_handle{hi},'Visible','off')
                                    end

                                    considered_count_handle{hi} = text(0.075,0.85,['X_',num2str(hi) '=' num2str(x_normalized4commandwindow(hi),'%7.6f') ' - Denormalized best value'],...
                                        'FontSize',8,'Parent',hAxes(hi),'Units','normalized');

                                    if strcmp(get(title1_handle(1),'Visible'),'on')
                                        set(title1_handle(1),'Visible','off')
                                    end

                                    if strcmp(get(title1_handle(2),'Visible'),'on')
                                        set(title1_handle(2),'Visible','off')
                                    end

                                    if strcmp(get(Best_fitness_handle,'Visible'),'on')
                                        set(Best_fitness_handle,'Visible','off')
                                    end

                                    if strcmp(get(Best_particle_no_handle,'Visible'),'on')
                                        set(Best_particle_no_handle,'Visible','off')
                                    end

                                    if parameter.n_random_ini==1
                                        if (considered(ipp_gbest,1)==1)
                                            title1_handle(1) = title(['s_{11}' '=' num2str(s14fig_vec(1),'%4.3f') ...
                                                ' || ' 'mean=' num2str(x_bar4fig_vec(1),'%4.3f') ...
                                                ' || ' 's_{12}' '=' num2str(s24fig_vec(1),'%4.3f')],...
                                                'FontSize',8,'Parent',hAxes(1));
                                        else
                                            title1_handle(2) = title(['s_{21}' '=' num2str(s14fig_vec(1),'%4.3f') ...
                                                ' || ' 'mean=' num2str(x_bar4fig_vec(1),'%4.3f') ...
                                                ' || ' 's_{22}' '=' num2str(s24fig_vec(1),'%4.3f')],...
                                        'FontSize',8,'Parent',hAxes(2));
                                        end
                                    
                                    else
                                        
                                        title1_handle(1) = title(['s_{11}' '=' num2str(s14fig_vec(1),'%4.3f') ...
                                            ' || ' 'mean=' num2str(x_bar4fig_vec(1),'%4.3f') ...
                                            ' || ' 's_{12}' '=' num2str(s24fig_vec(1),'%4.3f')],...
                                            'FontSize',8,'Parent',hAxes(1));
                                        
                                        title1_handle(2) = title(['s_{21}' '=' num2str(s14fig_vec(2),'%4.3f') ...
                                            ' || ' 'mean=' num2str(x_bar4fig_vec(2),'%4.3f') ...
                                            ' || ' 's_{22}' '=' num2str(s24fig_vec(2),'%4.3f')],...
                                            'FontSize',8,'Parent',hAxes(2));
                                    end

                                    Best_fitness_handle = text(0.985,1.15,['Best Fitness = ' num2str(table.objective(1,1,ipp_gbest),'%7.6f')],...
                                        'Units','normalized','FontWeight','bold','Parent',main_axes);
                                    
                                    na_pg_handle{hi} = text(0.985,1.25,['Active particles = ' num2str(na_pg,'%7.0f')],...
                                         'Units','normalized','FontWeight','bold','Parent',main_axes);

                                     
                                    if na_pg > 1
                                        Best_particle_no_handle = text(0.985,1.35,['No. of best particle = ' num2str(ipp_gbest,'%7.0f')],...
                                            'Units','normalized','FontWeight','bold','Parent',main_axes);
                                    end

                                    if strcmp(get(title1_handle(1),'Visible'),'off')
                                        set(title1_handle(1),'Visible','on')
                                    end

                                    if strcmp(get(title1_handle(2),'Visible'),'off')
                                        set(title1_handle(2),'Visible','on')
                                    end

                                    if strcmp(get(Best_fitness_handle,'Visible'),'off')
                                        set(Best_fitness_handle,'Visible','on')
                                    end

                                    if strcmp(get(Best_particle_no_handle,'Visible'),'off')
                                        set(Best_particle_no_handle,'Visible','on')
                                    end

                                    if rem(hi,n_var) == 0
                                        drawnow
                                    end  
                                end     
                            end
                        end
                    end
                    
                    x12fict=parameter.x_min + parameter.scaling .* table.bests(1,:,ipp_gbest);
                    x12plot = best_history(:,1);
                    x22plot = best_history(:,2);

                    look4fitnesswindow = findobj('type','figure','Name','Search History');

                    if ((Show_Fit_Value > 1 && fnccall_cnt == nnz(considered(ipp_gbest,:)) ) || fe_count1 == parameter.MaxEval)  
                        FitnessPlot_handle = figure(3);

                        set(FitnessPlot_handle,'NumberTitle','off','ToolBar','none')
                        set(FitnessPlot_handle,'Name','Search History')

                        subplot(2,2,1:2,'Parent',FitnessPlot_handle)
                        plot(iterf,best_results,'LineWidth',2)
                        xlabel('Iteration-No. or Function-Evaluation-No. \rightarrow')
                        ylabel('Fitness-Value \rightarrow')
                        title(['Best Fitness = ' num2str(table.fitness(1,1,ipp_gbest),'%7.6f')],...
                            'FontSize',10,'FontWeight','bold');

                        subplot(2,2,3,'Parent',FitnessPlot_handle)
                        plot(iterf,x12plot,'LineWidth',2)
                        xlabel('Iteration-No. or Function-Evaluation-No. \rightarrow')
                        ylabel('x_1 \rightarrow')
                        title(['X_1 = ' num2str(x12fict(1),'%7.6f') '  - Actual denormalized best value'],...
                            'FontSize',8');

                        subplot(2,2,4,'Parent',FitnessPlot_handle)
                        plot(iterf,x22plot,'LineWidth',2)
                        xlabel('Iteration-No. or Function-Evaluation-No. \rightarrow')
                        ylabel('x_2 \rightarrow')
                        title(['X_2 = ' num2str(x12fict(2),'%7.6f') '  - Actual denormalized best value'],...
                            'FontSize',8);

                        Show_Fit_Value = 1;

                        if fe_count1 ~= parameter.MaxEval
                            uicontrol ('style','pushbutton', ...
                            'units','normalized', ...
                            'position', [0.02 0.0325 0.095 0.045], ...
                            'string','Refresh','Callback',@Show_Fit);
                        end
                        drawnow
                    end
                    

                    look4archivewindow = findobj('type','figure','Name','Archive');

                    if (mod(fe_count1,printff) == 0 && Show_Archive_Value > 1 && fnccall_cnt == nnz(considered(ipp_gbest,:)))
                        archivedata(1:parameter.n_tosave,1) = table.objective(:,1,ipp_gbest);

                        archivedata(1:parameter.n_tosave,2:end) = table.bests(:,:,ipp_gbest);
                        archivedata(parameter.n_tosave+1,2:end) = meann;
                        archivedata(parameter.n_tosave+2,2:end) = variance(ipp_gbest,:);

                        if isempty(look4archivewindow)
                            Archive_fig_handle = figure(4);
                            set(Archive_fig_handle,'Units','normalized');
                            set(Archive_fig_handle,'DoubleBuffer','off')
                            set(Archive_fig_handle,'NumberTitle','off','ToolBar','none')
                            set(Archive_fig_handle,'Name','Archive')
                            set(Archive_fig_handle,'MenuBar','none')
                            set(Archive_fig_handle,'Position',[0.05 0.05 0.255 0.025*(parameter.n_tosave+3)])
                            col_names = {'Fitness','x_1','x_2'};

                            for i = 1:parameter.n_tosave
                                row_names{i} = num2str(i);
                            end

                            row_names{i+1}='mean_i';
                            row_names{i+2}='variance_i';
                        end

                        if ~isempty(look4archivewindow)
                            uitable('Data',archivedata,'ColumnName',col_names,... 
                               'RowName',row_names,'Units','normalized','Position',[0 0 1 1],...
                               'Parent',look4archivewindow);
                            drawnow
                        end
                        
                        if exist('Archive_fig_handle')
                            set(Archive_fig_handle,'CloseRequestFcn',@closefnc)
                        end
                    end

                    if fnccall_cnt >= nnz(considered(ipp_gbest,:))
                        fnccall_cnt = 0;
                    end

                    fnccall_cnt = fnccall_cnt + 1;

                if fe_count1 == parameter.MaxEval
                    set(Freeze,'Visible','off')
                    set(Cancel_button,'Visible','off')
                    Freeze = uicontrol ('style','togglebutton', ...
                        'units','normalized', ...
                        'position', [0.02 0.0325 0.095 0.045], ...
                        'callback','uiresume(gcbf)', ...
                        'string','New Problem/Quit','FontWeight','bold');

                    uiwait(gcf)

                    close all
                    close all hidden
                    clear all
                    clear variables
                    clear classes
                    clc 

                    MVMO_SH_run
                end
                elseif (n_var>2)
                    if isempty(fnccall_cnt)  
                        fnccall_cnt = 1;
                    end

                    look4fitnesswindow = findobj('type','figure','Name','Search History');
                    if ((Show_Fit_Value > 1 && fnccall_cnt == nnz(considered(ipp_gbest,:)) ) || fe_count1 == parameter.MaxEval)
                        FitnessPlot_handle = figure(3);
                        set(FitnessPlot_handle,'NumberTitle','off','ToolBar','none')
                        set(FitnessPlot_handle,'Name','Search History')
                        plot(iterf,best_results,'LineWidth',2)
                        xlabel('Iteration-No. or Function-Evaluation-No. \rightarrow')
                        ylabel('Fitness-Value \rightarrow')
                        title(['Best Fitness = ' num2str(table.fitness(1,1,ipp_gbest),'%7.6f')],...
                            'FontSize',10,'FontWeight','bold');
                        Show_Fit_Value = 1;
                        if fe_count1 ~= parameter.MaxEval
                            uicontrol ('style','pushbutton', ...
                            'units','normalized', ...
                            'position', [0.02 0.0325 0.095 0.045], ...
                            'string','Refresh','Callback',@Show_Fit);
                        end
                        drawnow
                    end

                    if fnccall_cnt >= nnz(considered(ipp_gbest,:))
                        fnccall_cnt = 0;
                    end

                    fnccall_cnt = fnccall_cnt + 1;

                if fe_count1 == parameter.MaxEval
                    set(Freeze,'Visible','off')
                    set(Cancel_button,'Visible','off')
                    Freeze = uicontrol ('style','togglebutton', ...
                        'units','normalized', ...
                        'position', [0.02 0.0325 0.095 0.045], ...
                        'callback','uiresume(gcbf)', ...
                        'string','New Problem/Quit','FontWeight','bold');

                    uiwait(gcf)

                    close all
                    close all hidden
                    clear all
                    clear variables
                    clear classes
                    clc 

                    MVMO_SH_run
                end 
                
            end
        end
    end

    %% Callback to close the archive-window
    function closefnc(hh,jj)
        global Show_Archive_Value
        Show_Archive_Value = 1;
        close 'Archive'
    end

    %% Creation of the respective problem-surface for spatially plotting
    function [f,problemname] = fitness_evaluation_surfplot(x_in_probemspace,testcase)
        n = 2;
        global alpha12 beta
        %%
        x1 = x_in_probemspace;

        %% Two-dimensional Test-Functions
        switch testcase
            case 1
                %% De Jong’s function in 2D, f(x; y) = x.^2 + y.^2
                % -5.12 <= xi <= 5.12
                % Minimum fmin==0 at x==0,y==0

                f = 0;
                for i = 1:n
                    f = f + x1{i}.^2;
                end

                problemname = 'De Jong’s function" - f_{gmin}=0 for Xi=0, i=1,...,n';

            case 2
                %% Axis parallel hyper-ellipsoid function
                % -5.12 <= xi <= 5.12
                % Minimum fmin==0 at x==0,y==0

                f = 0;
                for i = 1:n
                    f = f + i.*x1{i}.^2;
                end

                problemname = 'Axis parallel hyper-ellipsoid function" - f_{gmin}=0 for Xi=0, i=1,...,n';

            case 3
                %% Rotated hyper-ellipsoid function
                % -65536 <= xi <= 65536
                % Minimum fmin==0 at x==0,y==0
                f = 0;
                for i = 1:n
                    for j = 1:i
                        f = f + x1{j}.^2;
                    end
                end

                problemname = 'Rotated hyper-ellipsoid function" - f_{gmin}=0 for Xi=0, i=1,...,n';

            case 4
                %% Rosenbrock’s valley
                % -2.048 <= xi <= 2.048
                % Minimum fmin==0 at x==1,y==1

                f = 0;
                for i = 1:n-1
                    f = f + ( 100.*(x1{i+1} - x1{i}.^2).^2 + (1 - x1{i}).^2 );
                end

                problemname = 'Rosenbrock’s valley" - f_{gmin}=0 for Xi=1, i=1,...,n';

            case 5
                %% Rastrigin’s function
                % -5.12 <= xi <= 5.12
                % Minimum fmin==0 at x==0,y==0

                f = 0;
                for i = 1:n
                    f = f + (x1{i}.^2 - 10 .* cos(2.*pi.*x1{i}));
                end
                f = 10 .* n + f;

                problemname = 'Rastrigin’s function" - f_{gmin}=0 for Xi=0, i=1,...,n';

            case 6
                %% Schwefel’s function
                % -500 <= xi <= 500
                % Minimum fmin==-418.9829.*n (here, n == 2) at x==420.9687,y==420.9687

                f = 0;
                for i = 1:n
                    f = f + ( -x1{i} .* sin(sqrt(abs(x1{i}))) );
                end

                problemname = 'Schwefel’s function" - f_{gmin}=-837.9658 for Xi=420.9687, i=1,...,n';

            case 7
                %% Griewangk’s function 
                % -600 <= xi <= 600
                % Minimum fmin==0 at x==0,y==0

                f_prod = 0; f_sum = 0; f = 0;
                for i = 1:n
                    f_prod = f_prod .* (cos(x1{i}./sqrt(i)) + 1);
                end

                for i = 1:n
                    f_sum = f_sum + x1{i}.^2;
                end

                f = 1./4000 .* f_sum - f_prod;

                problemname = 'Griewangk’s function" - f_{gmin}=0 for Xi=0, i=1,...,n';

            case 8
                %% Sum of different power functions
                % -1 <= xi <= 1
                % Minimum fmin==0 at x==0,y==0

                f = 0;
                for i = 1:n
                    f = f + abs(x1{i}).^(i+1);
                end            

                problemname = 'Sum of different power functions" - f_{gmin}=0 for Xi=0, i=1,...,n';

            case 9
                %% Ackley’s function
                % -32.768 <= xi <= 32.768
                % Minimum fmin==0 at x==0,y==0

                n = 2; a = 20; b = 0.2; c = 2.*pi; s1 = 0; s2 = 0; f = 0;
                for i = 1:n;
                   s1 = s1+x1{i}.^2;
                   s2 = s2+cos(c.*x1{i});
                end
                f = -a.*exp(-b.*sqrt(1./n.*s1))-exp(1./n.*s2)+a+exp(1);

                problemname = 'Ackley’s function" - f_{gmin}=0 Xi=0, i=1,...,n';

            case 10
                %% Michalewicz’s function
                % 0 <= xi <= pi
                % Minimum fmin==-1.8013 for m == 10 at x==2.2029,y==1.5708

                f = 0; m = 10;
                for i = 1:n
                    f = f + ((sin(x1{i}) .* (sin((i.*x1{i}.^2)./pi)).^(2.*m)));
                end
                f = -f;

                problemname = 'Michalewicz’s function" - f_{gmin}=-1.8013 for X_1=2.2029, X_2=1.5708'; %ACHTUNG: Verify for n=30

            case 11
                %% Branin’s function
                % -5 <= x <= 10 and  0 <= y <= 15
                % Three equal global minima fmin==0.397887 at: 
                % (1)   x==-pi,y==12.275
                % (2)   x==pi, y==2.275
                % (3)   x==9.42478,y==2.475

                a = 1;  b = 5.1./(4.*pi.^2);   c = 5./pi;   d = 6;  e = 10; fi = 1./(8.*pi);
                f = a .* (x1{2} - b .* x1{1} .^ 2 + c .* x1{1} - d).^2 + e .* (1 - fi) .* cos(x1{1}) + e;

                problemname = 'Branin’s function" - f_{gmin}=0.397887 for three different (X_1,X_2)'; %ACHTUNG: Verify for n=30

            case 12
                %% Easom’s function
                % -100 <= xi <= 100
                % Minimum fmin==-1 at x==pi,y==pi

                f = -cos(x1{1}) .* cos(x1{2}) .* exp(-(x1{1}-pi).^2 - (x1{2}-pi).^2);

                problemname = 'Easom’s function" - f_{gmin}=-1 for X_1=pi, X_2=pi'; %ACHTUNG: Verify for n=30

            case 13
                %% Goldstein-Price’s function
                % -2 <= xi <= 2
                % Minimum fmin==3 at x==0,y==-1

                a = 1+(x1{1}+x1{2}+1).^2.*(19-14.*x1{1}+3.*x1{1}.^2-14.*x1{2}+6.*x1{1}.*x1{2}+3.*x1{2}.^2);
                b = 30+(2.*x1{1}-3.*x1{2}).^2.*(18-32.*x1{1}+12.*x1{1}.^2+48.*x1{2}-36.*x1{1}.*x1{2}+27.*x1{2}.^2);
                f = a .* b;

                problemname = 'Goldstein-Price’s function" - f_{gmin}=3 for X_1=0, X_2=-1';  %ACHTUNG: Verify for n=30

            case 14
                %% “Drop wave” function
                % -5.12 <= xi <= 5.12

                f = - ( (1 + cos(12.*sqrt(x1{1}.^2 + x1{2}.^2))) ./ (0.5.* (x1{1}.^2 + x1{2}.^2)+ 2) );   

                problemname = 'Drop wave function" - f_{gmin}=-1 for X_1=0, X_2=0'; %ACHTUNG: Verify for n=30

            case 15
                %% Shubert function
                % -5.12 <= xi <= 5.12
                % 18 global minima  at fmin = -186.7309
                s1 = 0;     s2 = 0;
                for i = 1:5;   
                    s1 = s1+i.*cos((i+1).*x1{1}+i);
                    s2 = s2+i.*cos((i+1).*x1{2}+i);
                end
                f = s1.*s2;   

                problemname = 'Shubert function" - f_{gmin}=-186.7309, several different solutions'; %ACHTUNG: Verify for n=30
                
            case 16
                %% Deceptive function Type III
                % 0 <= xi <= 1
                % Minimum fmin==0 at x==random,y==random (both random values from [0 1])
                alpha12 = rand(1,n);
                gg = x1; 
                beta = 1;
                sum1 = 0;
                
                for i = 1:n
                    for j = 1:numel(x1{i})
                        if x1{i}(j) >= 0 && x1{i}(j) <= ((4/5) * alpha12(i))
                            gg{i}(j)=-((x1{i}(j)) ./ alpha12(i)) + 4/5;
                        elseif x1{i}(j) > ((4/5) * alpha12(i)) && x1{i}(j) <= alpha12(i)
                            gg{i}(j)=((5 .* (x1{i}(j))) ./ alpha12(i)) - 4;
                        elseif x1{i}(j) > alpha12(i) && x1{i}(j) <= ((1 + 4 .* alpha12(i)) / 5)
                            gg{i}(j)= ((5 .* ((x1{i}(j)) - alpha12(i))) ./ (alpha12(i) - 1)) + 1;
                        elseif x1{i}(j) > ((1 + 4 .* alpha12(i)) / 5) && x1{i}(j) <= 1
                            gg{i}(j)=((x1{i}(j) - 1) ./ (1 - alpha12(i))) + 4/5;
                        end
                    end
                end
                
                for i = 1:n
                    sum1 = sum1 + gg{i};
                end
                
                f = ((sum1 ./ n)).^beta;
                
                problemname = ['Deceptive function Type III" - f_{gmin}=0 for random (X_1,X_2), chosen \beta = ' num2str(beta) '(linear)'];
        
            case 17
                %% Deceptive function Type III
                % 0 <= xi <= 1
                % Minimum fmin==0 at x==random,y==random (both random values from [0 1])
                alpha12 = rand(1,n);
                gg = x1; 
                beta = 0.1;
                sum1 = 0;
                
                for i = 1:n
                    for j = 1:numel(x1{i})
                        if x1{i}(j) >= 0 && x1{i}(j) <= ((4/5) * alpha12(i))
                            gg{i}(j)=-((x1{i}(j)) ./ alpha12(i)) + 4/5;
                        elseif x1{i}(j) > ((4/5) * alpha12(i)) && x1{i}(j) <= alpha12(i)
                            gg{i}(j)=((5 .* (x1{i}(j))) ./ alpha12(i)) - 4;
                        elseif x1{i}(j) > alpha12(i) && x1{i}(j) <= ((1 + 4 .* alpha12(i)) / 5)
                            gg{i}(j)= ((5 .* ((x1{i}(j)) - alpha12(i))) ./ (alpha12(i) - 1)) + 1;
                        elseif x1{i}(j) > ((1 + 4 .* alpha12(i)) / 5) && x1{i}(j) <= 1
                            gg{i}(j)=((x1{i}(j) - 1) ./ (1 - alpha12(i))) + 4/5;
                        end
                    end
                end
                
                for i = 1:n
                    sum1 = sum1 + gg{i};
                end
                
                f = ((sum1 ./ n)).^beta;
                
                problemname = ['Deceptive function Type III" - f_{gmin}=0 for random (X_1,X_2), chosen \beta = ' num2str(beta) '(non-linear)'];
                
        end
    end

    %% Problem domain
    function boundaries(testcase)
        global parameter
        switch testcase
            case {1,2,5,14,15}
                min = -5.12;        max = 5.12;
            case 3
                min = -65536;       max = 65636;
            case 4
                min = -2.048;       max = 2.048;
            case 6
                min = -500;         max = 500;
            case 7
                min = -600;         max = 600;
            case 8
                min = -1;           max = 1;
            case 9
                min = -32.768;      max = 32.768;
            case 10
                min = 0;            max = pi;
            case 11
                min = -5;           max = 15;
            case 12
                min = -100;         max = 100;
            case 13
                min = -2;           max = 2;
            case {16,17}
                min = 0;            max = 1;
        end

        parameter.x_min = min * ones(1,parameter.n_var);
        parameter.x_max = max * ones(1,parameter.n_var);

        parameter.v_min = min * ones(1,parameter.n_var);
        parameter.v_max = max * ones(1,parameter.n_var);
    end

    %% Fitness evaluation
    function [f,xn] = fitness_evaluation(xn,testcase)
        global parameter
        global alpha12 beta
        n = parameter.n_var;%2;

        %% Retransformation to the original space for fitness-evaluation
        x1 = parameter.x_min + parameter.scaling .* xn;
        
        %% Two-dimensional Test-Functions
        switch testcase
            case 1
                %% De Jong’s function in 2D, f(x; y) = x.^2 + y.^2
                % -5.12 <= xi <= 5.12
                % Minimum fmin==0 at x==0,y==0

                f = 0;
                for i = 1:n
                    f = f + x1(i).^2;
                end

            case 2
                %% Axis parallel hyper-ellipsoid function
                % -5.12 <= xi <= 5.12
                % Minimum fmin==0 at x==0,y==0

                f = 0;
                for i = 1:n
                    f = f + i.*x1(i).^2;
                end

            case 3
                %% Rotated hyper-ellipsoid function
                % -65536 <= xi <= 65536
                % Minimum fmin==0 at x==0,y==0
                f = 0;
                for i = 1:n
                    for j = 1:i
                        f = f + x1(j).^2;
                    end
                end

            case 4
                %% Rosenbrock’s valley
                % -2.048 <= xi <= 2.048
                % Minimum fmin==0 at x==1,y==1

                f = 0;
                for i = 1:n-1
                    f = f + ( 100.*(x1(i+1) - x1(i).^2).^2 + (1 - x1(i)).^2 );
                end

            case 5
                %% Rastrigin’s function
                % -5.12 <= xi <= 5.12
                % Minimum fmin==0 at x==0,y==0

                f = 0;
                for i = 1:n
                    f = f + (x1(i).^2 - 10 .* cos(2.*pi.*x1(i)));
                end
                f = 10 .* n + f;

            case 6
                %% Schwefel’s function
                % -500 <= xi <= 500
                % Minimum fmin==-418.9829.*n (here, n == 2) at x==420.9687,y==420.9687

                f = 0;
                for i = 1:n
                    f = f + ( -x1(i) .* sin(sqrt(abs(x1(i)))) );
                end

            case 7
                %% Griewangk’s function 
                % -600 <= xi <= 600
                % Minimum fmin==0 at x==0,y==0

                f_prod = 0; f_sum = 0; f = 0;
                for i = 1:n
                    f_prod = f_prod .* (cos(x1(i)./sqrt(i)) + 1);
                end

                for i = 1:n
                    f_sum = f_sum + x1(i).^2;
                end

                f = 1./4000 .* f_sum - f_prod;

            case 8
                %% Sum of different power functions
                % -1 <= xi <= 1
                % Minimum fmin==0 at x==0,y==0

                f = 0;
                for i = 1:n
                    f = f + abs(x1(i)).^(i+1);
                end            

            case 9
                %% Ackley’s function
                % -32.768 <= xi <= 32.768
                % Minimum fmin==0 at x==0,y==0

                n = 2; a = 20; b = 0.2; c = 2.*pi; s1 = 0; s2 = 0; f = 0;
                for i = 1:n;
                   s1 = s1+x1(i).^2;
                   s2 = s2+cos(c.*x1(i));
                end
                f = -a.*exp(-b.*sqrt(1./n.*s1))-exp(1./n.*s2)+a+exp(1);

            case 10
                %% Michalewicz’s function
                % 0 <= xi <= pi
                % Minimum fmin==-1.8013 for m == 10 at x==2.2029,y==1.5708

                f = 0; m = 10;
                for i = 1:n
                    f = f + ((sin(x1(i)) .* (sin((i.*x1(i).^2)./pi)).^(2.*m)));
                end
                f = -f;

            case 11
                %% Branin’s function
                % -5 <= x <= 10 and  0 <= y <= 15
                % Three equal global minima fmin==0.397887 at: 
                % (1)   x==-pi,y==12.275
                % (2)   x==pi, y==2.275
                % (3)   x==9.42478,y==2.475

                a = 1;  b = 5.1./(4.*pi.^2);   c = 5./pi;   d = 6;  e = 10; fi = 1./(8.*pi);
                f = a .* (x1(2) - b .* x1(1) .^ 2 + c .* x1(1) - d).^2 + e .* (1 - fi) .* cos(x1(1)) + e;

            case 12
                %% Easom’s function
                % -100 <= xi <= 100
                % Minimum fmin==-1 at x==pi,y==pi

                f = -cos(x1(1)) .* cos(x1(2)) .* exp(-(x1(1)-pi).^2 - (x1(2)-pi).^2);

            case 13
                %% Goldstein-Price’s function
                % -2 <= xi <= 2
                % Minimum fmin==3 at x==0,y==-1

                a = 1+(x1(1)+x1(2)+1).^2.*(19-14.*x1(1)+3.*x1(1).^2-14.*x1(2)+6.*x1(1).*x1(2)+3.*x1(2).^2);
                b = 30+(2.*x1(1)-3.*x1(2)).^2.*(18-32.*x1(1)+12.*x1(1).^2+48.*x1(2)-36.*x1(1).*x1(2)+27.*x1(2).^2);
                f = a .* b;

            case 14
                %% “Drop wave” function
                % -5.12 <= xi <= 5.12

                f = - ( (1 + cos(12.*sqrt(x1(1).^2 + x1(2).^2))) ./ (0.5 .* (x1(1).^2 + x1(2).^2) + 2) );   

            case 15
                %% Shubert function
                % -5.12 <= xi <= 5.12
                % Minimum fmin = -186.7309
                s1 = 0;     s2 = 0;
                for i = 1:5;   
                    s1 = s1+i.*cos((i+1).*x1(1)+i);
                    s2 = s2+i.*cos((i+1).*x1(2)+i);
                end
                f = s1.*s2;
                
            case {16,17}
                %% Deceptive function Type III
                % 0 <= xi <= 1
                % Minimum fmin==0 at x==alpha12(1),y==alpha12(2) (both random values from [0 1])
                gg = zeros(1,n); sum1 = 0;
                for i = 1:n
                    %for j = 1:numel(x1{i})
                        if x1(i) >= 0 && x1(i) <= ((4/5) * alpha12(i))
                            gg(i)=-((x1(i)) / alpha12(i)) + 4/5;
                        elseif x1(i) > ((4/5) * alpha12(i)) && x1(i) <= alpha12(i)
                            gg(i)=((5 * (x1(i))) / alpha12(i)) - 4;
                        elseif x1(i) > alpha12(i) && x1(i) <= ((1 + 4 * alpha12(i)) / 5)
                            gg(i)= ((5 * ((x1(i)) - alpha12(i))) / (alpha12(i) - 1)) + 1;
                        elseif x1(i) > ((1 + 4 * alpha12(i)) / 5) && x1(i) <= 1
                            gg(i)=((x1(i) - 1) / (1 - alpha12(i))) + 4/5;
                        end
                    %end
                end
                
                for i = 1:n
                    sum1 = sum1 + gg(i);
                end
                
                f = ((sum1 ./ n)).^beta;
        end
    end

    %% Control for freezing the figure
    function Freeze_Figure(obj,ll)
        val = get(obj,'Value');

        if val == 1
            set(obj,'String','Continue Evaluation')
            pause
        elseif val == 0
            inputemu('key_normal','\ENTER')
            set(obj,'String','Freeze Figure')
        end
    end

    %% Cancel execution
    function my_finishdlg(ii,hh)
        button = questdlg('Do you want to cancel execution?', ...
                  'Cancel','Yes','No','No');
        switch button
          case 'Yes'
            close all 
          case 'No',
            quit cancel;
        end
    end

    %% Setting the MVMO-SH-inherent parameters via input dialog
    function Algorithm_Settings_Dialogue(hObject,eventdata)
       
        global parameter
        global ratio_gute_max ratio_gute_min

        prompt = {'Particle setting [Number, Initial percentage ''good'', Final percentage ''good'']: ',...
                  'Archive size (2-5): ',...
                  'Selection mode (1-5): ',...
                  'No. of variables selected for mutation [Initial(<total), Final(<initial)]: ',...
                  '[fs_factor_start(<=1), fs_factor_end(>=1)]: ',...
                  '[dddd(1-5), delta_dddd_start(0.01-0.4), delta_dddd_end(<delta_dddd_start)]: ',...
                  'Decrement strategy of variables selected for mutation [1,2,3,41,42,43]:'};
        
        dlg_title = 'Algorithm Settings';
        num_lines = 1;
        def = {'[1,0.7,0.2]','4','4','[2,1]','[1,1]','[1,0.05,0.05]','41'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        varaux1 = str2num(answer{1});
        parameter.n_par=varaux1(1);
        ratio_gute_max=varaux1(2);
        ratio_gute_min=varaux1(3);
        
        while ratio_gute_max < ratio_gute_min
            errordlg('Final percentage ''good'' must be smaller than initial percentage ''good''.','Error','normal')
            uiwait(gcf);
            Algorithm_Settings_Dialogue;
        end
        
        while round(ratio_gute_min*parameter.n_par) < 3 && parameter.n_par > 1
            errordlg('''round(Final percentage ''good'' x Number of particles)'' must be greater than 3.','Multi-parent operation','normal')
            uiwait(gcf);
            Algorithm_Settings_Dialogue;
        end
        
        parameter.n_tosave = str2double(answer{2});
        
        parameter.mode = str2double(answer{3});
        varaux2=str2num(answer{4});
        parameter.n_random_ini = varaux2(1);
        parameter.n_random_last = varaux2(2);
        varaux3=str2num(answer{5}); 
        parameter.fs_factor_start = varaux3(1);
        parameter.fs_factor_end = varaux3(2);
        varaux4=str2num(answer{6});
        parameter.dddd = varaux4(1);
        parameter.delta_dddd_start = varaux4(2);
        parameter.delta_dddd_end = varaux4(3);
        parameter.nrand_method = str2double(answer{7});
    end

%% Optionally setting the problem dimensionality (if>2: -> without graphical representation)
    function Procedure_Settings_Dialogue(hObject,eventdata)
        global parameter
        global printff
        global testcase local_search0_percentage
        
        prompt = {'Maximum number of function evaluations: ',...
                  'Print step size: ',...
                  'Probability of local search (percentage / number of optimization variables): ',...
                  'Number of optimization variables (2 up to 30): '};
        dlg_title = 'Procedure Settings';
        num_lines = 1;
        def = {'1000','25','0','2'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);

        parameter.MaxEval = str2double(answer{1});
        printff = str2double(answer{2});
        local_search0_percentage = str2double(answer{3});
        parameter.n_var = str2double(answer{4});
 
        while parameter.MaxEval > 100000
            errordlg('For reasons of didactics, in this demo the maximum number of function evaluations is limited to 100000!.','Error','normal')
            uiwait(gcf);
            Procedure_Settings_Dialogue;
        end

        while parameter.MaxEval < printff
            errordlg('The step size has to be smaller than the maximum number of function evaluations!','Error','normal')
            uiwait(gcf);
            Procedure_Settings_Dialogue;
        end

        while parameter.MaxEval / printff > 200
            errordlg('The granularity is too large. Adjust the step size and/or the maximum number of function evaluations!','Error','normal')
            uiwait(gcf);
            Procedure_Settings_Dialogue;
        end
        
        while ((parameter.n_var < 2) || (parameter.n_var > 30)) %|| ((parameter.n_var > 2) & (parameter.n_var < 30))
            errordlg('Problem dimension should be between 2-30!.','Error','normal')
            uiwait(gcf);
            Procedure_Settings_Dialogue;
        end 
        boundaries(testcase);
    end

%% Popup menu items for problem selection
    function Prob_Selection(hObj,event)
        global testcase
        val = get(hObj,'Value');
        global parameter
        global hAxes_surfplot
        global X Y Z
        
        switch val
            case 3 % 'De Jong’s function'
                testcase = 1;
            case 4 % 'Axis parallel hyper-ellipsoid function'
                testcase = 2;
            case 5 % 'Rotated hyper-ellipsoid function'
                testcase = 3;
            case 6 % 'Rosenbrock’s valley'
                testcase = 4;
            case 7 % 'Rastrigin’s function'
                testcase = 5;
            case 8 % 'Schwefel’s function'
                testcase = 6;
            case 9 % 'Griewangk’s function'
                testcase = 7;
            case 10 % 'Sum of different power functions'
                testcase = 8;
            case 11 % 'Ackley’s function'
                testcase = 9;
            case 12 % 'Michalewicz’s function'
                testcase = 10;
            case 13 % 'Branin’s function'
                testcase = 11;
            case 14 % 'Easom’s function'
                testcase = 12;
            case 15 % 'Goldstein-Price’s function'
                testcase = 13;
            case 16 % '“Drop wave” function'
                testcase = 14;
            case 17 % 'Shubert function'
                testcase = 15;
            case 18 % 'Deceptive function Type III linear'
                testcase = 16;
            case 19 % 'Deceptive function Type III non-linear'
                testcase = 17;
            otherwise
                testcase = 1;
        end
 
            boundaries(testcase);
            x_surfplot = linspace(parameter.x_min(1),parameter.x_max(1),200);
            y_surfplot = x_surfplot;
            [X Y] = meshgrid(x_surfplot,y_surfplot);
            
            [Z,pname] = fitness_evaluation_surfplot({X Y},testcase);
            sp = surf(hAxes_surfplot,X,Y,Z);
            xlabel('X_1 \rightarrow')
            ylabel('X_2 \rightarrow')
            zlabel('Fitness or Objective Value \rightarrow')
            title(['Problem: ' '"' pname],'Fontsize',12,'Fontweight','bold')
            set(sp,'AlphaDataMapping','none')
            set(sp,'EdgeColor','none','Facecolor','texture')
            
            set(sp,'FaceLighting','phong',...
              'FaceColor','interp',...
              'EdgeColor','none',...
              'BackFaceLighting','reverselit',...
              'CDataMapping','direct')
            
            light('Position',[1 3 2]);
            light('Position',[-3 -1 3]);
            material shiny
            alpha(0.3)
            
            text(0.75,1.225,'MVMO-SH Search Progress Visualization',...
                    'Units','normalized','FontSize',16,'FontWeight','bold') 
                
            text(0.38,1.15,' Application of n-dimensional Continuous and Discontinuous Unconstrained Minimization Problems',...
                    'Units','normalized','FontSize',12,'FontWeight','bold');
    end

    %% Popup menu items for Colormap alignment
    function setmap(hObj,event)       
        global myColormap    
        val = get(hObj,'Value');
        switch val
            case 3
                myColormap = colormap('jet');
            case 4
                myColormap = colormap('hsv');
            case 5
                myColormap = colormap('hot');
            case 6
                myColormap = colormap('cool');
            case 7
                myColormap = colormap('autumn');
            case 8
                myColormap = colormap('winter');
            otherwise
                myColormap = colormap('autumn');
        end
        myColormap = flipud(myColormap);
    end

    %% Control search history-window containing fitness- and shape-factor-charts
    function Show_Fit(hObject,eventdata)
        global Show_Fit_Value
        Show_Fit_Value = Show_Fit_Value + 1;
    end

    %% Control archive-window
    function Show_Archive(hObject,eventdata)
        global Show_Archive_Value
        Show_Archive_Value = Show_Archive_Value + 1;
    end

    %% Emulation of 'Return'-key to terminate MATLAB-'Pause'-mode
    function inputemu(varargin)
        %INPUTEMU   Java-Based Mouse/Keyboard Emulator
        %   INPUTEMU(ACTION,PARAM) emulates the human input via mouse and
        %   keyboard. This utility uses Java java.awt.Robot class.
        %
        %   INPUTEMU(MOUSE_ACTION,[X Y WHEEL]) performs mouse emulation with 5
        %   action types:
        %
        %      'move'|'none'|'wheel' - No mouse click (default)
        %      'normal'              - Click left mouse button
        %      'extend'              - Shift-click left mouse button
        %      'alternate'           - Control-click left mouse button
        %      'open'                - Double click any mouse button
        %
        %   The mouse cursor is first moved to the coordinate [X Y] (with respect
        %   to the lower-left corner of the primary screen), followed by button
        %   click action (if specified), then it turns the mouse wheel by WHEEL
        %   notches. The number of "notches" to move the mouse wheel negative
        %   values indicate movement up/away from the user, positive values
        %   indicate movement down/towards the user. If PARAM (or any of its
        %   elements) is missing, the default position is the current cursor
        %   position and WHEEL = 0.
        %
        %   INPUTEMU('wheel',WHEEL) may be used to turn wheel without moving mouse
        %   cursor.
        %
        %   INPUTEMU(KEYBOARD_ACTION,'text') performs keyboard emulation with
        %   4 action types:
        %
        %      'key_normal'  - Normal keypress (shift-key pressed as needed)
        %      'key_ctrl'    - CONTROL-key keypress
        %      'key_alt'     - ALT-key keypress
        %      'key_win'     - WIN-key keypress
        %
        %   The 'text' to be typed may contain special keys as escape characters
        %   with '\' prefix:
        % 
        %      '\BACKSPACE'   '\TAB'         '\ENTER'       '\SHIFT'
        %      '\CTRL'        '\ALT'         '\PAUSE'       '\CAPSLOCK'
        %      '\ESC'         '\PAGEUP'      '\PAGEDOWN'    '\END'
        %      '\HOME'        '\LEFT'        '\UP'          '\RIGHT'
        %      '\DOWN'        '\PRINTSCREEN' '\INSERT'      '\DELETE'
        %      '\WINDOWS'     '\NUMLOCK'     '\SCROLLLOCK'  
        %
        %   These escape characters are NOT case sensitive while regular characters
        %   are. For regular backslash, use '\\' unless it is the only character
        %   then '\' may be used. 
        %
        %   In addition to above action types, there are 8 low-level actions to
        %   specify button or key to be down (pressed) or up (released).
        %
        %      'left_down'/'left_up' | 'drag_on'/'drag_off' - Left mouse button
        %      'right_down'/'right_up'                      - Right mouse button
        %      'middle_down'/'middle_up'                    - Middle mouse button
        %      'key_down'/'key_up'                          - Keyboard key
        %
        %   INPUTEMU(CMDS) can be used to perform multiple commands. CMDS is an
        %   2-by-N cell array where the n-th column {ACTION_n;PARAM_n} defines the
        %   n-th input action. CMDS may be given in an N-by-2 array (if N~=2)
        %
        %   INPUTEMU(CMDS,T) performs the CMDS action sequence with T second
        %   update period.
        %
        %   INPUTEMU(CMDS,Tints) where time is a N-element vector, specifies
        %   update interval (in s) for individual action. time(1) specifies the lag
        %   before performing the first action.
        %
        %   See also jmouseemu and inputlog.

        % Version 1.0 (Aug. 31, 2010)
        % Written by: Takeshi Ikuma
        % Created: Aug. 31, 2010
        % Revision History:
        %  - (Aug. 31, 2010) : initial release

        % To Do List
        % * Support for pause action
        % * Allow a choice between block/non-block modes
        % * Internationalization
        %   - support for user specified locale / OemChars def
        %   - support for special keys not on the US keyboard
        %   - support for non-European alphabet support (possible?)

        % Java environment check
        error(javachk('jvm'));
        error(javachk('awt'));

        % default (i.e., US|us locale) OEM character map
        import java.awt.event.KeyEvent
        defaultOemChars = {... % {'ch w/o shift' 'ch w/shift' Java virtual keycode}
           '`' '~'   KeyEvent.VK_BACK_QUOTE
           '1' '!'   KeyEvent.VK_1
           '2' '@'   KeyEvent.VK_2
           '3' '#'   KeyEvent.VK_3
           '4' '$'   KeyEvent.VK_4
           '5' '%'   KeyEvent.VK_5
           '6' '^'   KeyEvent.VK_6
           '7' '&'   KeyEvent.VK_7
           '8' '*'   KeyEvent.VK_8
           '9' '('   KeyEvent.VK_9
           '0' ')'   KeyEvent.VK_0
           '-' '_'   KeyEvent.VK_MINUS
           '=' '+'   KeyEvent.VK_EQUALS
           '[' '{'   KeyEvent.VK_OPEN_BRACKET
           ']' '}'   KeyEvent.VK_CLOSE_BRACKET
           '\' '|'   KeyEvent.VK_BACK_SLASH
           ';' ':'   KeyEvent.VK_SEMICOLON
           '''' '"'  KeyEvent.VK_QUOTE
           ',' '<'   KeyEvent.VK_COMMA
           '.' '>'   KeyEvent.VK_PERIOD
           '/' '?'   KeyEvent.VK_SLASH
        };

        % parse input arguments
        [cmds,deltaT,idx_mouse,oemChars,msg] = parseInput(defaultOemChars,nargin,varargin);
        if ~isempty(msg), error(msg); end

        % expand upper-level actions
        [actions,params,dt,msg] = convertToLowLevelActions(cmds,deltaT,idx_mouse,oemChars);
        if ~isempty(msg), error(msg); end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        import java.awt.Robot
        import java.awt.event.InputEvent

        % initialize the java.robot
        robot = Robot();

        for n = 1:size(actions)
            % cue the action
            switch actions{n}
                case 'move'
                    robot.mouseMove(params{n}(1),params{n}(2));
                case 'wheel'
                    robot.mouseWheel(params{n});
                case 'left_down'
                    robot.mousePress(InputEvent.BUTTON1_MASK);
                case 'left_up'
                    robot.mouseRelease(InputEvent.BUTTON1_MASK);
                case 'right_down'
                    robot.mousePress(InputEvent.BUTTON3_MASK);
                case 'right_up'
                    robot.mouseRelease(InputEvent.BUTTON3_MASK);
                case 'middle_down'
                    robot.mousePress(InputEvent.BUTTON2_MASK);
                case 'middle_up'
                    robot.mouseRelease(InputEvent.BUTTON2_MASK);
                case 'key_down'
                    robot.keyPress(params{n});
                case 'key_up'
                    robot.keyRelease(params{n});
            end

            % delay
            robot.delay(dt(n));
        end

        % wait until robot's done
        robot.waitForIdle();

%         end % jinputemu

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Parse Input Arguments
        function [cmds,deltaT,idx_mouse,oemChars,msg] = parseInput(defaultOemChars,n,argin)

        cmds = {};
        deltaT = [];
        idx_mouse = [];
        oemChars = defaultOemChars;

        msg = nargchk(1,2,n); %#ok
        if ~isempty(msg), return; end

        if ~iscell(argin{1}) % single-command mode
            switch n
                case 1
                    cmds = [argin(1);cell(1)];
                case 2
                    cmds = argin;
                otherwise
                    cmds = {}; % never gets here
            end
        else % multiple-command mode
            cmds = argin{1};
            if n>1
                deltaT = argin{2};
            end
        end

        % get dimension
        dim = size(cmds);

        % use the first size-2 dimension as the command
        idx = find(dim==2,1,'first');
        if ~isempty(idx) && idx~=1
            idx = [idx 1:idx-1 idx+1:ndims(cmds)];
            cmds = permute(cmds,idx);
            dim(:) = dim(idx);
        end

        % check cmds
        if dim(1)~=2
            msg = 'One of the sizes of the CMDS matrix must must be 2.';
            return;
        end

        mouse_actions = {'none','move','wheel','normal','extend','alternate',...
            'open','drag_on','drag_off','left_down','left_up','middle_down',...
            'middle_up','right_down','right_up'};
        keyboard_actions = {'key_normal','key_alt','key_ctrl','key_win','key_up','key_down'};

        if any(cellfun(@(action)~isempty(action) && ~ischar(action),cmds(1,:)))
            msg = 'Command actions (CMDS{1,:}) must be a string of characters.';
            return;
        end

        idx_mouse = cellfun(@(action)any(strcmpi(action,mouse_actions)),cmds(1,:));
        idx_keyboard = cellfun(@(action)any(strcmpi(action,keyboard_actions)),cmds(1,:));

        if any(~any(idx_mouse|idx_keyboard))
            msg = 'At least one command action is invalid.';
            return;
        end

        if any(cellfun(@(x)~isempty(x) && (~isnumeric(x) || ~any(numel(x)==[1 2 3]) || any(isinf(x)) || any(isnan(x))),cmds(2,idx_mouse)))
            msg = 'Position must be a 2- or 3-element vector with finite values';
            return;
        end

        if any(cellfun(@(x)~ischar(x),cmds(2,idx_keyboard)))
            msg = 'Text must be a character string';
            return;
        end

        % format mouse parameters
        scrnsz = get(0,'ScreenSize');
        ptrloc = get(0,'PointerLocation');
        for n = find(idx_mouse)
            switch numel(cmds{2,n})
                case 0
                    cmds{2,n}=[ptrloc 0];
                case 1
                    if strcmpi(cmds{1,n},'wheel')
                        cmds{2,n}=[ptrloc cmds{2,n}];
                    else
                        cmds{2,n}=[cmds{2,n} ptrloc(2) 0];
                    end
                case 2
                    cmds{2,n}=[cmds{2,n}(:)' 0];
            end
            if all(cmds{2,n}(1:2)==ptrloc)
                cmds{2,n}(1:2) = nan; % indicating no change
            else
              % move origin from lower-left corner to upper-left corner
              cmds{2,n}(2) = scrnsz(4)-cmds{2,n}(2);

              ptrloc = cmds{2,n}(2);  % update current mouse location
            end

        end

        % check time
        if isempty(deltaT)
            deltaT = zeros([1 dim(2:end)]);
        elseif any(numel(deltaT)~=1 || ~isnumeric(deltaT) || deltaT<0 || isnan(deltaT) || isinf(deltaT))
            msg = 'deltaT must be a non-negative scalar.';
        elseif numel(deltaT)==1
            deltaT = repmat(round(deltaT*1000),[1 dim(2:end)]);
        elseif numel(deltaT)~=numel(cmds)/2
            msg = 'deltaT must be a scalar or match the number of commands.';
        else
           deltaT = round(deltaT*1000);
        end

        end % parseinput()

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % expand all actions to include only 'XXX_down','XXX_up','move','wheel'
        function [actions,params,dt,msg] = convertToLowLevelActions(cmds,deltaT,idx_mouse,oemChars)

        import java.awt.event.KeyEvent

        actions = {};
        params = {};
        dt = [];
        msg = '';

        Ncmds = numel(deltaT);

        % determine total number of actions
        Nactions = zeros(Ncmds,1);
        for n = 1:Ncmds
            if idx_mouse(n) % mouse command
              % determine number of LL actions needed for each mouse command
                Nactions(n) = sum([...
                    cmds{2,n}(1)>=0 || cmds{2,n}(2)>=0% move
                    cmds{2,n}(3)~=0 % wheel
                    ~any(strcmpi(cmds{1,n},{'none','move','wheel'})) % click action
                    strcmpi(cmds{1,n},'normal')   % both down & up
                    3*any(strcmpi(cmds{1,n},{'open','extend','alternate'}))]); % down & up & modifier key down & up

                if strcmpi(cmds{1,n},'drag_on')
                    cmds{1,n}='left_down';
                elseif strcmpi(cmds{1,n},'drag_off')
                    cmds{1,n}='left_up';
                end
            else % keyboard command
                [cmds{1,n},cmds{2,n},msg] = parseKeyboardText(cmds{1,n},cmds{2,n},oemChars);
                if ~isempty(msg), return; end
                Nactions(n) = size(cmds{1,n},1);
            end
        end

        % create actions & params
        Ntotal = sum(Nactions);
        actions = cell(Ntotal,1);
        params = cell(Ntotal,1);
        dt = zeros(Ntotal,1);
        I = 1;
        for n = 1:Ncmds
            dt(I) = deltaT(n);
            if idx_mouse(n)
                if ~isnan(cmds{2,n}(1))
                    actions{I}='move';
                    params{I} = cmds{2,n}([1 2]);
                    I = I + 1;
                end
                switch lower(cmds{1,n})
                    case 'normal'
                        actions{I}='left_down';
                        I = I + 1;
                        actions{I}='left_up';
                        I = I + 1;
                    case 'extend'
                        actions{I}='key_down';
                        params{I} = KeyEvent.VK_SHIFT;
                        I = I+1;
                        actions{I}='left_down';
                        I = I + 1;
                        actions{I}='left_up';
                        I = I + 1;
                        actions{I}='key_up';
                        params{I} = KeyEvent.VK_SHIFT;
                        I = I+1;
                    case 'alternate'
                        actions{I}='key_down';
                        params{I} = KeyEvent.VK_CONTROL;
                        I = I+1;
                        actions{I}='left_down';
                        I = I + 1;
                        actions{I}='left_up';
                        I = I + 1;
                        actions{I}='key_up';
                        params{I} = KeyEvent.VK_CONTROL;
                        I = I+1;
                    case 'open'
                        actions{I}='left_down';
                        I = I + 1;
                        actions{I}='left_up';
                        I = I + 1;
                        actions{I}='left_down';
                        I = I + 1;
                        actions{I}='left_up';
                        I = I + 1;
                    case {'left_down','left_up',...
                            'middle_down','middle_up',...
                            'right_down','right_up'} % single button action
                        actions{I} = lower(cmds{1,n});
                        I = I + 1;
                end
                if cmds{2,n}(3)~=0
                    actions{I}='wheel';
                    params{I} = cmds{2,n}(3);
                 I = I + 1;
                end
            else % keyboard action
                idx = I:(I+Nactions(n)-1);
                actions(idx) = cmds{1,n};
                params(idx) = cmds{2,n};
                I = I + Nactions(n);
            end
        end

        end % convertToLowLevelActions

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % expand keyboard action to a sequence of 'key_down' and 'key_up'
        function	[actions,params,msg] = parseKeyboardText(cmd,text,oemChars)

        import java.awt.event.KeyEvent

        actions = {};
        params = {};
        msg = '';

        % virtual key codes for escape keys
        spKeys = [{'\\'      KeyEvent.VK_BACK_SLASH
                '\BACKSPACE'   KeyEvent.VK_BACK_SPACE
                '\TAB'         KeyEvent.VK_TAB
                '\ENTER'       KeyEvent.VK_ENTER
                '\SHIFT'       KeyEvent.VK_SHIFT
                '\CTRL'        KeyEvent.VK_CONTROL
                '\ALT'         KeyEvent.VK_ALT
                '\PAUSE'       KeyEvent.VK_PAUSE
                '\CAPSLOCK'    KeyEvent.VK_CAPS_LOCK
                '\ESC'         KeyEvent.VK_ESCAPE
                '\PAGEUP'      KeyEvent.VK_PAGE_UP
                '\PAGEDOWN'    KeyEvent.VK_PAGE_DOWN
                '\END'         KeyEvent.VK_END
                '\HOME'        KeyEvent.VK_HOME
                '\LEFT'        KeyEvent.VK_LEFT
                '\UP'          KeyEvent.VK_UP
                '\RIGHT'       KeyEvent.VK_RIGHT
                '\DOWN'        KeyEvent.VK_DOWN
              '\PRINTSCREEN' KeyEvent.VK_PRINTSCREEN
                '\INSERT'      KeyEvent.VK_INSERT
                '\DELETE'      KeyEvent.VK_DELETE
                '\WINDOWS'     KeyEvent.VK_WINDOWS
                '\NUMLOCK'     KeyEvent.VK_NUM_LOCK
                '\SCROLLLOCK'  KeyEvent.VK_SCROLL_LOCK}
              [strcat({'\F'},num2str((1:12)','%02d')) num2cell(KeyEvent.VK_F1+(0:11)')
              strcat({'\F'},num2str((13:24)','%02d')) num2cell(KeyEvent.VK_F13+(0:11)')]];

        % scan the text and determine regular vs. special key and total number of characters
        chCodes = zeros(numel(text),1);
        idx = 0;

        % ID type of each character in text: 
        %    regular char -> ASCII code
        %    escape char -> negative of row index to spKeys followed by zeros
        if numel(text)==1 % single character (to account for '\' case)
           chCodes = text;
        else % string of characters
           while ~isempty(text)
              if text(1)=='\' % first character = special character
                 tok = '';
                 rest = text;
              else % first character = regular character
                 [tok,rest] = strtok(text,'\');
              end
              chCodes(idx+(1:numel(tok))) = tok; % regular characters
              idx = idx+numel(tok);
              if isempty(rest), break; end

              I = find(cellfun(@(sch)strncmpi(rest,sch,numel(sch)),spKeys(:,1)));
              if isempty(I)
                 msg = 'Invalid key.';
                 return;
              end
              chCodes(idx+1)=-I;
              idx = idx + numel(spKeys{I});
              text = rest(numel(spKeys{I})+1:end);
           end
        end

        % remove the zero fillers for escape key string
        chCodes = chCodes(chCodes~=0);

        % total number of characters to type
        Nchars = numel(chCodes);

        % determine the virtual key codes for each key & need for shift key press
        keycodes = zeros(Nchars,1);
        shifts = false(Nchars,1);
        I = chCodes>0;

        %  - regular keys first
        vkc = char(chCodes(I));
        s = false(size(vkc));

        %  * upper case
        idx = vkc>='A' & vkc<='Z';
        s(idx) = true;

        %  * lower case
        idx = vkc>='a' & vkc<='z';
        vkc(idx) = vkc(idx)-'a'+'A';

        %  * OEM keys
        [idx1,J1] = ismember(vkc,oemChars(:,1)); % no shift
        [idx2,J2] = ismember(vkc,oemChars(:,2)); % with shift
        vkc(idx1) = cell2mat(oemChars(J1(idx1),3));
        vkc(idx2) = cell2mat(oemChars(J2(idx2),3));
        s(idx2) = true;

        %  - update regular keys
        keycodes(I) = vkc;
        shifts(I) = s;

        %  - then escape keys (no shift)
        I = ~I;
        keycodes(I) = cell2mat(spKeys(-chCodes(I),2));

        %determine how to press keys
        idx_type = strcmpi(cmd,{'key_normal','key_ctrl','key_alt','key_win'});
        mod = any(idx_type(2:4)); % true to press modifier key during
        press = ~strcmpi(cmd,'key_up'); % true to press key
        release = ~strcmpi(cmd,'key_down'); % true to release key
        shifts = [shifts(1);diff(shifts)]; % if 1, shift down, if -1 shift up (for each key)
        release_shift = sum(shifts)>0; % if last character needs shift down, release it at the end

        %determine total # of key down/up actions to perform
        Nacts = 2*mod ... % if needs modiefier key, must press & release (2)
           + Nchars*(press+release) ... % for each character, needs 2 (if down&up) or 1 (if down or up)
           + sum(shifts~=0) ... % for each shift press/release actions
            + release_shift; % to release the shift at the end

        % define the keyboard action command cell
        actions = cell(Nacts,1);
        params = cell(Nacts,1);
        I = 1; % action index
        if mod % press the modifier key
            actions{I}='key_down';
            switch find(idx_type,1)
                case 2 %'key_ctrl'
                    mod_param = KeyEvent.VK_CONTROL;
                case 3 %'key_alt'
                    mod_param = KeyEvent.VK_ALT;
                case 4 %'key_win'
                    mod_param = KeyEvent.VK_WINDOWS;
              otherwise % should never get here
                    mod_param = [];
            end
            params{I} = mod_param;
            I = I + 1;
        end
        for n = 1:Nchars; % for each character
            if shifts(n)~=0 % shift key action needed
                if shifts(n)>0 % 1:press shift key
                    actions{I}='key_down';
              else           % -1:release shift key
                    actions{I}='key_up';
                end
                params{I} = KeyEvent.VK_SHIFT;
                I = I + 1;
            end
            if press
                actions{I}='key_down';
                params{I} = keycodes(n);
                I = I + 1;
            end
            if release
                actions{I}='key_up';
                params{I} = keycodes(n);
                I = I+1;
            end
        end
        if release_shift
            actions{I}='key_up';
            params{I} = KeyEvent.VK_SHIFT;
            I = I + 1;
        end
        if mod % release the modifier key
            actions{I}='key_up';
            params{I} = mod_param;
        end

        end % parseKeyboardText

        % Copyright (c)2010, Takeshi Ikuma
        % All rights reserved.
        %
        % Redistribution and use in source and binary forms, with or without
        % modification, are permitted provided that the following conditions are
        % met:
        %
        %   * Redistributions of source code must retain the above copyright
        %   notice, this list of conditions and the following disclaimer. *
        %   Redistributions in binary form must reproduce the above copyright
        %   notice, this list of conditions and the following disclaimer in the
        %   documentation and/or other materials provided with the distribution.
        %   * Neither the names of its contributors may be used to endorse or
        %   promote products derived from this software without specific prior
        %   written permission.
        %
        % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
        % IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
        % THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
        % PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
        % CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
        % EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
        % PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
        % PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
        % LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
        % NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
        % SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    end
end
