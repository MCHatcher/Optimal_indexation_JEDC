%% ========================================================================
% JEDC_IT_baseline.m
%
% Author: Michael Hatcher
% Date:   Finalised on 16 Feb 2014 
%
% The author benefited from an original code for running loops in Dynare++
% written by Joris De Wind and downloaded from Wouter Den Haan's personal webpage. 
%The author would like to thank Joris De Wind and Wouter Den Haan for making the code available.

%==========================================================================

%% ========================================================================
% 0. Preliminaries
%==========================================================================

vcov = [3.100000e-003 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;... 
0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;... 
0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;... 
0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004 0;...
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.100000e-004];
%Variance-covariance matrix of shocks to be input into mod file

mydynareppdirectory = sprintf('C:\\dynare++'); 
%specify here the location of Dynare++

%tinit = 0.109219111; %perfect indexation
%tinit = 0.1092080; Correlated shocks and lagged indexation
%ntax = 50; %perfect indexation

tinit = 0.1134315;   %initial tax
%ntax = 160;         %dimension of loop over taxes
%nshare = 101;       %no. of indexation shares for which model simulated
stack_IT = [];       %vector containing mean govt spending for each tax/share combination
stacku_IT = [];      %vector containing mean utility for each tax/share combinatio
stacki_IT = [];      %vector containing indexation shares
tstack_IT = [];      %vector containing taxes that meet the govt. spending targe
stackmcy_IT = [];
stackmco_IT = [];
stackvcy_IT = [];
stackvco_IT = [];
irpdiffstack_IT = [];
stackri_IT = [];
stackrn_IT = [];
stackvri_IT = [];
stackvrn_IT = [];
stackmk_IT = [];
stackmy_IT = [];
WG_IT = [];               %welfare gain (or loss) in % aggregate consumption relative to zero indexation 
RA = 15;               %coefficient of relative risk aversion (must be changed to match mod file)

nshare = 101;
ntax = 200;

%% ========================================================================
% 1. Setting up inner and outer loops
%==========================================================================

%Indexation share loop (outer loop)

    for i = 1:nshare;
            sharei = (1*i-1)/100;
            stacki_IT(i) = sharei*100;
        
    %Tax loop (inner loop)
        
        for j = 1:ntax;
        
    if i==1    
        t = tinit + 0.0000012*j;
    elseif i>1    
        t = tend - 0.00005/(i^0.025) + 0.0000012*j;
        %sets tax at start of loop close to value that met govt spending target for previous indexation share (far more efficient) 
    end
    
      
%% ========================================================================
% 2. Preparing the olg_zin_index_IT_total.mod file
%==========================================================================

delete olg_zin_index_IT_total.mod                                      
%Otherwise the new stuff keeps on getting appended to old stuff.
!type olg_zin_index_IT_block_begin.mod >> olg_zin_index_IT_total.mod;
%this puts the *_begin.mod file into *_total.mod
temp = sprintf('echo t = %d; >> olg_zin_index_IT_total.mod',t);
system(temp);
temp = sprintf('echo sharei = %d; >> olg_zin_index_IT_total.mod',sharei);
system(temp);
%this creates a line that sets the value of sharei and adds it to *_total.mod
!type olg_zin_index_IT_block_end.mod >> olg_zin_index_IT_total.mod;
%this puts the *_end.mod file into *_total.mod
temp = sprintf('echo vcov = [%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d; %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d]; >> olg_zin_index_IT_total.mod',vcov);
%clc
%clear
system(temp);
%this creates a line that sets the value of vcov and adds it to *_total.mod

%% ========================================================================
% 3. Running Dynare++
%==========================================================================

temp = sprintf('%s\\dynare++ --order 2 --sim 100 --per 1100 --burn 100 --no-irfs olg_zin_index_IT_total.mod',...
    mydynareppdirectory);
system(temp);
load 'olg_zin_index_IT_total.mat'

stack_IT(i,j) = dyn_mean(15);
%unconditional mean of govt spending ratio

        if stack_IT(i,j) > 0.11000004
        stack_IT(i,j) = 0;
        break    
    end    

stacku_IT(i,j) = dyn_mean(14);
%unconditional mean of lifetime utility (=Social Welfare)

    if stack_IT(i,j) > 0.11000004
        stacku_IT(i,j) = 0;
        break
    end    

delete olg_zin_index_IT_total.jnl
delete olg_zin_index_IT_total_f.m
delete olg_zin_index_IT_total_ff.m
delete olg_zin_index_IT_total.dump

        end

        %Indexation and real variables
        tend = t;   %last tax in loop (satisfies govt spending target)
        tstack_IT = [tstack_IT; tend];   %vector of taxes that satisfy govt spending target (one for each indexation share)
        stackmcy_IT = [stackmcy_IT; dyn_mean(3)];   %Corresponding mean consumption by the young
        stackmco_IT = [stackmco_IT; dyn_mean(29)];  %Corresponding mean consumption by the old
        stackvcy_IT = [stackvcy_IT; dyn_vcov(3,3)]; %Variance of consumption by the young 
        stackvco_IT = [stackvco_IT; dyn_vcov(29,29)];   %Variance of consumption by the old
        irpdiffstack_IT = [irpdiffstack_IT; dyn_mean(16)];  %Average difference between inflation risk premium on nominal and indexed bonds
        stackri_IT = [stackri_IT; dyn_mean(17)];    %Mean return on indexed bonds
        stackrn_IT = [stackrn_IT; dyn_mean(7)];     %Mean return on nominal bonds
        stackvri_IT = [stackvri_IT; dyn_vcov(17,17)];   %Variance of return on indexed bonds
        stackvrn_IT = [stackvrn_IT; dyn_vcov(7,7)];     %Variance of return on nominal bonds
        stackmk_IT = [stackmk_IT; dyn_mean(19)];    %Mean capital
        stackmy_IT = [stackmy_IT; dyn_mean(9)];    %Mean output
    end

g = stack_IT';
gmean_IT = max(g);
%Picks out mean govt spending with last tax in the loop
%This line uses the fact that mean(g) rises as the tax rate rises (this was
%verified for each simulation)
    
u = stacku_IT';    
utility_IT = min(u);
%Picks out mean utility (= Social Welfare) for the last tax in the loop
%This line uses the fact that mean(utility) falls as the tax rate rises (this was
%verified for each simulation

%% ========================================================================
% 4. Calculating the welfare gain in consimption term
%==========================================================================

for i=1:nshare

WG_IT(i) = 100*( (utility_IT(i)/utility_IT(1))^(1/(1-RA)) - 1 );
%welfare gain relative to the case of zero indexation

end


%% ========================================================================
% 5. Plotting the results
%==========================================================================
figure(1)
    plot(stacki_IT, WG_IT);
    title('Utility','fontsize',10)
    
figure (2)
    plot(stacki_IT, gmean_IT);
    title('Mean g','fontsize',10)
       
figure(3)
    hold on,
   
    subplot(3,4,1); plot(stacki_IT, stackmcy_IT);
    title('Mean C by young','fontsize',10)
    
    subplot(3,4,2); plot(stacki_IT, stackmco_IT);
    title('Mean C by old','fontsize',10)
    
    subplot(3,4,3); plot(stacki_IT, stackvcy_IT);
    title('Var C by young','fontsize',10)
    
    subplot(3,4,4); plot(stacki_IT, stackvco_IT);
    title('Var C by old','fontsize',10)
    
    subplot(3,4,5); plot(stacki_IT, stackri_IT);
    title('Mean return indexed','fontsize',10)
    
    subplot(3,4,6); plot(stacki_IT, stackrn_IT);
    title('Mean return nominal','fontsize',10)
    
    subplot(3,4,7); plot(stacki_IT, stackvri_IT);
    title('Var return indexed','fontsize',10)
    
    subplot(3,4,8); plot(stacki_IT, stackvrn_IT);
    title('Var return nominal','fontsize',10)
    
    subplot(3,4,9); plot(stacki_IT, irpdiffstack_IT);
    title('IRP diff','fontsize',10)
    
    subplot(3,4,10); plot(stacki_IT, tstack_IT);
    title('Taxes','fontsize',10)
    
    subplot(3,4,11); plot(stacki_IT, stackmk_IT);
    title('Capital','fontsize',10)
    
    subplot(3,4,12); plot(stacki_IT, stackmy_IT);
    title('Output','fontsize',10)

    hold off