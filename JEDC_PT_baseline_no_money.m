%% ========================================================================
% JEDC_PT_baseline_no_money.m
%
% Author: Michael Hatcher
% Date:   Finalised on 9 Feb 2014 
%
% The author benefited from an original code for running loops in Dynare++
% written by Joris De Wind and downloaded from Wouter Den Haan's personal webpage. 
% The author would like to thank Joris De Wind and Wouter Den Haan for making the code available.

%==========================================================================

%% ========================================================================
% 0. Preliminaries
%==========================================================================

vcov = [3.10000e-003 0  0;...
    0  0.00011  0.00000;...
    0  0.00000  0.00011];
%Variance-covariance matrix of shocks to be input into mod file

mydynareppdirectory = sprintf('C:\\dynare++'); 
%specify here the location of Dynare++

tinit = 0.12632250;   %initial tax
%ntax = 60;             %dimension of loop over taxes
%nshare = 101;          %no. of indexation shares for which model simulated
stack_no_money_PT = [];            %vector containing mean govt spending for each tax/share combination
stacku_no_money_PT = [];           %vector containing mean utility for each tax/share combination
stacki_no_money_PT = [];           %vector containing indexation shares
tstack_no_money_PT = [];           %vector containing taxes that meet the govt. spending target
stackmcy_no_money_PT = [];
stackmco_no_money_PT = [];
stackvcy_no_money_PT = [];
stackvco_no_money_PT = [];
irpdiffstack_no_money_PT = [];
stackri_no_money_PT = [];
stackrn_no_money_PT = [];
stackvri_no_money_PT = [];
stackvrn_no_money_PT = [];
stackmk_no_money_PT = [];
stackmy_no_money_PT = [];
WG_no_money_PT = [];               %welfare gain (or loss) in % aggregate consumption relative to zero indexation 
RA = 15;               %coefficient of relative risk aversion (must be changed to match mod file)

ntax = 300;
nshare = 101;


%% ========================================================================
% 1. Setting up inner and outer loops
%==========================================================================

%Indexation share loop (outer loop)
    for i = 1:nshare;
            sharei = (1*i-1)/100;
            stacki_no_money_PT(i) = sharei*100;
            
   %Tax loop (inner loop)         
        for j = 1:ntax;
        
    if i==1    
        t = tinit + 0.00000011*j;  %smaller step size for accuracy
    elseif i>1    
        t = tend - 0.0000005/(i^0.025) + 0.00000011*j;    
        %sets tax at start of loop close to value that met govt spending target for previous indexation share (far more efficient) 
    end
    
      
%% ========================================================================
% 2. Preparing the olg_zin_index_PT_no_money_total.mod file
%==========================================================================

delete olg_zin_index_PT_no_money_total.mod                                      
%Otherwise the new stuff keeps on getting appended to old stuff.
!type olg_zin_index_PT_no_money_block_begin.mod >> olg_zin_index_PT_no_money_total.mod;
%this puts the *_begin.mod file into *_total.mod
temp = sprintf('echo t = %d; >> olg_zin_index_PT_no_money_total.mod',t);
system(temp);
temp = sprintf('echo sharei = %d; >> olg_zin_index_PT_no_money_total.mod',sharei);
system(temp);
%this creates a line that sets the value of sharei and adds it to *_total.mod
!type olg_zin_index_PT_no_money_block_end.mod >> olg_zin_index_PT_no_money_total.mod;
%this puts the *_end.mod file into *_total.mod
temp = sprintf('echo vcov = [%d %d %d; %d %d %d; %d %d %d]; >> olg_zin_index_PT_no_money_total.mod',vcov);
system(temp);
%this creates a line that sets the value of vcov and adds it to *_total.mod

%% ========================================================================
% 3. Running Dynare++
%==========================================================================

temp = sprintf('%s\\dynare++ --order 2 --sim 100 --per 1100 --burn 100 --no-irfs olg_zin_index_PT_no_money_total.mod',...
    mydynareppdirectory);
system(temp);
load 'olg_zin_index_PT_no_money_total.mat'

stack_no_money_PT(i,j) = dyn_mean(14);
%unconditional mean of govt spending ratio
   
     if stack_no_money_PT(i,j) > 0.11000004
        stack_no_money_PT(i,j) = 0;
        break    
     end    

stacku_no_money_PT(i,j) = dyn_mean(13);
%unconditional mean of lifetime utility (=Social Welfare)

    if stack_no_money_PT(i,j) > 0.11000004
        stacku_no_money_PT(i,j) = 0;
        break
    end    

delete olg_zin_index_PT_no_money_total.jnl
delete olg_zin_index_PT_no_money_total_f.m
delete olg_zin_index_PT_no_money_total_ff.m
delete olg_zin_index_PT_no_money_total.dump

        end

        %Indexation and real variables
        tend = t;   %last tax in loop (satisfies govt spending target)
        tstack_no_money_PT = [tstack_no_money_PT; tend];   %vector of taxes that satisfy govt spending target (one for each indexation share)
        stackmcy_no_money_PT = [stackmcy_no_money_PT; dyn_mean(3)];   %Corresponding mean consumption by the young
        stackmco_no_money_PT = [stackmco_no_money_PT; dyn_mean(27)];  %Corresponding mean consumption by the old
        stackvcy_no_money_PT = [stackvcy_no_money_PT; dyn_vcov(3,3)]; %Variance of consumption by the young 
        stackvco_no_money_PT = [stackvco_no_money_PT; dyn_vcov(27,27)];   %Variance of consumption by the old
        irpdiffstack_no_money_PT = [irpdiffstack_no_money_PT; dyn_mean(16)];  %Average difference between inflation risk premium on nominal and indexed bonds
        stackri_no_money_PT = [stackri_no_money_PT; dyn_mean(15)];    %Mean return on indexed bonds
        stackrn_no_money_PT = [stackrn_no_money_PT; dyn_mean(7)];     %Mean return on nominal bonds
        stackvri_no_money_PT = [stackvri_no_money_PT; dyn_vcov(15,15)];   %Variance of return on indexed bonds
        stackvrn_no_money_PT = [stackvrn_no_money_PT; dyn_vcov(7,7)];     %Variance of return on nominal bonds
        stackmk_no_money_PT = [stackmk_no_money_PT; dyn_mean(18)];    %Mean capital
        stackmy_no_money_PT = [stackmy_no_money_PT; dyn_mean(8)];    %Mean output

    end

g = stack_no_money_PT';
gmean_no_money_PT = max(g);
%Picks out mean govt spending with last tax in the loop
%This line uses the fact that mean(g) rises as the tax rate rises (this was
%verified for each simulation)
    
u = stacku_no_money_PT';    
utility_no_money_PT = min(u);
%Picks out mean utility (= Social Welfare) for the last tax in the loop
%This line uses the fact that mean(utility) falls as the tax rate rises (this was
%verified for each simulation) 


%% ========================================================================
% 4. Calculating the welfare gain in consimption term
%==========================================================================

for i=1:nshare

WG_no_money_PT(i) = 100*( (utility_no_money_PT(i)/utility_no_money_PT(1))^(1/(1-RA)) - 1 );
%welfare gain relative to the case of zero indexation

end

%% ========================================================================
% 5. Plotting the results
%==========================================================================
  

figure(1)
    plot(stacki_no_money_PT, WG_no_money_PT);
    title('Utility','fontsize',10)
    
figure (2)
    plot(stacki_no_money_PT, gmean_no_money_PT);
    title('Mean g','fontsize',10)
       
figure(3)
    hold on,
   
    subplot(3,4,1); plot(stacki_no_money_PT, stackmcy_no_money_PT);
    title('Mean C by young','fontsize',10)
    
    subplot(3,4,2); plot(stacki_no_money_PT, stackmco_no_money_PT);
    title('Mean C by old','fontsize',10)
    
    subplot(3,4,3); plot(stacki_no_money_PT, stackvcy_no_money_PT);
    title('Var C by young','fontsize',10)
    
    subplot(3,4,4); plot(stacki_no_money_PT, stackvco_no_money_PT);
    title('Var C by old','fontsize',10)
    
    subplot(3,4,5); plot(stacki_no_money_PT, stackri_no_money_PT);
    title('Mean return indexed','fontsize',10)
    
    subplot(3,4,6); plot(stacki_no_money_PT, stackrn_no_money_PT);
    title('Mean return nominal','fontsize',10)
    
    subplot(3,4,7); plot(stacki_no_money_PT, stackvri_no_money_PT);
    title('Var return indexed','fontsize',10)
    
    subplot(3,4,8); plot(stacki_no_money_PT, stackvrn_no_money_PT);
    title('Var return nominal','fontsize',10)
    
    subplot(3,4,9); plot(stacki_no_money_PT, irpdiffstack_no_money_PT);
    title('IRP diff','fontsize',10)
    
    subplot(3,4,10); plot(stacki_no_money_PT, tstack_no_money_PT);
    title('Taxes','fontsize',10)
    
    subplot(3,4,11); plot(stacki_no_money_PT, stackmk_no_money_PT);
    title('Capital','fontsize',10)
    
    subplot(3,4,12); plot(stacki_no_money_PT, stackmy_no_money_PT);
    title('Output','fontsize',10)

    hold off