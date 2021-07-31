%JEDC IT-PT plotter

figure(1)
    hold on,
    plot(stacki_IT, WG_IT);
    plot(stacki_IT, WG_PT);
    title('Utility','fontsize',10)
    hold off
    
figure (2)
    hold on,
    plot(stacki_IT, gmean_IT);
    plot(stacki_IT, gmean_PT);
    title('Mean g','fontsize',10)
    hold off   
    
figure(3)
    subplot(3,4,1); hold on, plot(stacki_IT, stackmcy_IT),
    plot(stacki_IT, stackmcy_PT);
    title('Mean C by young','fontsize',10)
    hold off
    
    subplot(3,4,2); hold on, plot(stacki_IT, stackmco_IT);
    plot(stacki_IT, stackmco_PT);
    title('Mean C by old','fontsize',10)
    hold off
    
    subplot(3,4,3); hold on, plot(stacki_IT, stackvcy_IT);
    plot(stacki_IT, stackvcy_PT);
    title('Var C by young','fontsize',10)
    hold off
    
    subplot(3,4,4); hold on, plot(stacki_IT, stackvco_IT);
    plot(stacki_IT, stackvco_PT)
    title('Var C by old','fontsize',10)
    hold off
    
    subplot(3,4,5); hold on, plot(stacki_IT, stackri_IT);
    plot(stacki_IT, stackri_PT);
    title('Mean return indexed','fontsize',10)
    hold off
    
    subplot(3,4,6); hold on, plot(stacki_IT, stackrn_IT);
    plot(stacki_IT, stackrn_PT);
    title('Mean return nominal','fontsize',10)
    hold off
    
    subplot(3,4,7); hold on, plot(stacki_IT, stackvri_IT);
    plot(stacki_IT, stackvri_PT);
    title('Var return indexed','fontsize',10)
    hold off
    
    subplot(3,4,8); hold on, plot(stacki_IT, stackvrn_IT);
    plot(stacki_IT, stackvrn_PT);
    title('Var return nominal','fontsize',10)
    hold off
    
    subplot(3,4,9); hold on, plot(stacki_IT, irpdiffstack_IT);
    plot(stacki_IT, irpdiffstack_PT);
    title('IRP diff','fontsize',10)
    hold off
   
    subplot(3,4,10); hold on, plot(stacki_IT, tstack_IT);
    plot(stacki_IT, tstack_PT);
    title('Taxes','fontsize',10)
    hold off
    
    subplot(3,4,11); hold on, plot(stacki_IT, stackmk_IT);
    plot(stacki_IT, stackmk_PT);
    title('Capital','fontsize',10)
    hold off
    
    subplot(3,4,12); hold on, plot(stacki_IT, stackmy_IT);
    plot(stacki_IT, stackmy_PT);
    title('Output','fontsize',10)
    hold off