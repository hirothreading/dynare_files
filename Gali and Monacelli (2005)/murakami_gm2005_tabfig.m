% Generate Table 1 of Gali & Monacelli (2005) and plots for different
% regimes

if exist('murakami_gm2005_OP_1st.mat') && exist('murakami_gm2005_DITR_1st.mat') && exist('murakami_gm2005_CITR_1st.mat') && exist('murakami_gm2005_PEG_1st.mat')
    OP1=load('murakami_gm2005_OP_1st.mat');
    DITR1=load('murakami_gm2005_DITR_1st.mat');
    CITR1=load('murakami_gm2005_CITR_1st.mat');
    PEG1=load('murakami_gm2005_PEG_1st.mat');
    var_string={'pih','x','pic','tot','r','der'};
    
    %find output gap and inflation in covariance matrix
    x_pos=strmatch('x',OP1.M_.endo_names,'exact'); %can use across regimes
    pih_pos=strmatch('pih',OP1.M_.endo_names,'exact');
    
    %%% Cost push shock plot
    figure
    for fig_iter=1:length(var_string)
        subplot(3,2,fig_iter)
        plot(1:20,OP1.oo_.irfs.([var_string{fig_iter},'_epsu']),'r-',1:20,DITR1.oo_.irfs.([var_string{fig_iter},'_epsu']),'b--',1:20,CITR1.oo_.irfs.([var_string{fig_iter},'_epsu']),'g-x',1:20,PEG1.oo_.irfs.([var_string{fig_iter},'_epsu']),'m-s');
        %Need to annualise inflation and interest rates 
        %also add legend to first subplot
        if fig_iter==1 || fig_iter==3 || fig_iter==5
            % need to annualise rates
            plot(1:20,4*OP1.oo_.irfs.([var_string{fig_iter},'_epsu']),'r-',1:20,4*DITR1.oo_.irfs.([var_string{fig_iter},'_epsu']),'b--',1:20,4*CITR1.oo_.irfs.([var_string{fig_iter},'_epsu']),'g-x',1:20,4*PEG1.oo_.irfs.([var_string{fig_iter},'_epsu']),'m-s');
        end
        grid on
        title(OP1.M_.endo_names_long(strmatch(var_string{fig_iter},OP1.M_.endo_names,'exact'),:))
    end
    hlegend = legend('Optimal Policy','Domestic-Taylor','CPI-Taylor','FOREX Peg', 'Orientation', 'horizontal', 'Location', 'southoutside');
    hlegend.Position(1)=0.1775; %centre at bottom
    hlegend.Position(2)=0.009; %centre at bottom
    set(gcf,'Color','w'); %set background colour
    export_fig murakami_gm2005_costpushshock.pdf
    
    
    %%% Natural rate of output shock 
    figure
    for fig_iter=1:length(var_string)
        subplot(3,2,fig_iter)
        plot(1:20,OP1.oo_.irfs.([var_string{fig_iter},'_epsz']),'r-',1:20,DITR1.oo_.irfs.([var_string{fig_iter},'_epsz']),'b--',1:20,CITR1.oo_.irfs.([var_string{fig_iter},'_epsz']),'g-x',1:20,PEG1.oo_.irfs.([var_string{fig_iter},'_epsz']),'m-s');
        %Need to annualise inflation and interest rates 
        %also add legend to first subplot
        if fig_iter==1 || fig_iter==3 || fig_iter==5
            plot(1:20,4*OP1.oo_.irfs.([var_string{fig_iter},'_epsz']),'r-',1:20,4*DITR1.oo_.irfs.([var_string{fig_iter},'_epsz']),'b--',1:20,4*CITR1.oo_.irfs.([var_string{fig_iter},'_epsz']),'g-x',1:20,4*PEG1.oo_.irfs.([var_string{fig_iter},'_epsz']),'m-s');
        end
        grid on
        title(OP1.M_.endo_names_long(strmatch(var_string{fig_iter},OP1.M_.endo_names,'exact'),:))
    end
    hlegend = legend('Optimal Policy','Domestic-Taylor','CPI-Taylor','FOREX Peg', 'Orientation', 'horizontal', 'Location', 'southoutside');
    hlegend.Position(1)=0.1775; 
    hlegend.Position(2)=0.009;
    set(gcf,'Color','w');
    export_fig murakami_gm2005_ynshock.pdf
    
    
    %%% Natural interest rate shock 
    figure
    for fig_iter=1:length(var_string)
        subplot(3,2,fig_iter)
        plot(1:20,OP1.oo_.irfs.([var_string{fig_iter},'_epsrn']),'r-',1:20,DITR1.oo_.irfs.([var_string{fig_iter},'_epsrn']),'b--',1:20,CITR1.oo_.irfs.([var_string{fig_iter},'_epsrn']),'g-x',1:20,PEG1.oo_.irfs.([var_string{fig_iter},'_epsrn']),'m-s');
        %Need to annualise inflation and interest rates 
        %also add legend to first subplot
        if fig_iter==1 || fig_iter==3 || fig_iter==5
            plot(1:20,4*OP1.oo_.irfs.([var_string{fig_iter},'_epsu']),'r-',1:20,4*DITR1.oo_.irfs.([var_string{fig_iter},'_epsrn']),'b--',1:20,4*CITR1.oo_.irfs.([var_string{fig_iter},'_epsrn']),'g-x',1:20,4*PEG1.oo_.irfs.([var_string{fig_iter},'_epsrn']),'m-s');
        end
        grid on
        title(OP1.M_.endo_names_long(strmatch(var_string{fig_iter},OP1.M_.endo_names,'exact'),:))
    end
    hlegend = legend('Optimal Policy','Domestic-Taylor','CPI-Taylor','FOREX Peg', 'Orientation', 'horizontal', 'Location', 'southoutside');
    hlegend.Position(1)=0.1775; 
    hlegend.Position(2)=0.009;
    set(gcf,'Color','w');
    export_fig murakami_gm2005_rnshock.pdf
    
    
    %%% Foreign interest rate shock
    figure
    for fig_iter=1:length(var_string)
        subplot(3,2,fig_iter)
        plot(1:20,OP1.oo_.irfs.([var_string{fig_iter},'_epsrstar']),'r-',1:20,DITR1.oo_.irfs.([var_string{fig_iter},'_epsrstar']),'b--',1:20,CITR1.oo_.irfs.([var_string{fig_iter},'_epsrstar']),'g-x',1:20,PEG1.oo_.irfs.([var_string{fig_iter},'_epsrstar']),'m-s');
        %Need to annualise inflation and interest rates 
        %also add legend to first subplot
        if fig_iter==1 || fig_iter==3 || fig_iter==5
            plot(1:20,4*OP1.oo_.irfs.([var_string{fig_iter},'_epsrstar']),'r-',1:20,4*DITR1.oo_.irfs.([var_string{fig_iter},'_epsrstar']),'b--',1:20,4*CITR1.oo_.irfs.([var_string{fig_iter},'_epsrstar']),'g-x',1:20,4*PEG1.oo_.irfs.([var_string{fig_iter},'_epsrstar']),'m-s');
        end
        grid on
        title(OP1.M_.endo_names_long(strmatch(var_string{fig_iter},OP1.M_.endo_names,'exact'),:))
    end
    hlegend = legend('Optimal Policy','Domestic-Taylor','CPI-Taylor','FOREX Peg', 'Orientation', 'horizontal', 'Location', 'southoutside');
    hlegend.Position(1)=0.1775; 
    hlegend.Position(2)=0.009;
    set(gcf,'Color','w');
    export_fig murakami_gm2005_rstarshock.pdf 
   
   
    fprintf('\nTABLE 1: Cyclical Properties of Alternative Policy Regimes\n')
    fprintf('%20s \t %3s \t %3s \t %3s \t %3s \n','','OP','DITR','CITR','PEG')
    fprintf('%20s \t %3s \t %3s \t %3s \t %3s \n','','sd%','sd%','sd%','sd%')
for var_iter=1:length(var_string)
        fprintf('%20s \t %3.2f \t %3.2f \t %3.2f \t %3.2f \n',OP1.M_.endo_names_long{strmatch(var_string{var_iter},OP1.M_.endo_names,'exact')},OP1.cyc_moments(var_iter),DITR1.cyc_moments(var_iter),CITR1.cyc_moments(var_iter),PEG1.cyc_moments(var_iter))
end


    VOP
    VDITR
    VCITR
    VPEG

end   