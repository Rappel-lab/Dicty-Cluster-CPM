close all;
clear; 


wavespeedall = [];
wavespeedvarall = [];

num = 200;
wavespeed = [];
for ens = 1:1

    foldername = ['Viscous_'  num2str(num)];

    %% ----------- Uncomment when you run multiple ensumble and named in this way ---------------
    % entire = load([foldername '/entire_' num2str(5) 'ens_' num2str(ens) '.mat'], "entire").entire;
    
    entire = load([foldername '/entire_' num2str(num/20) '.mat' ], "entire").entire;
    if max(entire.S, [], "all") > 1e-5
        wavespeed(end+1) = 60*2*pi/mean(entire.corrperiod_diffr((entire.reff>0.6) & (entire.reff<0.8)) );
    end   
end

wavespeedall(end+1) = mean(wavespeed);
wavespeedvarall(end+1) = std(wavespeed);





dt = 0.01;


%% ----------- Three methods for measuring angular velocity of wave ---------------
%% Method (1): Meausre at fixed range of normalized distance.  
%% Method (2): Meausre at fixed range of real distance.
%% Method (3): Meausre at fixed maximum tangential velocity.

method1omega = zeros(5, 1);
method2omega = zeros(5, 1);
method3omega = zeros(5, 1);

method1omegastd = zeros(5, 1);
method2omegastd = zeros(5, 1);
method3omegastd = zeros(5, 1);

plotallcell = figure;

for num = 200:50:200

    %% ----------- Uncomment when you run multiple ensumble and named in this way ---------------
    % conf = load(['Viscous_' num2str(num) '/entire_' num2str(num/20) 'ens_1.mat']).entire;
    conf = load(['Viscous_' num2str(num) '/entire_' num2str(num/20) '.mat']).entire;
    
    %%
    
    
    p.ncell = num;
    x_umpergrid = 1.25;     % physical conversion,  micron per grid
    t_sperMCS = 1.55;          % physical conversion,  sec per MCS
    p.area = round((pi*5^2/x_umpergrid^2)); % cell area
    hexa_packing = 0.91;
    Aeff = num*p.area/hexa_packing;
    
    Reff = sqrt(Aeff/pi);
    rfirst = 2*sqrt(p.area./pi);
    
    
    
    dr = round(2*sqrt(p.area./pi)/2);
    rmin = rfirst;
    rmax = sqrt(p.area*p.ncell/pi)*1.5;
    N = ceil((rmax - rmin)/dr);
    ravg = zeros(1,N);
    allvel = zeros(N, 10);

    for ens = 1:10

        %% ----------- Uncomment when you run multiple ensumble and named in this way ---------------
        % conf = load(['Viscous_' num2str(num) '/entire_' num2str(num/20) 'ens_' num2str(ens) '.mat']).entire;
        conf = load(['Viscous_' num2str(num) '/entire_' num2str(num/20)  '.mat']).entire;
        
        velstep = 5;
        cmx_all = conf.cmx;
        cmy_all = conf.cmy;
        
        allstep = size(cmx_all, 1);
        startstep = allstep/2;
        measureN = round((allstep-startstep)/velstep);
        allvel_T = zeros(measureN, N);
        rclus = sqrt(p.area*p.ncell/pi);
        
        ravg = zeros(1,N);
        measureN = round((allstep-startstep)/velstep)-1;
        allvel_T = zeros(measureN, N);

        
        
        for measureindex = 1:measureN
        
            cmx = cmx_all(startstep+(measureindex-1)*velstep,:);
            cmy = cmy_all(startstep+(measureindex-1)*velstep,:);
            nextcmx = cmx_all(startstep+measureindex*velstep,:);
            nextcmy = cmy_all(startstep+measureindex*velstep,:);
            
            vx = nextcmx - cmx;
            vy = nextcmy - cmy;
            [theta, rho] = cart2pol((cmx - mean(cmx)),( cmy - mean(cmy)));
            [Vtheta,Vrho] = cart2pol(vx,vy);
            Vradial = Vrho.*cos(Vtheta - theta) * x_umpergrid /(dt*200*velstep) * 60;     % um/min
            Vtangent = Vrho.*sin(Vtheta - theta) * x_umpergrid /(dt*200*velstep) * 60;    % um/min
            velR_dr_avg = zeros(1,N);
            velT_dr_avg = zeros(1,N);
            for i = 1:N
               rlb = rmin + dr*(i-1);
               rub = rmin + dr*(i);
               ravg(i) =(rlb+rub)/2;
            
               valid = (rho >= rlb) & (rho <= rub);
            
               valid = single(valid);
               valid(valid == 0) = NaN; % exclude unmasked data
               
               velR_dr_avg(:,i) = mean(Vradial.*valid, 2, 'omitnan');
               velT_dr_avg(:,i) = mean(Vtangent.*valid, 2, 'omitnan');
              
            
            end
            allvel_T(measureindex,:) = velT_dr_avg;
        end
        mostnan = sum(isnan(allvel_T), 1);
        

        effvTr = mean(allvel_T, 1, 'omitnan');
        effvTr(mostnan>size(allvel, 1)*0.5) = nan;
        effr = ravg;
        
        vTmaxindex = find(effvTr == max(effvTr));
        allvel(:,ens) = effvTr;
    end

    normravg = ravg./sqrt(p.area*p.ncell/pi);
    allvT = mean(allvel, 2)';
    allomega = mean(allvel, 2)'./(ravg.*x_umpergrid);
    stdomega = std(allvel, 1, 2)'./(ravg.*x_umpergrid);
    if mod(num, 100) == 0
        figure(plotallcell)
        hold on 
        errorbar(ravg./sqrt(p.area*p.ncell/pi),allomega, stdomega, DisplayName=['$N=$' num2str(num) ', $\beta=5(s^{-1})$' ], LineWidth=2)
        
        hold off
    end
    method1omega(num/50-1) = mean(allomega(normravg>0.7&normravg<0.9), 'all');
    method2omega(num/50-1) = mean(allomega(ravg>25&ravg<40), 'all', 'omitnan');
    method3omega(num/50-1) = mean(allomega(allvT == max(allvT, [],'all')), 'all', 'omitnan');
    

    method1omegastd(num/50-1) = mean(stdomega(normravg>0.7&normravg<0.9), 'all');
    method2omegastd(num/50-1) = mean(stdomega(ravg>25&ravg<40), 'all', 'omitnan');
    method3omegastd(num/50-1) = mean(stdomega(allvT == max(allvT, [],'all')), 'all', 'omitnan');
    
    
    

end





legendFont = 15;
axFont = 15;
labelFont = 30;
cellspeedall = method3omega;
cellspeedvarall = method3omegastd;


area = [100, 150, 200, 250, 300];
area = area.*p.area*(x_umpergrid^2);

figure(plotallcell)
hold on
box on
ax = gca;   

axes(ax);
ax.LineWidth = 2;
ax.FontSize = axFont;
xlabel('$r/R_{max}$ ', 'Interpreter', 'latex','fontsize',labelFont)
ylabel('$\omega_{cell}$(rad/min)', 'Interpreter', 'latex','fontsize',labelFont)
legend('Interpreter','latex',FontSize=legendFont)

hold off


%% ----------- Uncomment when you run multiple ensumble  ---------------

%{

figure
hold on
errorbar(area, wavespeedall, wavespeedvarall, 'b', LineWidth=2)
box on
ax = gca;   
axes(ax);
ax.LineWidth = 2;
ax.FontSize = axFont;
xlabel('Area $\times 10^4$  ($\mu m^2$)', 'Interpreter', 'latex','fontsize',labelFont)
ylabel('$\omega_{wave}$(rad/min)', 'Interpreter', 'latex','fontsize',labelFont)
xticks([5000, 10000, 15000, 20000, 25000])
xticklabels([0.5, 1, 1.5, 2, 2.5])
ylim([1.1, 1.3])
saveas(gca, './figure/wavesize','svg')


figure

hold on 
errorbar(area, cellspeedall, cellspeedvarall, 'b', LineWidth=2)
box on
ax = gca;   
axes(ax);
ax.LineWidth = 2;
ax.FontSize = axFont;
xlabel('Area $\times 10^4$  ($\mu m^2$)', 'Interpreter', 'latex','fontsize',labelFont)
ylabel('$\omega_{cell}$ (rad/min)', 'Interpreter', 'latex','fontsize',labelFont)
xticks([5000, 10000, 15000, 20000, 25000])
xticklabels([0.5, 1, 1.5, 2, 2.5])
hold off
saveas(gca, './figure/cellsize','svg')


figure

hold on
errorbar(cellspeedall,wavespeedall, wavespeedvarall, Color='b', LineWidth=2, LineStyle='none',HandleVisibility='off')
errorbar(cellspeedall,wavespeedall, cellspeedvarall, 'horizontal', Color='b', LineWidth=2, LineStyle='none',HandleVisibility='off')

[P, S] = polyfit(cellspeedall,wavespeedall,1);
yfit = P(1)*cellspeedall+P(2);
plot(cellspeedall, yfit, 'k--', LineWidth=2, DisplayName=['$R^2$=' num2str(S.rsquared)])
ylim([1, 1.35])
legend(Location="northwest", Interpreter="latex", FontSize=legendFont)
box on
ax = gca;   
axes(ax);
ax.LineWidth = 2;
ax.FontSize = axFont;
xlabel('$\omega_{cell}$ (rad/min)', 'Interpreter', 'latex','fontsize', labelFont)
ylabel('$\omega_{wave}$(rad/min)', 'Interpreter', 'latex','fontsize',labelFont)

saveas(gca,'./figure/both','svg')



%}

