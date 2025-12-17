clear; 
close all;

Nrun = 1;


st = 1;
for M = st:Nrun
    % close all;
    clearvars -except Nrun M 
    % close all

    rng('shuffle');
     
    
    % time settings
    tN = 25000; % total timesteps
    dt = 0.01;  % time step


    %--------------------
    % Excitable system parameters
    %--------------------
    
    
    q.A = 0.0167;
    q.Kg = 0.9393;
    q.Kr = 0.8125;
    q.a = 0.1863;
    finaltau = 100;   % influences dr update

    
    %========================
    % Small simulation of single-cell excitable dynamics
    %========================
   
    
    g = 0;
    r = 0;
    s = 0;
    recordg = [];
    recordr = [];
    records = [];

    for i = 2:tN
        % Set an imposed signal s for testing early transient
        if i == 2
            s = 100;
        elseif i >= 3000
            s = 0;
        end

        % Trigger
        I = q.A*log10(1 + s);
        

        % Excitable dynamics 
        dg = (-q.Kg*g.*(g-1).*(g-q.a) - q.Kr*r + I);
        
        dr = ((0.5.*g-r)./finaltau);

        rnew = r + dr*dt;
        gnew = g + dg*dt;
        g = gnew;
        r = rnew;

        recordg(end+1) = g;
        recordr(end+1) = r;
    end
    

end

% Plot results of single-cell excitable dynamics
figure;
plot(recordg)
hold on
plot(recordr)

hilbertresult = hilbert(recordg);
phase = angle(hilbertresult);

figure
hold on
plot(phase, recordg)
yyaxis right 
plot(phase, recordr)

% Save hilbert result for later interpolation use
hil.phase = phase;
hil.g = recordg;
hil.r = recordr;
save("hilresult", "hil")
