clear; close all;

Nrun = 1;
st = 1;

for M = st:Nrun
    close all;
    clearvars -except Nrun M
    
    %% -------------------- Initialization --------------------
    ncell = 300;
    foldername = ['Viscous_' num2str(ncell)];
    if ~isfolder(foldername)
        mkdir(foldername)
    end

    % Physical conversions
    x_umpergrid = 1.25;   % micron per grid
   
    p.area = round(pi*5^2/x_umpergrid^2);  % cell area
    MeshN = 256;

    % Cell cluster initialization (stacked layers)
    ringin = 1;
    layer = 15;
    testcellcount = 1;
    x00 = MeshN/2; y00 = MeshN/2;
    nsig = zeros(MeshN, MeshN, 'int16');  % cell label
    ntyp = zeros(MeshN, MeshN, 'int16');  % cell type
    [x, y] = meshgrid(1:MeshN, 1:MeshN);

    for i = 1:layer
        radius = ringin + (i-1)*2*ceil(sqrt(p.area/pi));
        dtheta = 2*sqrt(p.area/pi)/radius;
        for j = 1:floor(2*pi/dtheta)+1
            if testcellcount > ncell
                break
            end
            theta = (j-1)*dtheta;
            x0 = x00 + radius*cos(theta);
            y0 = y00 + radius*sin(theta);

            cellmap = (sqrt((x-x0).^2 + (y-y0).^2) <= sqrt(p.area/pi)) & (nsig==0);
            nsig(cellmap) = testcellcount;
            ntyp(cellmap) = 1;
            testcellcount = testcellcount + 1;
        end
    end
    ncell = testcellcount - 1;
    disp([M, ncell])
    p.ncell = ncell;

    %% -------------------- Simulation parameters --------------------
    rng('shuffle');
    savemovie = 1;          % 1 to save AVI
    everyN = 200;             % save every N frame
    framerate = 10;         
    subplotN = 3;           
    maketiff = 0;           
    recordstep = 200;       
    tN = 1e5;               
    dx = x_umpergrid;       
    dt = 0.01;

    % Activator-inhibitor parameters
    factor = 1.25;
    q.J = ncell/20;
    % q.J = 5;
    q.thres = 0.7;
    q.rho = 1; 
    q.D = 8000;     
    DcAMP = 200/(x_umpergrid^2); 
    q.A = factor*2*5*0.0178*0.5*6/2/20; 
    q.Kg = factor*2*0.85*17.68*0.5/20; 
    q.Kr = factor*2*13*0.5/20; 
    q.a = 0.1863; 
    finaltau = 100;
    q.Kv = 2/40; 
    q.Kp = 4/40; 

    % Cell mechanical properties
    p.areamin = p.area;
    p.xlam = 50;
    p.vis = 0;
    p.mmt = -50;
    p.ecm = 5000;
    p.temp = 50;
    p.rseed = rand();
    p.dx = 2;               
    p.PottsIterations = 1;  
    p.stng = 0.04;          
    p.rThresh = 0.15;

    jbg = 0;
    p.jma = zeros(3,3);
    p.jma(1,2) = 700 + jbg;
    p.jma(2,1) = 700 + jbg;
    p.jma(2,2) = 50 + jbg;

    % Initialize arrays for cell metrics
    [nar, cmx, cmy, nl, Px, Py] = deal(zeros(1, p.ncell, 'single'));
    ECM = zeros(MeshN, MeshN, 'single');

    % Neighbor offsets
    nbxupd = [1,1,0,-1,-1,-1,0,1];
    nbyupd = [0,1,1,1,0,-1,-1,-1];

    % Center of mass of whole system
    CM = zeros(tN/recordstep,2);
    CM(1,:) = [sum(sum(x.*(nsig>0),1),2)/sum(nar), sum(sum(y.*(nsig>0),1),2)/sum(nar)];

    % FFT mesh
    kx = [0:MeshN/2 -MeshN/2+1:-1]*2*pi/MeshN;
    ky = kx;
    [KX, KY] = meshgrid(kx, ky);
    K2 = KX.^2 + KY.^2;

    % Movie setup
    if savemovie
        v = VideoWriter([foldername '/mov' num2str(M) '_' num2str(ncell/20) '.avi']);
        v.Quality = 99;
        v.FrameRate = framerate;
        open(v);
    end

    % Initialize concentrations 
    [g, r, s] = deal(zeros(1, p.ncell, 'single'));
    S = zeros(MeshN, MeshN);
    Sall = zeros(MeshN, MeshN, tN/recordstep);
    CMx_all = zeros(tN/recordstep, p.ncell);
    CMy_all = zeros(tN/recordstep, p.ncell);

    % Initial CM velocity noise
    Amp = 0;
    cmxOld = cmx - Amp*(rand(1,p.ncell)-0.5);
    cmyOld = cmy - Amp*(rand(1,p.ncell)-0.5);

    %% -------------------- Initial cell metrics --------------------
    for i = 1:p.ncell
        nar(i) = sum(nsig==i,'all');
        [rows, cols] = find(nsig==i);
        for j = 1:length(rows)
            tmp = sum([nsig(rows(j)+1,cols(j)), nsig(rows(j)-1,cols(j)), nsig(rows(j),cols(j)+1), nsig(rows(j),cols(j)-1)] ~= nsig(rows(j),cols(j)));
            nl(i) = nl(i) + tmp;
        end
        cmx(i) = sum(sum(x .* (nsig == i)))/nar(i);
        cmy(i) = sum(sum(y .* (nsig == i)))/nar(i);
    end

    %% -------------------- Main simulation loop --------------------
    FIG = figure('Position',[10 10 subplotN*512 512]);
    visible = true;
    hil = load('./hilresult.mat').hil;
    [cmphase, cmr] = cart2pol(cmx-mean(cmx,'all'), cmy-mean(cmy,'all'));
    intpg = interp1(hil.phase,hil.g,cmphase);
    intpr = interp1(hil.phase,hil.r,cmphase);
    
    g(:) = intpg;
    r(:) = intpr;
    S(:,:) = 0;

    

    for i = 2:tN

        % Concentration gradients

        Sx=circshift(S,[0,-1])-S;
        Sy=circshift(S,[-1,0])-S;
        
        for j = 1:p.ncell
            s(j) = sum(S(nsig == j));       % this line is slow
        end
        
        I = q.A*log10(1 + s./nar);
        
        % Activator-inhibitor dynamics
       
        nsig = nsig + 1;
        g = [0 g];
        dS = (g(nsig) > q.thres)*q.rho*q.D;
        dSfft = fft2(dS);
        Sk = fft2(S);
        Sknew = (Sk+dSfft*dt)./(1+DcAMP*K2*dt+q.J*dt);
        Snew = real(ifft2(Sknew));
       
        g(1) = [];
        nsig = nsig -1;
    
        dg = (-q.Kg*g.*(g-1).*(g-q.a) - q.Kr*r + I);
        dr = ((0.5*g-r)./finaltau);
    
        gnew = g + dg*dt;
        rnew = r + dr*dt;
    
        [thetaP, rhoP] = cart2pol((cmx - MeshN/2),(cmy - MeshN/2));
        CPM_step = 45;
        
    
        % This first 5000 iteration make cell form counter-clockwise rotation 
        if i < 5000
        
            P = 0.1;
            Px = -P.*sin(thetaP);
            Py = P.*cos(thetaP);
        end
    
        % This first 200 iteration make initial lattice like cell form single cluster 
    
        if i < 200
            
            velx=(cmx-cmxOld)/dt; vely=(cmy-cmyOld)/dt; vel=(velx.^2 + vely.^2).^(1/2);
            cmxOld=cmx; cmyOld=cmy;
            [nar, nsig, ntyp, nl, cmx, cmy] = CPM_evolve(cmx, cmy, velx, vely, p, dt, MeshN, ECM, S, Sx, Sy, r, nbxupd, nbyupd, Px, Py, nsig, nar, ntyp, nl);
        
        elseif mod(i, CPM_step)==0 

            velx=(cmx-cmxOld)/dt; vely=(cmy-cmyOld)/dt; vel=(velx.^2 + vely.^2).^(1/2);
            cmxOld=cmx; cmyOld=cmy;
            [nar, nsig, ntyp, nl, cmx, cmy] = CPM_evolve(cmx, cmy, velx, vely, p, dt, MeshN, ECM, S, Sx, Sy, r, nbxupd, nbyupd, Px, Py, nsig, nar, ntyp, nl);
            
            P = (Px.^2 + Py.^2).^(1/2);
            dPx = -q.Kp*Px + q.Kv*velx./(vel+1e-10).*(1-P);
            dPy = -q.Kp*Py + q.Kv*vely./(vel+1e-10).*(1-P);
    
            Px = Px + dPx*(CPM_step*dt);
            Py = Py + dPy*(CPM_step*dt);
    
        end

        S = Snew; g = gnew; r = rnew;
        S(S<0) = 0;
        
        
        if mod(i,200) == 0
            
            if subplotN > 1
                subplot(1,3,1)
                imagesc(Snew, [0, q.D/q.J])
    
                hold on
                set(gca,'YDir','normal')
    
                firstcell = nsig ==14;
                contour(firstcell, 'r')
                hold off
                title(num2str(i))
            end
            
            if subplotN > 1
                subplot(1,3,2)
                imagesc(nsig*(2^15*0.8/p.ncell)+int16((nsig>0)*2^15*0.2)+int16(ECM*2^15*0.05),[0 2^15])
                hold on
                set(gca,'YDir','normal')
                hold on 
                quiver([1 0 cmx],[0 0 cmy],[1 0 Px],[0 0 Py],0.5,'color','green','LineWidth',1)
                subplot(1,3,3)
                nsig = nsig + 1;
                r = [0 r];
                imagesc( r(nsig))
                r(1) = [];
                nsig = nsig -1;
                set(gca,'YDir','normal')
                hold off
            else
                imagesc(nsig*(2^15*0.8/p.ncell)+int16((nsig>0)*2^15*0.2)+int16(ECM*2^15*0.05),[0 2^15])
                hold on 
                quiver([1 0 cmx],[0 0 cmy],[1 0 velx],[0 0 vely],0.5,'color','white','LineWidth',1)
                hold off
            end
            title(num2str(i))
            
            if savemovie && (mod(i,everyN) == 0)
                drawnow; %Force the figure to render
                frame = getframe(FIG); %Convert the figure to a movie frame of figure
                writeVideo(v, frame); %Write the frame to the movie file
            end
        end

        % Record data
    
        if mod(i, recordstep) == 0
            Sall(:,:,i/recordstep) = S;
            CMx_all(i/recordstep,:) = cmx;
            CMy_all(i/recordstep,:) = cmy;
            CM(i/recordstep,1) = mean(cmx);
            CM(i/recordstep,2) = mean(cmy);
        end
        
    
    end
    
    if savemovie
        close(v);
    end

    %% -------------------- Post-processing --------------------
    Reff = sqrt(p.ncell*p.area/0.91/pi);
    entire.S = S;
    
    entire.cmx = CMx_all;
    entire.cmy = CMy_all;
    entire.recordstep = recordstep;
    
    %% -------------------- Calculating spiral arm number --------------------
    dr = 3;
    rmin = 0;
    rmax = sqrt(p.area*p.ncell/pi)*1.5;
    N = ceil((rmax - rmin)/dr);
    rclus = sqrt(p.area*p.ncell/pi);
    ravg = zeros(1,N);
    
    for i = 1:N
       rlb = rmin + dr*(i-1);
       rub = rmin + dr*(i);
       ravg(i) =(rlb+rub)/2;
    end

    entire.reff = ravg*x_umpergrid/(rclus*x_umpergrid);

    MeshN = size(S, 1);
    [x, y] = meshgrid(1:MeshN, 1:MeshN);
    [theta, r] = cart2pol(x-MeshN/2, y-MeshN/2);
    rrange = (r > 0.75*Reff) & ( r < 0.85*Reff);
    Srange = S(rrange);
    thetarange = theta(rrange);
    thetasplit = 20;
    thetacoarse = -pi:2*pi/thetasplit:pi;
    Scoarse = zeros(1, thetasplit);
    for k = 1:thetasplit
        thetaid = (thetarange>thetacoarse(k)) & (thetarange<thetacoarse(k+1));
        Scoarse(1, k) = mean(Srange(thetaid), 'all');
    end
    
    fftScoarse =fft(Scoarse-mean(Scoarse));
    fftScoarsehalf = fftScoarse(1:size(fftScoarse, 2)/2);
    
    peakid = find(fftScoarsehalf == max(fftScoarsehalf,[], 'all'));
    spiralnum = peakid-1;

    %% -------------------- Calculating period --------------------
    hexa_packing = 0.91;
    Aeff = p.ncell*p.area/hexa_packing;
    Reff = sqrt(Aeff/pi);
    [theta, R] = cart2pol(x-MeshN/2, y-MeshN/2);
    corrperiod_diffr = zeros(1, size(ravg, 2));
    for ri = 1:size(ravg, 2)
        if ri == 1
            rin = 0;
        else
            rin = (ravg(ri)+ravg(ri-1))/2;
        end
    
        if ri == size(ravg,2)
            rout = Reff;
        else 
            rout = (ravg(ri)+ravg(ri+1))/2;
        end
    
        idr2 = (R > rin) & (R < rout); 
        corrperiod = [];
        for i = 1:MeshN
            for j =1:MeshN
                if (nsig(i, j)~=0) && (idr2(i,j) == 1)
                    [R2, lags] = xcorr(squeeze(Sall(i, j, 1:end)) - mean(squeeze(Sall(i, j , 1:end))), 'coeff');
                    
                    R2pos = R2(lags > 0);
                    [~, locs] = findpeaks(R2pos);
                    for k = 1:size(locs, 1)
                        if R2pos(locs(k)) > 0
                            corrperiod(end+1) = locs(k) *dt*recordstep;
                            break
                        end
                    end
                end
            end
        end
        if size(corrperiod, 2) ~= 0 
            corrperiod_diffr(1, ri) = mean(corrperiod, 'omitnan');
        end
    end

    entire.corrperiod_diffr = corrperiod_diffr;
    disp([p.ncell, M, mean(corrperiod_diffr((entire.reff>0.6) & (entire.reff<0.8)))])
    entire.S = S;
    entire.spiralnum = spiralnum;

    save([foldername '/entire_' num2str(ncell/20) '.mat'], "entire")
end
