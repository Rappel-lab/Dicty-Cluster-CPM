    
function [nar, nsig, ntyp, nl, cmx, cmy] = CPM_evolve(cmxOld, cmyOld, velx, vely, p, dt, MeshN, ECM, S, Sx, Sy, r, nbxupd, nbyupd, Px, Py, oldnsig, oldnar, oldntyp, oldnl)
    
    
    [nar,cmx,cmy] = deal(zeros([1 p.ncell],'single'));
    nsig = deal(zeros([MeshN MeshN],'int16')); % cell label, cell type
    ntyp = deal(zeros([MeshN MeshN],'int16')); % cell label, cell type
    nl = deal(zeros([1 p.ncell],'single'));
    nar(:,:) = oldnar;
    nsig(:,:) = oldnsig;
    ntyp(:,:) = oldntyp;
    nl(:,:) = oldnl;
    cmx(:,:) = cmxOld;
    cmy(:,:) = cmyOld;
        
    conx = - p.stng*interp2(1:MeshN,1:MeshN,Sx,cmx,cmy) .* (r<p.rThresh);
    cony = - p.stng*interp2(1:MeshN,1:MeshN,Sy,cmx,cmy) .* (r<p.rThresh);
     
    RAN = ceil(rand([2 MeshN^2*p.PottsIterations])*(MeshN-4))+2; % wall
    RAN8 = ceil(rand([1 MeshN^2*p.PottsIterations])*8);
   
    for ii = 1:MeshN^2*p.PottsIterations
       % Randomly pick a site to flip       
       xran = RAN(1,ii);
       yran = RAN(2,ii);
       
       % Randomly pick a neighbour site
       nbidx = RAN8(ii);
       xnb = xran + nbxupd(nbidx);
       ynb = yran + nbyupd(nbidx);
       
       % Monte-Carlo update 
   
       
           if (nsig(yran,xran) ~= nsig(ynb,xnb)) && (ntyp(ynb,xnb)~=2) 
                
                dE = energies_dicty(p, xran, yran, xnb, ynb, MeshN, nsig, ntyp, nar, nl, cmx, cmy, velx, vely, Px, Py, conx, cony, nbxupd, nbyupd, S, ECM);
    	        prob = exp(-dE/p.temp);
		        ran = rand();
                if (ran <= prob)
                    
				    % update center of mass and area %
				    if (ntyp(yran,xran) ~= 0) 
					    term1 = nar(nsig(yran,xran))/(nar(nsig(yran,xran))-1)*cmx(nsig(yran,xran));
					    term2 = xran/(nar(nsig(yran,xran))-1);
					    cmx(nsig(yran,xran)) = term1 - term2;
					    term1 = nar(nsig(yran,xran))/(nar(nsig(yran,xran))-1)*cmy(nsig(yran,xran));
					    term2 = yran/(nar(nsig(yran,xran))-1);
					    cmy(nsig(yran,xran)) = term1 - term2;
                        
                        nar(nsig(yran,xran)) = nar(nsig(yran,xran)) - 1;
                        nl(nsig(yran,xran)) = nl(nsig(yran,xran)) + (sum([nsig(yran+1,xran),nsig(yran-1,xran), nsig(yran,xran+1), nsig(yran,xran-1)] ~= nsig(yran,xran))-2)*(-2);
                    end
                    
				    if (ntyp(ynb,xnb) ~= 0)
					    term1 = nar(nsig(ynb,xnb))/(nar(nsig(ynb,xnb))+1)*cmx(nsig(ynb,xnb));
					    term2 = xran/(nar(nsig(ynb,xnb))+1);
					    cmx(nsig(ynb,xnb)) = term1 + term2;
					    term1 = nar(nsig(ynb,xnb))/(nar(nsig(ynb,xnb))+1)*cmy(nsig(ynb,xnb));
					    term2 = yran/(nar(nsig(ynb,xnb))+1);
					    cmy(nsig(ynb,xnb)) = term1 + term2;
                        
                        nar(nsig(ynb,xnb)) = nar(nsig(ynb,xnb)) + 1;
                        nl(nsig(ynb,xnb)) = nl(nsig(ynb,xnb)) + (sum([nsig(yran+1,xran),nsig(yran-1,xran), nsig(yran,xran+1), nsig(yran,xran-1)] ~= nsig(ynb,xnb))-2)*2;
                    end
                    
				    % Flip it %
				    nsig(yran,xran)=nsig(ynb,xnb);
				    ntyp(yran,xran)=ntyp(ynb,xnb);
                end
        
		        
           end
    end
end