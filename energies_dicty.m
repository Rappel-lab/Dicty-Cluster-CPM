function dE = energies_dicty(p, xran, yran, xnb, ynb, MeshN, nsig, ntyp, nar, nl, cmx, cmy, velx, vely, Px, Py, conx, cony, nbxupd, nbyupd, S, ECM)
    
	nold = ntyp(yran,xran);
	nnew = ntyp(ynb,xnb);
    
    % Cell-Cell Adhesion
    enold = 0;
    for i=1:8
        if nsig(yran,xran) ~= nsig(yran+nbyupd(i),xran+nbxupd(i))
           enold = enold + p.jma(nold+1,ntyp(yran+nbyupd(i),xran+nbxupd(i))+1);
        end
    end
    
    ennew = 0;
    for i=1:8
        if nsig(ynb,xnb) ~= nsig(yran+nbyupd(i),xran+nbxupd(i))
           ennew = ennew + p.jma(nnew+1,ntyp(yran+nbyupd(i),xran+nbxupd(i))+1);
        end
    end

	% Area constraint
	if (nold~=1)
		fac2 = nar(nsig(ynb,xnb)) - p.area;
		etold = fac2*fac2;
		fac2 = nar(nsig(ynb,xnb)) + 1 - p.area;
		etnew = fac2*fac2;
    else if (nnew~=1)
			fac1 = nar(nsig(yran,xran)) - p.area;
			etold = fac1*fac1;
			fac1 = nar(nsig(yran,xran)) - 1 - p.area;
			etnew = fac1*fac1;
        else
			fac1= nar(nsig(yran,xran)) - p.area;
			fac2 = nar(nsig(ynb,xnb)) - p.area;
			etold = fac1*fac1 + fac2*fac2;
			fac1 = nar(nsig(yran,xran)) - 1 - p.area;
			fac2 = nar(nsig(ynb,xnb)) + 1 - p.area;
			etnew = fac1*fac1 + fac2*fac2;
        end
    end
 
	% chemotaxis
    if (nsig(yran,xran) ~= 0)
        epold = conx(nsig(yran,xran))*(xran - cmx(nsig(yran,xran))) + cony(nsig(yran,xran))*(yran - cmy(nsig(yran,xran)));

        if (nsig(ynb,xnb) ~= 0)
            epnew = conx(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + cony(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
        else
            epnew = 0;
        end
    else
        epold = 0;
        epnew = conx(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + cony(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
    end
    
   

      
    % Polarization
    if (nsig(yran,xran) ~= 0)
        emnew1 = -Px(nsig(yran,xran))*(xran - cmx(nsig(yran,xran))) - Py(nsig(yran,xran))*(yran - cmy(nsig(yran,xran)));
        if (nsig(ynb,xnb) ~= 0)
            emnew2 = Px(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + Py(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
        else
            emnew2 = 0;
        end
    else
        emnew1 = 0;
        emnew2 = Px(nsig(ynb,xnb))*(xran - cmx(nsig(ynb,xnb))) + Py(nsig(ynb,xnb))*(yran - cmy(nsig(ynb,xnb)));
    end



	eoldtot = enold + p.xlam*etold + epold;
	enewtot = ennew + p.xlam*etnew + epnew + p.mmt*(emnew1 + emnew2);


	dE = enewtot - eoldtot;
end