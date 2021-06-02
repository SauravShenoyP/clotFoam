
WSS = nu*rho*mag(fvc::grad(U));

//Coefficient function calculation at walls

const fvPatchList& patches = mesh.boundary();

forAll(patches, patchI)
{
    const fvPatch& curPatch = patches[patchI];
    if(isType<wallFvPatch>(curPatch))
    {
        forAll(curPatch, faceI)
        {
            label faceCellI = curPatch.faceCells()[faceI];
            check[faceCellI] = 1.0;
            if (WSS[faceCellI] >= tau_high_wall.value())
            {
                P_TSP[faceCellI] = 0.0;
            }
            else if (WSS[faceCellI] <= tau_low_wall.value())
            {
                P_TSP[faceCellI] = 1.0;
            }
            else
            {
                P_TSP[faceCellI] = (WSS[faceCellI]/tau_high_wall.value())*(Foam::exp((WSS[faceCellI] - tau_low_wall.value())/(tau_high_wall.value() - tau_low_wall.value())) - 1.0);
                P_TSP[faceCellI] /= (Foam::constant::mathematical::e - 1.0);
                P_TSP[faceCellI] = 1.0 - P_TSP[faceCellI];
            }
               
            if (WSS[faceCellI] >= tau_bd_wall.value())
            {
                beta_eps[faceCellI] = B.value();
            }
            else
            {
                beta_eps[faceCellI] = 0.0;
            }
        }
    }
}

Info<< "Calculating surface shear stress\n" << endl;

forAll(epsilon, cellI)
{
    
    //TSP calculation for within/at the surface of the thrombus
    
    if (check[cellI] == 1.0)
        continue;
    
    if (epsilon[cellI] >= eps_t.value())
    {
        if (WSS[cellI] >= tau_high_thr.value())
        {
            P_TSP[cellI] = 0.0;
        }
        else if (WSS[cellI] <= tau_low_thr.value())
        {
            P_TSP[cellI] = 1.0;
        }
        else
        {
            P_TSP[cellI] = (WSS[cellI]/tau_high_thr.value())*(Foam::exp((WSS[cellI] - tau_low_thr.value())/(tau_high_thr.value() - tau_low_thr.value())) - 1.0);
            P_TSP[cellI] /= (Foam::constant::mathematical::e - 1.0);
            P_TSP[cellI] = 1.0 - P_TSP[cellI];
        }
    }
    else    
    {
        const labelList& adjCells = mesh.cellCells()[cellI];
        forAll(adjCells, i)
        {
            if (epsilon[adjCells[i]] >= eps_t.value())
            {
                if (WSS[cellI] >= tau_high_thr.value())
                {
                    P_TSP[cellI] = 0.0;
                }
                else if (WSS[cellI] <= tau_low_thr.value())
                {
                    P_TSP[cellI] = 1.0;
                }
                else
                {
                    P_TSP[cellI] = (WSS[cellI]/tau_high_thr.value())*(Foam::exp((WSS[cellI] - tau_low_thr.value())/(tau_high_thr.value() - tau_low_thr.value())) - 1.0);
                    P_TSP[cellI] /= (Foam::constant::mathematical::e - 1.0);
                    P_TSP[cellI] = 1.0 - P_TSP[cellI];
                }
                break;
            }
        }
    }
    
    //Breakdown
    
    if (WSS[cellI] >= tau_bd_thr.value())
    {
        beta_eps[cellI] = B.value();
    }
    else
    {
        beta_eps[cellI] = 0.0;
    }
}

//Solve

Info<< "Solving epsilon equation\n" << endl;

solve(fvm::ddt(epsilon) == alpha_eps*P_TSP*(phi_a*phi_a) - fvm::Sp(beta_eps, epsilon));

Info<< "Max/min epsilon: " << max(epsilon).value() << "/" << min(epsilon).value() << nl << endl;








