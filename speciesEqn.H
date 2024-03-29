/******************************************************************************
The activation of platelets due to the activity of ADP is quantified. 
A threshold ADP concentration must be achieved before the activation begins.
******************************************************************************/
Info<< "\nCalculating chemical activation\n" << endl;

forAll(A_C, i)
{
    if(ADP[i] < ADP_t.value())
        A_C[i] = 0.0;
    else
        A_C[i] = ADP[i]/ADP_t.value()/t_ADP.value();
}

Info<< "Max/min A_C: " << max(A_C).value() << "/" << min(A_C).value() << nl << endl;

/******************************************************************************
Calculate shear stress by multiplying the dynamic viscocity with two times the 
symmetric part of the local velocity gradient. The magniture of the tensor is 
then calculated and this scalar is used in a power-law approximation to 
estimate the contribution on shear to platelet activation.
******************************************************************************/

Info<< "\nCalculating mechanical activation\n" << endl;

volTensorField tauTensor = nu*rho*fvc::grad(U);
tauTensor += tauTensor.T();

forAll(tau, cellI)
{
    tau[cellI] = tauTensor[cellI].xx()*tauTensor[cellI].xx() + tauTensor[cellI].yy()*tauTensor[cellI].yy() + tauTensor[cellI].zz()*tauTensor[cellI].zz() - tauTensor[cellI].xx()*tauTensor[cellI].yy() - tauTensor[cellI].xx()*tauTensor[cellI].zz() - tauTensor[cellI].yy()*tauTensor[cellI].zz() + 3.0*(tauTensor[cellI].xy()*tauTensor[cellI].xy() + tauTensor[cellI].xz()*tauTensor[cellI].xz() + tauTensor[cellI].yz()*tauTensor[cellI].yz());
    
    tau[cellI] = Foam::sqrt(tau[cellI]/3.0);
}  

Info<< "Max/min shear stress: " << max(tau).value() << "/" << min(tau).value() << nl << endl;

phi_f = phi_a/(phi_a + phi_n);    //Fraction of activated platelets

forAll(A_M, cellI)
{
    A_M[cellI] = Foam::pow(C.value(), 1.0/beta.value())*beta.value()*Foam::pow(phi_f[cellI], (beta.value() - 1.0)/beta.value())*Foam::pow(tau[cellI], alpha.value()/beta.value());
}

Info<< "Max/min A_M: " << max(A_M).value() << "/" << min(A_M).value() << nl << endl;

/******************************************************************************
Solve the three species equations.
******************************************************************************/

//Inactivated platelets
Info<< "\nSolving non-activated platelet equation" << endl;
solve(fvm::ddt(phi_n) + fvm::div(phi, phi_n) == fvm::laplacian(D_n, phi_n) - (A_C + A_M)*phi_n);
Info<< "Max/min non-activated platelets: " << max(phi_n).value() << "/" << min(phi_n).value() << nl << endl;

//Activated platelets
Info<< "\nSolving activated platelet equation" << endl;
solve(fvm::ddt(phi_a) + fvm::div(phi, phi_a) == fvm::laplacian(D_a, phi_a) + (A_C + A_M)*phi_n);
Info<< "Max/min activated platelets: " << max(phi_a).value() << "/" << min(phi_a).value() << nl << endl;

//ADP
Info<< "\nSolving ADP equation" << endl;
solve(fvm::ddt(ADP) + fvm::div(phi, ADP) == fvm::laplacian(D_ADP, ADP) + R_ADP*(A_C + A_M)*phi_n);
Info<< "Max/min ADP: " << max(ADP).value() << "/" << min(ADP).value() << nl << endl;

