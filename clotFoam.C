/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    clotFoam

Description
    Simulate thrombosis as described by the model of Taylor et al. 

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "kinematicMomentumTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    //Declaring fields for use in loops
    
    volScalarField A_C
    (
        IOobject
        (
            "A_C",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("A_C", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
    
    volScalarField A_M
    (
        IOobject
        (
            "A_M",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("A_M", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );
    
    volScalarField phi_f
    (
        IOobject
        (
            "phi_f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("phi_f", dimensionSet(0,0,0,0,0,0,0), 0.0)
    );
    
    volScalarField tau
    (
        IOobject
        (
            "tau",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tau", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    );
    
    volScalarField P_TSP
    (
        IOobject
        (
            "P_TSP",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("P_TSP", dimensionSet(0,0,0,0,0,0,0), 0.0)
    );

volScalarField beta_eps
(
    IOobject
    (
        "beta_eps",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("beta_eps", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Calculating wall shear stress\n" << endl;

volScalarField WSS
(
    IOobject
    (
        "WSS",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("WSS", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

volScalarField check
(
    IOobject
    (
        "check",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("check", dimensionSet(0,0,0,0,0,0,0), 0.0)
);
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //    
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            #include "UEqn.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        }

        laminarTransport.correct();
        turbulence->correct();

	//Species
	#include "speciesEqn.H"
	
	//Aggregation intensity
	#include "epsEqn.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
