Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field phi_n\n" << endl;
volScalarField phi_n
(
    IOobject
    (
        "phi_n",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field phi_a\n" << endl;
volScalarField phi_a
(
    IOobject
    (
        "phi_a",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field ADP\n" << endl;
volScalarField ADP
(
    IOobject
    (
        "ADP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field epsilon\n" << endl;
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

dimensionedScalar D_a
(
    "D_a",
    transportProperties.lookup("D_a")
);

dimensionedScalar D_ADP
(
    "D_ADP",
    transportProperties.lookup("D_ADP")
);

dimensionedScalar D_n
(
    "D_n",
    transportProperties.lookup("D_n")
);

dimensionedScalar C
(
    "C",
    transportProperties.lookup("C")
);

dimensionedScalar alpha
(
    "alpha",
    transportProperties.lookup("alpha")
);

dimensionedScalar beta
(
    "beta",
    transportProperties.lookup("beta")
);

dimensionedScalar kCon
(
    "kCon",
    transportProperties.lookup("kCon")
);

dimensionedScalar A
(
    "A",
    transportProperties.lookup("A")
);

dimensionedScalar alpha_eps
(
    "alpha_eps",
    transportProperties.lookup("alpha_eps")
);

dimensionedScalar B
(
    "B",
    transportProperties.lookup("B")
);

dimensionedScalar ADP_t
(
    "ADP_t",
    transportProperties.lookup("ADP_t")
);

dimensionedScalar t_ADP
(
    "t_ADP",
    transportProperties.lookup("t_ADP")
);

dimensionedScalar tau_low_wall
(
    "tau_low_wall",
    transportProperties.lookup("tau_low_wall")
);

dimensionedScalar tau_high_wall
(
    "tau_high_wall",
    transportProperties.lookup("tau_high_wall")
);

dimensionedScalar tau_bd_wall
(
    "tau_bd_wall",
    transportProperties.lookup("tau_bd_wall")
);

dimensionedScalar tau_low_thr
(
    "tau_low_thr",
    transportProperties.lookup("tau_low_thr")
);

dimensionedScalar tau_high_thr
(
    "tau_high_thr",
    transportProperties.lookup("tau_high_thr")
);

dimensionedScalar tau_bd_thr
(
    "tau_bd_thr",
    transportProperties.lookup("tau_bd_thr")
);

dimensionedScalar rho
(
    "rho",
    transportProperties.lookup("rho")
);

dimensionedScalar eps_t
(
     "eps_t",
     transportProperties.lookup("eps_t")
);

dimensionedScalar phi_ai
(
     "phi_ai",
     transportProperties.lookup("phi_ai")
);

dimensionedScalar phi_ni
(
     "phi_ni",
     transportProperties.lookup("phi_ni")
);

dimensionedScalar R_ADP
(
     "R_ADP",
     transportProperties.lookup("R_ADP")
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
