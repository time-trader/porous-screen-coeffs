/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "porosityModel3D.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porosityModel3D, 0);
    defineRunTimeSelectionTable(porosityModel3D, mesh);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

//tentativo di modifica
/********************************************************************************\
void Foam::porosityModel3D::adjustNegativeResistance(dimensionedTensor& resist)
{

    scalar maxCmpt = max(0, cmptMax(resist.value()));

    if (maxCmpt < 0)
    {
        FatalErrorInFunction
            << "Negative resistances are invalid, resistance = " << resist
            << exit(FatalError);
    }
    else
    {
        tensor& val = resist.value();
        for (label cmpt = 0; cmpt < tensor::nComponents; cmpt++)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}

\**************************************************************/

Foam::label Foam::porosityModel3D::fieldIndex(const label i) const
{
    label index = 0;
    if (!coordSys_.R().uniform())
    {
        index = i;
    }
    return index;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModel3D::porosityModel3D
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    active_(true),
    zoneName_(cellZoneName),
    cellZoneIDs_(),
    coordSys_(coordinateSystem::New(mesh, coeffs_))
{

    if (zoneName_ == word::null)
    {
        dict.readIfPresent("active", active_);
        dict_.lookup("cellZone") >> zoneName_;
    }

    cellZoneIDs_ = mesh_.cellZones().findIndices(zoneName_);
     Info << "\nsono in porosityModel3D.C\n"<< endl;
    Info<< "    creating porous zone: " << zoneName_ << endl;

    bool foundZone = !cellZoneIDs_.empty();
    reduce(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorInFunction
            << "cannot find porous cellZone " << zoneName_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModel3D::~porosityModel3D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModel3D::transformModelData()
{
    if (!mesh_.upToDatePoints(*this))
    {   Info << "\nsono in porosityModel3D.C transformModelData  \n" << endl;
        calcTransformModelData();

        // set model up-to-date wrt points
        mesh_.setUpToDatePoints(*this);
    }
}


Foam::tmp<Foam::vectorField> Foam::porosityModel3D::porosityModel3D::force
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    transformModelData();

    tmp<vectorField> tforce(new vectorField(U.size(), Zero));
    
    Info << "\nsono in force" << endl;
    if (!cellZoneIDs_.empty())
    {
        this->calcForce(U, rho, mu, tforce.ref());
    }

    return tforce;
}


void Foam::porosityModel3D::addResistance(fvVectorMatrix& UEqn)
{  
    Info << "\nsono in porosityModel3D.C AddResistance 1 \n" << endl;

    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    Info << "\nsono in porosityModel3D.C AddResistance 1 dopo trasformModeldata  \n" << endl;
    this->correct(UEqn);
    Info << "\nsono in porosityModel3D.C AddResistance 1 dopo correct \n" << endl;
}


void Foam::porosityModel3D::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{  Info << "\nsono in porosityModel3D.C AddResistance 2\n" << endl;
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, rho, mu);
}


void Foam::porosityModel3D::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
)
{   Info << "\nsono in porosityModel3D.C AddResistance 3 \n" << endl;
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, AU);

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


bool Foam::porosityModel3D::writeData(Ostream& os) const
{
    return true;
}


bool Foam::porosityModel3D::read(const dictionary& dict)
{
    dict.readIfPresent("active", active_);

    coeffs_ = dict.optionalSubDict(type() + "Coeffs");

    dict.lookup("cellZone") >> zoneName_;
    cellZoneIDs_ = mesh_.cellZones().findIndices(zoneName_);

    return true;
}


// ************************************************************************* //
