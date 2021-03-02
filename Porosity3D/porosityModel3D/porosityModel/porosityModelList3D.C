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

#include "porosityModelList3D.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModelList3D::porosityModelList3D
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<porosityModel3D>(),
    mesh_(mesh)
{
    reset(dict);

    active(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModelList3D::~porosityModelList3D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::porosityModelList3D::active(const bool warn) const
{
    bool a = false;
    forAll(*this, i)
    {
        a = a || this->operator[](i).active();
    }

    if (warn && this->size() && !a)
    {
        Info<< "No porosity models active" << endl;
    }

    return a;
}


void Foam::porosityModelList3D::reset(const dictionary& dict)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& modelDict = iter().dict();

            this->set
            (
                i++,
                porosityModel3D::New(name, mesh_, modelDict)
            );
        }
    }
}


bool Foam::porosityModelList3D::read(const dictionary& dict)
{
    bool allOk = true;
    forAll(*this, i)
    {
        porosityModel3D& pm = this->operator[](i);
        bool ok = pm.read(dict.subDict(pm.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


bool Foam::porosityModelList3D::writeData(Ostream& os) const
{
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    return os.good();
}


void Foam::porosityModelList3D::addResistance
(
    fvVectorMatrix& UEqn
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn);
    }
}


void Foam::porosityModelList3D::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, rho, mu);
    }
}


void Foam::porosityModelList3D::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC         
)
{
    forAll(*this, i)
    {
        this->operator[](i).addResistance(UEqn, AU, correctAUprocBC);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const porosityModelList3D& models
)
{
    models.writeData(os);
    return os;
}


// ************************************************************************* //
