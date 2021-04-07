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

Class
    Foam::IOporosityModelList3D

Description
    List of porosity models with IO functionality

SourceFiles
    IOporosityModelList3D.C

\*---------------------------------------------------------------------------*/

#ifndef IOporosityModelList_H
#define IOporosityModelList_H

#include "IOdictionary.H"
#include "porosityModelList3D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class IOporosityModelList3D Declaration
\*---------------------------------------------------------------------------*/

class IOporosityModelList3D
:
    public IOdictionary,
    public porosityModelList3D
{
private:

    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Disallow default bitwise copy construct
        IOporosityModelList3D(const IOporosityModelList3D&);

        //- Disallow default bitwise assignment
        void operator=(const IOporosityModelList3D&);


public:

    // Constructors

        //- Construct from mesh
        IOporosityModelList3D(const fvMesh& mesh);


        //- Destructor
        virtual ~IOporosityModelList3D()
        {}


    // Member Functions

        //- Read dictionary
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
