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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModel3Ds::myDarcy::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& dZones = D_[zoneI];
        const tensorField& fZones = F_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            const label j = this->fieldIndex(i);
            const tensor Cd =
                mu[celli]*dZones[j] + (rho[celli]*mag(U[celli]))*fZones[j];


            //Info << "\nmu ="<< mu << endl;
            //Info << "\nrho =" << rho[celli]  << endl;
            //Info << "\nmagU="<< mag(U[celli]) << endl;
            //Info << "\ndZones  size ="<< dZones.size() << endl;
            //Info << "\ndZones ="<< dZones[j] << endl;
            //Info << "\nfZones ="<< fZones[j] << endl;


            //Info << "Cd size =" << Cd.size() << endl;

	    //Info << "Cd =" << Cd << endl;

            const scalar isoCd = tr(Cd); //mag(diag(Cd)); //

 	   // Info << "\n traccia Cd ="<< isoCd << endl;

            Udiag[celli] += V[celli]*isoCd;

            //Info << "\n Udiag ciclo =" << Udiag << endl;

            Usource[celli] -= V[celli]*((Cd - I*isoCd) & U[celli]);

            //Info << "\nU ciclo ="<< U << endl;

        }
    }
//Info << "Udiag ="<< Udiag << endl;
//Info << "volume cella ="<< V << endl;
//Info << "Usource ="<< Usource << endl;


Info << "Udiag[0] ="<< Udiag[0] << endl;
//Info << "Udiag[1] ="<< Udiag[1] << endl;
//Info << "Udiag[2] ="<< Udiag[2] << endl;
//Info << "Udiag[3] ="<< Udiag[3] << endl;
//Info << "Udiag[4] ="<< Udiag[4] << endl;
//Info << "Udiag[5] ="<< Udiag[5] << endl;
Info << "Usource[0] ="<< Usource[0] << endl;


}


template<class RhoFieldType>
void Foam::porosityModel3Ds::myDarcy::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& dZones = D_[zoneI];
        const tensorField& fZones = F_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            Info << "\nsono nel secondo apply" << endl;
            const label celli = cells[i];
            const label j = this->fieldIndex(i);
            const tensor D = dZones[j];
            const tensor F = fZones[j];

            AU[celli] += mu[celli]*D + (rho[celli]*mag(U[celli]))*F;
        }
    }
}


// ************************************************************************* //
