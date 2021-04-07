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

#include "addToRunTimeSelectionTable.H"
#include "myDarcy.H"
#include "geometricOneField.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModel3Ds
    {
        defineTypeNameAndDebug(myDarcy, 0);
        addToRunTimeSelectionTable(porosityModel3D, myDarcy, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModel3Ds::myDarcy::myDarcy
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel3D(name, modelType, mesh, dict, cellZoneName),
    dXYZ_("d", dimless/sqr(dimLength), coeffs_),
    
    fXYZ_("f", dimless/dimLength, coeffs_),
    D_(cellZoneIDs_.size()),
    F_(cellZoneIDs_.size()),
    rhoName_(coeffs_.lookupOrDefault<word>("rho", "rho")),
    muName_(coeffs_.lookupOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffs_.lookupOrDefault<word>("nu", "nu"))
{
    //adjustNegativeResistance(dXYZ_);
    //adjustNegativeResistance(fXYZ_);
    //Info << "\nrho =" << rhoName_ << endl;

    calcTransformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModel3Ds::myDarcy::~myDarcy()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModel3Ds::myDarcy::calcTransformModelData()
{
    if (coordSys_.R().uniform())
    {
	Info << "\nsono in myDarcy\n" << endl;

        forAll(cellZoneIDs_, zoneI)
        {
            D_[zoneI].setSize(1);
            F_[zoneI].setSize(1);

            D_[zoneI][0] = Zero;
            D_[zoneI][0].xx() = dXYZ_.value().xx();     //!!!!!!!!!!!!!!  modifico assegnazioni
            D_[zoneI][0].xy() = dXYZ_.value().xy(); 
            D_[zoneI][0].xz() = dXYZ_.value().xz(); 
	    
            D_[zoneI][0].yx() = dXYZ_.value().yx();
            D_[zoneI][0].yy() = dXYZ_.value().yy(); 
            D_[zoneI][0].yz() = dXYZ_.value().yz(); 

            D_[zoneI][0].zx() = dXYZ_.value().zx();
            D_[zoneI][0].zy() = dXYZ_.value().zy(); 
            D_[zoneI][0].zz() = dXYZ_.value().zz(); 

            D_[zoneI][0] = coordSys_.R().transformTensor(D_[zoneI][0]);

            // leading 0.5 is from 1/2*rho
            
	    F_[zoneI][0] = Zero;
            
            F_[zoneI][0].xx() = 0.5*fXYZ_.value().xx();
	    F_[zoneI][0].xy() = 0.5*fXYZ_.value().xy();
	    F_[zoneI][0].xz() = 0.5*fXYZ_.value().xz();
	    
            F_[zoneI][0].yx() = 0.5*fXYZ_.value().yx();
	    F_[zoneI][0].yy() = 0.5*fXYZ_.value().yy();
	    F_[zoneI][0].yz() = 0.5*fXYZ_.value().yz();
            
            F_[zoneI][0].zx() = 0.5*fXYZ_.value().zx();
	    F_[zoneI][0].zy() = 0.5*fXYZ_.value().zy();
	    F_[zoneI][0].zz() = 0.5*fXYZ_.value().zz();

            F_[zoneI][0] = coordSys_.R().transformTensor(F_[zoneI][0]);
        
        }
	    Info << "D =" << D_[0] <<  endl;
	    Info << "F =" << F_[0] <<  endl;
    }
    else
    {
	Info << "non-uniform" << endl;

        forAll(cellZoneIDs_, zoneI)
        {
            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            D_[zoneI].setSize(cells.size());
            F_[zoneI].setSize(cells.size());

            forAll(cells, i)
            {
                D_[zoneI][i] = Zero;
	        
	        D_[zoneI][i].xx() = dXYZ_.value().xx();     //!!!!!!!!!!!!!!  modifico assegnazioni mettendo anche i al posto di 0
                D_[zoneI][i].xx() = dXYZ_.value().xy(); 
                D_[zoneI][i].xx() = dXYZ_.value().xz(); 
	    
                D_[zoneI][i].yy() = dXYZ_.value().yx();
                D_[zoneI][i].xx() = dXYZ_.value().yy(); 
                D_[zoneI][i].xx() = dXYZ_.value().yz(); 

                D_[zoneI][i].zz() = dXYZ_.value().zx();
                D_[zoneI][i].xx() = dXYZ_.value().zy(); 
                D_[zoneI][i].xx() = dXYZ_.value().zz(); 
                

                // leading 0.5 is from 1/2*rho
                F_[zoneI][i] = Zero;

                F_[zoneI][i].xx() = 0.5*fXYZ_.value().xx();
	        F_[zoneI][i].xy() = 0.5*fXYZ_.value().xy();
	        F_[zoneI][i].xz() = 0.5*fXYZ_.value().xz();
	    
                F_[zoneI][i].yx() = 0.5*fXYZ_.value().yx();
	        F_[zoneI][i].yy() = 0.5*fXYZ_.value().yy();
	        F_[zoneI][i].yz() = 0.5*fXYZ_.value().yz();
            
                F_[zoneI][i].zx() = 0.5*fXYZ_.value().zx();
	        F_[zoneI][i].zy() = 0.5*fXYZ_.value().zy();
	        F_[zoneI][i].zz() = 0.5*fXYZ_.value().zz();
                
            }

            const coordinateRotation& R = coordSys_.R(mesh_, cells);

            D_[zoneI] = R.transformTensor(D_[zoneI], cells);
            F_[zoneI] = R.transformTensor(F_[zoneI], cells);
		
        }
    }

    if (debug && mesh_.time().writeTime())
    {
	Info << "else-uniform" << endl;

        volTensorField Dout
        (
            IOobject
            (
                typeName + ":D",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("0", dXYZ_.dimensions(), Zero)
        );
        volTensorField Fout
        (
            IOobject
            (
                typeName + ":F",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor("0", fXYZ_.dimensions(), Zero)
        );

        UIndirectList<tensor>(Dout, mesh_.cellZones()[cellZoneIDs_[0]]) = D_[0];
        UIndirectList<tensor>(Fout, mesh_.cellZones()[cellZoneIDs_[0]]) = F_[0];

        Dout.write();
        Fout.write();
    }
}
      


void Foam::porosityModel3Ds::myDarcy::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
	Info << "1" << endl;
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();

    apply(Udiag, Usource, V, rho, mu, U);

    force = Udiag*U - Usource;

}


void Foam::porosityModel3Ds::myDarcy::correct
(
    fvVectorMatrix& UEqn
) const
{
    Info << "\nsono in myDarcy correct" << endl;

    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);
	Info << "2.1" << endl;
        if (mesh_.foundObject<volScalarField>(muName))
        {
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho, mu, U);
        }
        else
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho, rho*nu, U);
        }
    }
    else
    {
	Info << "\ncorrect 1 \n" << endl;

        if (mesh_.foundObject<volScalarField>(nuName))
        {  
            Info << "\ncorrect 2\n" << endl;
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);
 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
/*	int i=0;
        for(i=0;i<100;i++)
	{
		Info << "\nUdiag =" << Udiag[i] << endl;  
                Info << "\nUsource =" << Usource[i] << endl;  
        }   
*/

	 apply(Udiag, Usource, V, geometricOneField(), nu, U);
            
           
        }
        else
        {   Info << "\ncorrect 3\n" << endl;
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), mu/rho, U);
        }
    }
}


void Foam::porosityModel3Ds::myDarcy::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
	Info << "3" << endl;
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho, mu, U);
}


void Foam::porosityModel3Ds::myDarcy::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
	Info << "4" << endl;
    const volVectorField& U = UEqn.psi();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(rhoName);
        const volScalarField& mu = mesh_.lookupObject<volScalarField>(muName);

        apply(AU, rho, mu, U);
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const volScalarField& nu =
                mesh_.lookupObject<volScalarField>(nuName);

            apply(AU, geometricOneField(), nu, U);
        }
        else
        {
            const volScalarField& rho =
                mesh_.lookupObject<volScalarField>(rhoName);
            const volScalarField& mu =
                mesh_.lookupObject<volScalarField>(muName);

            apply(AU, geometricOneField(), mu/rho, U);
        }
    }
}


bool Foam::porosityModel3Ds::myDarcy::writeData(Ostream& os) const
{
	Info << "5" << endl;
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}





// ************************************************************************* //
