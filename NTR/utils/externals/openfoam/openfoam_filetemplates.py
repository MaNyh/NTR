
template_dict_cascadecase = {}

template_dict_cascadecase["U"] = """
/*--------------------------------*- C++ -*----------------------------------*
| =========                 |                                                 |
| \	    /  F ield         | OpenFOAM: The Open Source CFD Toolbox         |
|  \    /   O peration     | Version:  v1612+                                 |
|   \  /    A nd           | Web:	www.OpenFOAM.com                      |
|    \/     M anipulation  |                                                  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format	ascii;
    class	volVectorField;
    location    "0";
    object	U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


dimensions	[0 1 -1 0 0 0 0];

internalField   uniform (__uniform_initial_vectorfield__);


boundaryField
{
    INLET
    {
        __inletboundarycondition__
    }

    OUTLET
    {
     	type            waveTransmissive;
        gamma           1.4;
        value           uniform (__initial_outlet_value__); //startwert only

    }

    BLADE
    {
     	type            noSlip;
    }

    Z_lower
    {
        type            cyclic;
    }

    Z_upper
    {
        type            cyclic;
    }

    Y_lower
    {
        type            cyclic;
    }

    Y_upper
    {
        type            cyclic;
    }

}


// ************************************************************************* //

"""

template_dict_cascadecase["T"]="""
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1612+                                |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   uniform __uniform_initial_scalarfield__;

boundaryField
{
    INLET
    {
        type            waveTransmissive;
        gamma           1.4;
        fieldInf        __Tinlet__;
        lInf            0.01;
        value           uniform __Tinlet__;
    }
    OUTLET
    {
        type            waveTransmissive;
        gamma           1.4;
        value           uniform __Toutlet__;
    }
    BLADE
    {
        type            zeroGradient;
    }

    Z_lower
    {
        type            cyclic;
    }

    Z_upper
    {
        type            cyclic;
    }

    Y_lower
    {
        type            cyclic;
    }

    Y_upper
    {
        type            cyclic;
    }

}


// ************************************************************************* //
"""

template_dict_cascadecase["p"] = """
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1612+                                |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      p;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform __uniform_initial_scalarfield__;

boundaryField
{
    INLET
    {
        type            waveTransmissive;
        gamma           1.4;
        value           uniform __pinlet__;
        fieldInf        __pinlet__;
        lInf             __lInf__;             //a measure of how far away the far-field condition should be

    }
    OUTLET
    {
        type            waveTransmissive;
        gamma           1.4;
	    value		    uniform __initial_outlet_value__;
    }
    BLADE
    {
        type            zeroGradient;
    }

    Z_lower
    {
        type            cyclic;
    }

    Z_upper
    {
        type            cyclic;
    }

    Y_lower
    {
        type            cyclic;
    }

    Y_upper
    {
        type            cyclic;
    }


}
// ************************************************************************* //
"""

template_dict_cascadecase["alphat"] = """
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1612+                                |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alphat;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform __uniform_initial_scalarfield__;

boundaryField
{
    INLET
    {
        type            calculated;
        value           uniform __uniform_initial_scalarfield__;
    }
    OUTLET
    {
        type            calculated;
        value           uniform __uniform_initial_scalarfield__;
    }
    BLADE
    {
     	type            zeroGradient;
    }
    Z_lower
    {
        type            cyclic;
    }

    Z_upper
    {
        type            cyclic;
    }

    Y_lower
    {
        type            cyclic;
    }

    Y_upper
    {
        type            cyclic;
    }

}


// ************************************************************************* //
"""

template_dict_cascadecase["nut"] = """
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1612+                                |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      nut;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform __uniform_initial_scalarfield__;

boundaryField
{
    INLET
    {
        type            calculated;
        value           uniform __uniform_initial_scalarfield__;
    }
    OUTLET
    {
        type            calculated;
        value           uniform __uniform_initial_scalarfield__;
    }
    BLADE
    {
     	type            zeroGradient;
    }
    Z_lower
    {
        type            cyclic;
    }

    Z_upper
    {
        type            cyclic;
    }

    Y_lower
    {
        type            cyclic;
    }

    Y_upper
    {
        type            cyclic;
    }

}


// ************************************************************************* //
"""


template_dict_boundaryconditions  = {}

template_dict_boundaryconditions["turbulentDFSEMWaveTransmissiveInletMod"] = """        type            turbulentDFSEMWaveTransmissiveInletMod;
        boxScale        __boxScale__;
        d               __d__;
        coeffWakeU      (0 0 0);
        coeffWakeV      (0 0 0);
        nCellPerEddy    2;
        //readEddies      1;
        lastUPrime      uniform (1.6e-6 1e-06 1e-06);
        gamma           1.4;
        lInf            0.01;
        coeffU          3{0};
        coeffV          3{0};
        y_min           0;
        y_max           1;
        relax           0.0001;
        qMean           nonuniform List<vector> 8((1.0 1.0 1.0) (1.0 1.0 1.0) (1.0 1.0 1.0) (1.0 1.0 1.0) (1.0 1.0 1.0) (1.0 1.0 1.0) (1.0 1.0 1.0) (1 1 1));
        R               uniform (__Rinlet__);
        L               uniform __Linlet__;
        U               uniform (__Uinlet__);
        mapMethod       nearestCell;
	    value 		uniform (__Uinlet__);
	    fieldInf        uniform (__Uinlet__);"""
