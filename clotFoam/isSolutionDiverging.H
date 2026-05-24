// Check if the pressure is negative
if ( min(p).value() < 0 && runTime.time().value() > 1.0 )
{
    
    Info<< "!!! Negative pressure, program terminating !!!" <<endl;
    Info<< "\n!!!!! min(p) = "<< min(p).value() <<" at t = " << 
    runTime.time().value() << " !!!!!"<< endl;
    FatalError<< nl << exit(FatalError);
}

// Check if Theta_T is larger than 1.01
if ( max(Theta_T).value() > 1.01 )
{
    
    Info<< "!!! Total Platelet Fraction Exceeds 1.01 !!!" <<endl;
    Info<< "\n!!!!! max(Theta_T) > 1.01 at t = " << 
    runTime.time().value() << " !!!!!"<< endl;
    FatalError<< nl << exit(FatalError);
}
