function mu = gravParams(body)

%{ 
    This is a helper function that returns the gravitational parameter in
    km^3 s^-2for each planet given the body's code:

    Planet/Sun Codes:
        Me - Mercury
        V  - Venus
        E  - Earth
        M  - Mars
        J  - Jupiter
        S  - Saturn
        U  - Uranus
        N  - Neptune
        Su - Sun

%} 

    keys   = {'Me','V','E','M','J','S','U','N','Su'};
    values = [22031.868551, 324858.592000, 398600.435507, 42828.375816, ...
              126712764.100000, 37940584.841800, 5794556.400000,...
              6836527.100580, 1.3271244004193938e11];
    d  = dictionary(keys,values);
    mu = d(body);

end