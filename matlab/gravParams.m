function mu = gravParams(body)

%{ 
    This is a helper function that returns the gravitational parameter in
    km^3 s^-2for each planet given the body's code:

    Planet/Sun Codes:
       1 Me - Mercury
       2 V  - Venus
       3 E  - Earth
       4 M  - Mars
       5 J  - Jupiter
       6 S  - Saturn
       7 U  - Uranus
       8 N  - Neptune
       0 Su - Sun

%} 

    % keys   = {'Me','V','E','M','J','S','U','N','Su'};
    keys   = [1, 2, 3, 4, 5, 6, 7, 8, 0];
    values = [22031.868551, 324858.592000, 398600.435507, 42828.375816, ...
              126712764.100000, 37940584.841800, 5794556.400000,...
              6836527.100580, 1.3271244004193938e11];
    d  = dictionary(keys,values);
    mu = d(body);

end