function flybyDateIdxs = findFlybyDates(departureIdx, arrivalIdx, depSwingby, swingbyArr, tol)

    %{
        Inputs:
            departureIdx - index that corresponds to the Earth departure
                date
            arrivalIdx   - index that corresponds to the Saturn arrival
                date
            depSwingby   - array of Earth departure - Jupiter swingby dates
                with the jupiter-centric arrival speed computed
            swingbyArr   - arr of the Jupiter departure - Saturn arrival
                with the jupiter- centric departure speed computed
            tol  - tolerance for difference between departure and arrival
            velocity

        Outputs:
            Indices of the swingby dates

    %}

        jupArrival   = depSwingby(departureIdx,:);
        jupDeparture = swingbyArr(:,arrivalIdx)';

        flybyDateIdxs = [];

        for i = 1:length(jupArrival)

            if abs(jupArrival(i) - jupDeparture(i)) <= tol

                flybyDateIdxs(end + 1) = i;
                
            end

        end

end