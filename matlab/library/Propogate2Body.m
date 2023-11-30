function DState = Propogate2Body(State, mu)

% Parse state vector:
x  = State(1);
y  = State(2);
z  = State(3);
vx = State(4);
vy = State(5);
vz = State(6);

% Calculate derivatives
r     = [x;y;z];
r_mag = norm(r);
ax    = ((-mu/(r_mag)^3)*r(1)); 
ay    = ((-mu/(r_mag)^3)*r(2));
az    = ((-mu/(r_mag)^3)*r(3));


% Package derivatives
DState    = zeros(size(State));
DState(1) = vx;
DState(2) = vy;
DState(3) = vz;
DState(4) = ax;
DState(5) = ay;
DState(6) = az;

end