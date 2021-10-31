# Fiber-Laser-model
Python script to simulate performance of an Ytterbium doped fiber laser/amplifier.

The model assumes a two-level system with the system operating in the steady state (dN2/dt = 0). 
By assuming steady state operation, it is possible to get an exact value for the inversion at every step along the active fiber.
This value is then used to determine the power(s) of the signal and pump respectively, from the laser rate equations. These coupled 
wave equations are solved via the Runge-Kutta 4 algorithm. 

The script MainProg.py is the main program which calls the RK4 routine and other additional functions.

![abs spec2](https://user-images.githubusercontent.com/93448334/139584712-50debaf4-b488-48a3-836e-123efaa24e7d.png)

![amp model3](https://user-images.githubusercontent.com/93448334/139585502-42d80781-74d2-40fa-94a2-c49e8ff7d363.png)

