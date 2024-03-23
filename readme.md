
# The minimum of a function

In the following code, I'm trying to evaluate the minimum of a function by the gradient method, by the momentum method, and by the Nesterov method.

The main parameters (tolerances, initial_alpha, initial point, max iteration) are given and setted in the main. Also the function to use is given, and is setted by a function wrapper. All those are part of a struct called "data".

The alpha can be updated by different strategies: Exponential decay, Inverse decay and Approximate line search (with the Armijo rule).
To try to make the code more efficient, the strategy choice is a templete parameter, so the user has to set the strategy that he wants to use by setting it directly by the code; in particular he has to set it in the main function and at lines 127, 156, and 195 of the main.cpp file.
Note that with the momentum method and the Nesterov method you can't use the Armijo rule since the direction dk cannot be guaranteed to be a descent direction.

The code also gives the possibility of choosing whether to evaluate the gradient by the exact gradient or by the finite differences method. The user can make this choice run time, by giving it as an input (see line 237 of the main).

Because of the Makefile, for the user who wants to use my code it is sufficient to run the command make in the bash.
