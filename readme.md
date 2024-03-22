
# The minimum of a function

In the following code I'm trying to evaluate the minimum of a function by the gradient method, by the momentum method, and by the Nesterov method.

The main parameters (tollerances,initial_alpha, initial point,max iteration) are given and setted in the main. Also the function to use is given, and is setted by a function. All those are part of a struct called "data".

The alpha can be updated by different strategies: Exponential decay, Inverse decay and Approximate line search(with the Armijo rule).
To try to make the code more efficient, the strategy choice is a templete parameter, so the user have to set the strategy that he want to use by setting it directly by the code; in particulare he have to set it in the main function and at the lines 127,156 and 195 of the main.
Note that with the momentum method and the Nesterov method you can't use the Armijo rule since the direction dk cannot be guaranteed to be a descent direction.

The code also give the possibility to choice if to evaluate the gradient by the exact gradient or by the finite differences method. The user can do this choice run time, by giving it as an imput(see line 237 of the main).

Because of the Makefile, for the user who want to use my code is sufficiently to run the command make in the bash.

