#include <iostream>
#include "Point.h"
#include <functional>
#include <vector>
#include <cmath>


//The enum class with the possibles strategies I can use
enum class Strategy{
    ExponentialDecay,
    InverseDecay,
    ApproximateLineSearch
};

//The enum class with the possibles ways to evaluate the derivative
enum class Derivative{
    FiniteDifferences,
    ExactGradient
};

//The struct with data
struct data{
    Point x0;
    double tolerance_r;
    double tolerance_s;
    std::function<double(Point)> fun;
    std::function<std::vector<double>(Point, Derivative)>  dfun;
    int max_iteration;
    double alpha_0;
};


//The function
//RMK: This function is the particular one requested for the challenge
std::function<double(Point)> fun=[](Point const & x)->double{
    return (x.get_coordinate(0)*x.get_coordinate(1)+4*pow(x.get_coordinate(0),4)+x.get_coordinate(1)*x.get_coordinate(1)+3*x.get_coordinate(0));
};

//Evualuate the gradient of the function(the one write before)
//You give as impute the choice of the user about the method to evluate the gradient
std::function<std::vector<double>(Point ,Derivative)> dfun=[](Point const & x, Derivative method)->std::vector<double>{
std::vector<double> grad(2);
double h=0.01;
Point x_plush;
Point x_minush;
 switch(method){
    case Derivative::FiniteDifferences:
        for(int i=0; i<2; ++i){
            x_plush.set_coordinates({x.get_coordinate(0),x.get_coordinate(1)});
            x_minush.set_coordinates({x.get_coordinate(0),x.get_coordinate(1)});
            x_plush.set_coordinate(i,x.get_coordinate(i)+h);
            x_minush.set_coordinate(i,x.get_coordinate(i)-h);
            grad[i]=(fun(x_plush)-fun(x_minush))/(2*h);
        }
        break;
   case Derivative::ExactGradient:
        grad[0]=x.get_coordinate(1)+16*std::pow(x.get_coordinate(0),3)+3;
        grad[1]=x.get_coordinate(0)+2*x.get_coordinate(1);
        break;
    }
return grad;
};

//Help function to evaluate the norm of the gradient
double compute_norm_grad(std::vector<double> const & grad){
    double n=0.0;
    for(double i:grad){
        n+=i*i;
    }
    return std::sqrt(n);
}

//Help function to evaluate the new alpha at the k-esim iteration
template <Strategy S>
 double compute_new_alpha(data d, int k,Point xnew, Derivative method){
    double new_alpha=d.alpha_0;
    
    if constexpr(S==Strategy::ExponentialDecay){
        new_alpha=d.alpha_0*std::exp(-0.2*k);
       }
    
    else if constexpr(S== Strategy::InverseDecay){
        new_alpha=d.alpha_0/(1+0.2*k);
        }
    

   else if constexpr(S==Strategy::ApproximateLineSearch){
        Point x_to_evaluate;
        for(std::size_t i=0; i<xnew.get_dimension(); ++i){
            x_to_evaluate.set_coordinate(i,xnew.get_coordinate(i)-d.alpha_0*dfun(xnew,method)[i]);
        }
        while((fun(xnew)-fun(x_to_evaluate))<0.3*new_alpha*pow(compute_norm_grad(dfun(xnew,method)),2)){
            new_alpha=new_alpha/2;
            for(std::size_t i=0; i<xnew.get_dimension(); ++i)
                x_to_evaluate.set_coordinate(i,xnew.get_coordinate(i)-new_alpha*dfun(xnew,method)[i]);
        }
    }
    else {
        std::cout<<"Error!"<<std::endl;
    }

    return new_alpha;
}



//Function to evaluate the min by the gradient method 
Point argmin_gradient( data d,Derivative derivative){
        Point xnew=d.x0;
        Point xold=d.x0;
        int k=1; 
        double alpha=d.alpha_0;
        std::size_t dim=xold.get_dimension();
        std::vector<double> grad=dfun(xold,derivative);
        //set x1
        for(std::size_t i=0; i<dim; ++i){
            xnew.set_coordinate(i,(xold.get_coordinate(i)-alpha*dfun(xold, derivative)[i]));
        }
        //Update x
        while(k<=d.max_iteration && xnew.distance(xold)>=d.tolerance_s && compute_norm_grad(dfun(xold,derivative))>=d.tolerance_r){
            xold=xnew;
            for(std::size_t i=0; i<dim; ++i){
                xnew.set_coordinate(i,(xold.get_coordinate(i)-alpha*dfun(xold,derivative)[i]));
            }
            ++k;
            //!! Here I have to set which method I want to use
            alpha= compute_new_alpha<Strategy::InverseDecay>(d,k,xnew, derivative);
        }
        return xnew;
}


//Function to evaluate the min by the gradient method 
template <Strategy S>
Point argmin_momentum( data d,Derivative derivative){
    double alpha=d.alpha_0;
    Point x_old_old=d.x0;
    Point x_old=d.x0;
    Point x_new=d.x0;
    int k=2;
    if constexpr(S==Strategy::ApproximateLineSearch){
        std::cout<<"You can't use the Armjo rule to evaluate the minimum by the method of the momentum"<<std::endl;
    }
    
    else if constexpr(S==Strategy::ExponentialDecay || S==Strategy::InverseDecay){
    //initialize x1
    x_new.set_coordinate(0,x_old.get_coordinate(0)-alpha*dfun(x_old,derivative)[0]);
    x_new.set_coordinate(1,x_old.get_coordinate(1)-alpha*dfun(x_old,derivative)[0]);
    while(k<=d.max_iteration && x_new.distance(x_old)>=d.tolerance_s && compute_norm_grad(dfun(x_old,derivative))>=d.tolerance_r){
        x_old_old.set_coordinates(x_old.get_coordinates());
        x_old.set_coordinates(x_new.get_coordinates());
        for(std::size_t i=0; i<d.x0.get_dimension(); ++i){
                x_new.set_coordinate(i,(x_old.get_coordinate(i)-alpha*dfun(x_old,derivative)[i])+0.9*(x_old.get_coordinate(i)-x_old_old.get_coordinate(i)));
        }
        ++k;
        alpha=compute_new_alpha<Strategy::InverseDecay>(d,k,x_new,derivative);
    }
    }

    else {
        std::cout<<"Error!"<<std::endl;
    }


return x_new;
}

template <Strategy S>
Point argmin_nesterov(data d,Derivative derivative){
    double alpha=d.alpha_0;
    Point x_old_old=d.x0;
    Point x_old=d.x0;
    Point x_new=d.x0;
    Point y;
    int k=2;
     if constexpr(S==Strategy::ApproximateLineSearch){
        std::cout<<"You can't use the Armjo rule to evaluate the minimum by the method of the momentum"<<std::endl;
        return d.x0;
    }
    else if constexpr(S==Strategy::ExponentialDecay || S==Strategy::InverseDecay){
    //initialize x1
    x_new.set_coordinate(0,x_old.get_coordinate(0)-alpha*dfun(x_old,derivative)[0]);
    x_new.set_coordinate(1,x_old.get_coordinate(1)-alpha*dfun(x_old,derivative)[0]);
    while(k<=d.max_iteration && x_new.distance(x_old)>=d.tolerance_s && compute_norm_grad(dfun(x_old,derivative))>=d.tolerance_r){
        x_old_old.set_coordinates(x_old.get_coordinates());
        x_old.set_coordinates(x_new.get_coordinates());
        y.set_coordinates({x_old.get_coordinate(0)+0.9*(x_old.get_coordinate(0)-x_old_old.get_coordinate(0)), x_old.get_coordinate(1)+0.9*(x_old.get_coordinate(1)-x_old_old.get_coordinate(1))});
        for(std::size_t i=0; i<d.x0.get_dimension(); ++i){
            x_new.set_coordinate(i,y.get_coordinate(i)-alpha*dfun(y,derivative)[i]);
        }
        ++k;
        //ATTENTION! 
        // In this method we cannot apply Armijo rule since the direction dk cannot be guaranteed to be a descent direction
        //So, if you impose as strategy ApproximateLineSearch the result is not attendible
        alpha=compute_new_alpha<Strategy::InverseDecay>(d,k,x_new,derivative);
    }
    }

    else {
        std::cout<<"Error!"<<std::endl;
    }

return x_new;
}






int main() {
    data parameters;

    // Set the initial point here:
    parameters.x0.set_coordinate(0,0);
    parameters.x0.set_coordinate(1,0);

    //Set the tolerance_r here:
    parameters.tolerance_r=1e-6;

    //Set the tolerance_s here:
    parameters.tolerance_s=1e-6;

    //Set the number of max iterations
    parameters.max_iteration=100;

    //Set alpha_0
    parameters.alpha_0=0.01;

    //Assigne the function to the member of the struct
    parameters.fun=fun;
    parameters.dfun=dfun;

    //Let's choose the way to evaluate the derivative
    Derivative derivate;
    int der_input;
    std::cout << "Enter the derivative method number (1: Finite Differences, 2: Exact Gradient): ";
    std::cin >> der_input;

    switch(der_input) {
        case 1:
           derivate =Derivative::FiniteDifferences;
           break;
       case 2:
        derivate=Derivative::ExactGradient;
        break;
    default:
        std::cout << "Invalid strategy number" << std::endl;
        return 1;
}
   
    // Let's find the minimun point
    //ATTENTION: You have to set the method you want to use to evluate alpha also in the definition of those funxtions:
    //see lines 127,156,195
    Point min_grad = argmin_gradient(parameters,derivate);
    Point min_mom=argmin_momentum<Strategy::InverseDecay>(parameters,derivate);
    Point min_nesterov=argmin_nesterov<Strategy::InverseDecay>(parameters,derivate);



    // Let's print the result:
    std::cout<<"The point where f is minimized by the gradient method is:";
    min_grad.print();
    std::cout<<std::endl;
    std::cout<<"The point where f is minimized by the momentum method is:";
    min_mom.print();
    std::cout<<std::endl;
    std::cout<<"The point where f is minimized by the Nesterov method is:";
    min_nesterov.print();
    std::cout<<std::endl;

    return 0;
}


