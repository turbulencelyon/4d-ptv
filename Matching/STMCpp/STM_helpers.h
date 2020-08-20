/*
 * Miscellaneous tools used in STM
 *
 */

#ifndef STM_helpers
#define STM_helpers

#include <cmath>
#include <string>

template <typename number>
number specialdivision(number a,number b){
    if(b==0.0){
        return(-1*INFINITY);
    }
    else{
        return(a/b);
    }
}

template <typename number> // I don't check explicitly if it is a number, but it presumed to be a number (int, long, float, double...)
number squarenorm(number vector[3]) {
    number out = vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2];
    return(out);
}

template <typename number>
number norm(number vector[3]) {
    number out;
    out = squarenorm(vector);
    out = sqrt(out);
    return(out);
}

template <typename number>
bool atface(number bmin,number bmax,number b){
    return(bmin <=b && b <= bmax);
}

template <typename number>
signed int Sgn(number x)
{
    if (x<0)
        return(-1);
    else
        return(1);
}

std::string timeasstring();

#endif
