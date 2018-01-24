#include "pfHome.h"

double pfHome::randNormal(){
    static int have_number = 0;
    static double nd2 = 0.0;

    double x1 = 0.0;
    double x2 = 0.0;
    double sqr = 0.0;

    if (! have_number){
        do{
            x1 = 2.0 * randUniform(0, 1) - 1.0;
            x2 = 2.0 * randUniform(0, 1) - 1.0;

            sqr = x1 * x1 + x2 * x2;
        }while(sqr >= 1.0 || sqr == 0.0);

        // Box Muller Transformation
        double cnst = sqrt(-2.0 * log(sqr) / sqr);

        nd2 = x2 * cnst;
        have_number = 1;

        return x1 * cnst;

    }else{
         have_number = 0;
         return nd2;
    }
}

double pfHome::randUniform(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    double rand = dis(gen);
    return rand;
}

double pfHome::randUniform(const double min, const double max){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    double rand = dis(gen);
    return rand;
}