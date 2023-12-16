#include "additionalsolver.h"


macroParam AdditionalSolver::ExacRiemanSolver(macroParam left, macroParam right, double Gamma)
{
    double maxIteration = 40; // макс число итераций
    double TOL=1e-8;
    double lambda = 0; // линия на грани КО
    macroParam ret(left.mixture);;
    double left_soundspeed=sqrt ( Gamma*left.pressure/left.density );
    double right_soundspeed=sqrt( Gamma*right.pressure/right.density);

    double p_star= 0.5*(left.pressure+right.pressure) +
            0.125 * ( left.velocity_tau-right.velocity_tau ) *
            ( left.density+right.density ) *
            ( left_soundspeed+right_soundspeed );
    p_star=std::max(p_star,TOL);
    double pMin=std::min(left.pressure,right.pressure);
    double pMax=std::max(left.pressure,right.pressure);

    if ( p_star>pMax )
    {
        double temp1= sqrt ( ( 2.0/ ( Gamma+1.0 ) /left.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );
        double temp2= sqrt ( ( 2.0/ ( Gamma+1.0 ) /right.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        p_star= (temp1*left.pressure+temp2*right.pressure+ ( left.velocity_tau-right.velocity_tau ) ) / ( temp1+temp2 );
        p_star=std::max(p_star,TOL);
    }
    else if ( p_star<pMin )
    {
       double temp1= ( Gamma-1.0 ) / ( 2.0*Gamma );
       p_star= pow(( left_soundspeed+right_soundspeed+0.5*(Gamma-1.0 )*
                   ( left.velocity_tau-right.velocity_tau ) ) /
                   (left_soundspeed/pow(left.pressure,temp1) +
                   right_soundspeed/pow(right.pressure,temp1)), 1.0/temp1);
    }
    double f1 = 0, f2 = 0, f_d = 0 ;
    for(double iteration = 1;iteration < maxIteration; iteration++)
    {
        //LEFT
        double temp1 = sqrt ( ( 2.0/ ( Gamma+1.0 ) /left.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );

        if (p_star<=left.pressure)
            f1=2.0/ ( Gamma-1.0 ) *left_soundspeed*
                    (pow(p_star/left.pressure,(Gamma-1.0 )/(2.0*Gamma))- 1.0) ;
        else
            f1= ( p_star-left.pressure ) *temp1;
        if (p_star<=left.pressure)
            f_d= pow( p_star/left.pressure,-(Gamma+1.0 )/( 2.0*Gamma ))/
                    ( left.density*left_soundspeed );
        else
            f_d=temp1* ( 1.0-0.5* ( p_star-left.pressure ) /
                         ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );
        //RIGHT
        temp1 = sqrt ( ( 2.0/ ( Gamma+1.0 ) /right.density ) /
                       ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        if (p_star<=right.pressure)
            f2=2.0/ ( Gamma-1.0 ) *right_soundspeed*
                    (pow(p_star/right.pressure,(Gamma-1.0 )/(2.0*Gamma))- 1.0) ;
        else
            f2= ( p_star-right.pressure ) *temp1;
        if (p_star<=right.pressure)
            f_d= f_d + pow( p_star/right.pressure,-(Gamma+1.0 )/( 2.0*Gamma ))/
                    ( right.density*right_soundspeed );
        else
            f_d=f_d + temp1* ( 1.0-0.5* ( p_star-right.pressure ) /
                         ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        double p_new = p_star - (f1+f2 - (left.velocity_tau - right.velocity_tau))/f_d;
        if(abs(p_new - p_star)/(0.5*abs(p_new + p_star)) < TOL)
            break;
        p_star = p_new;
    }
    // calculate star speed */
    double star_speed=0.5* ( left.velocity_tau + right.velocity_tau ) +0.5* ( f2-f1 );
    double left_star_density, left_tail_speed, left_head_speed,
            right_star_density, right_tail_speed,right_head_speed;
    //LEFT
    if ( p_star>=left.pressure ) {
            // SHOCK
        left_star_density = left.density * ( p_star / left.pressure + ( Gamma-1.0 ) / ( Gamma+1.0 ) ) /
                ( ( Gamma-1.0 ) / ( Gamma+1.0 ) * p_star / left.pressure + 1.0 );
        left_tail_speed = left.velocity_tau -left_soundspeed * sqrt ( ( Gamma+1.0 ) / ( 2.0*Gamma ) * p_star/left.pressure +
                ( Gamma-1.0 ) / ( 2.0*Gamma ) );
        left_head_speed = left_tail_speed;
    }
    else // % left_wave_ == kRarefaction
    {
        left_star_density = left.density * pow(p_star/left.pressure,1.0/Gamma);
        left_head_speed = left.velocity_tau - left_soundspeed;
        left_tail_speed = star_speed - sqrt ( Gamma*p_star/left_star_density );
    }
    //RIGHT
    if ( p_star>=right.pressure )
    {
        right_star_density = right.density *
                            ( p_star / right.pressure + ( Gamma-1.0 ) / ( Gamma+1.0 ) ) /
                            ( ( Gamma-1.0 ) / ( Gamma+1.0 ) * p_star / right.pressure + 1.0 );
        right_tail_speed = right.velocity_tau +
           right_soundspeed * sqrt ( ( Gamma+1.0 ) / ( 2.0*Gamma ) * p_star/right.pressure +
           ( Gamma-1.0 ) / ( 2.0*Gamma ) );
        right_head_speed = right_tail_speed;
    }
    else // % right_wave_ == kRarefaction
    {
        right_star_density = right.density *  pow(p_star/right.pressure, 1.0/Gamma );
        right_head_speed = right.velocity_tau + right_soundspeed;
        right_tail_speed = star_speed + sqrt ( Gamma*p_star/right_star_density );
    }

    bool is_left_of_contact = lambda  < star_speed;

    if ( is_left_of_contact )
    {// % the u is left of contact discontinuity
        if ( p_star>=left.pressure )  //the left wave is a shock
        {
            if ( lambda < left_head_speed )
            { // the u is before the shock
                ret.density  = left.density;
                ret.velocity_tau = left.velocity_tau;
                ret.pressure = left.pressure;
            }
            else  //% the u is behind the shock
            {
                ret.density  = left_star_density;
                ret.velocity_tau = star_speed;
                ret.pressure = p_star;
            }
        }
        else // % the left wave is a rarefaction
        {
            if ( lambda < left_head_speed )//  % the u is before the rarefaction
            {
                ret.density  = left.density;
                ret.velocity_tau = left.velocity_tau;
                ret.pressure = left.pressure;
            }
            else
            {
                if ( lambda < left_tail_speed )//  % the u is inside the rarefaction
                {//% left_rarefaction (4.56)}
                    double temp1 = 2.0/ ( Gamma+1.0 ) + ( Gamma-1.0 ) / ( Gamma+1.0 )/left_soundspeed *(left.velocity_tau - lambda);
                    ret.density = left.density *  pow(temp1, 2.0/( Gamma-1.0 ));
                    ret.pressure = left.pressure * pow(temp1, 2.0*Gamma/ ( Gamma-1.0));
                    ret.velocity_tau = 2.0/ ( Gamma+1.0 ) * ( left_soundspeed + ( Gamma-1.0 ) /2.0*left.velocity_tau + lambda);
                }
                else//  % the u is after the rarefaction
                {
                    ret.density  = left_star_density;
                    ret.velocity_tau = star_speed;
                    ret.pressure = p_star;
                }
            }
        }
    }
    else// % the queried u is right of contact discontinuity
        //%------------------------------------------------------------------------
    {
        if ( p_star>=right.pressure )  //% the right wave is a shock
        {
            if ( lambda > right_head_speed )  //% the u is before the shock
            {
                ret.density  = right.density;
                ret.velocity_tau = right.velocity_tau;
                ret.pressure = right.pressure;
            }
            else  //% the u is behind the shock
            {
                ret.density  = right_star_density;
                ret.velocity_tau = star_speed;
                ret.pressure = p_star;
            }
        }
        else // % the right wave is a rarefaction
        {
            if ( lambda > right_head_speed ) // % the u is before the rarefaction
            {
                ret.density  = right.density;
                ret.velocity_tau = right.velocity_tau;
                ret.pressure = right.pressure;
            }
            else
            {
                if ( lambda > right_tail_speed ) // % the u is inside the rarefaction
                {
                    double temp1 =2.0/ ( Gamma+1.0 ) - ( Gamma-1.0 ) / ( Gamma+1.0 ) /right_soundspeed *(right.velocity_tau - lambda);
                    ret.density = right.density *  pow(temp1, 2.0/ ( Gamma-1.0 ) );
                    ret.pressure = right.pressure * pow(temp1, 2.0*Gamma/ ( Gamma-1.0 ) );
                    ret.velocity_tau = 2.0/ ( Gamma+1.0 ) * ( -right_soundspeed + ( Gamma-1.0 ) /2.0*right.velocity_tau + lambda);
                }
                else // % the u is after the rarefaction
                {
                    ret.density  = right_star_density;
                    ret.velocity_tau = star_speed;
                    ret.pressure = p_star;
                }
            }
        }
    }
    return ret;
}
conservative  AdditionalSolver::SolveEvolutionExplFirstOrder(vector<Matrix> F1, Matrix F2, Matrix F3,  vector<Matrix> U1old, Matrix U2old, Matrix U3old, double dt, double delta_h)
{
    int len = F2.size();
    Matrix rf; // вектор расстояний до ячеек
    for(double i = 0; i < len; i++)
        rf.push_back(i);
    rf = rf * delta_h;
    vector<Matrix> U1new(F1.size());
    for(size_t i = 0; i < len; i++)
    {
        U1new[i].resize(len + 1);
    }
    Matrix U2new (len + 1),U3new(len + 1);
    auto temp_rfL = rf; temp_rfL.removeLast();
    auto temp_rfR = rf; temp_rfR.removeFirst();
    auto temp = Matrix::REVERSE(temp_rfR - temp_rfL) * dt;
    for(int i = 0 ; i < temp.size(); i++)
    {
        for(size_t j = 0; j < len; j++)
        {
            U1new[j][i] = U1old[j][i+1] - temp[i]*(F1[j][i+1] - F1[j][i]);
        }
        U2new[i+1] = U2old[i+1] - temp[i]*(F2[i+1] - F2[i]);
        U3new[i+1] = U3old[i+1] - temp[i]*(F3[i+1] - F3[i]);
    }
    conservative res = std::make_pair(U1new,vector<Matrix> {U2new,U3new});
    return res;
}


