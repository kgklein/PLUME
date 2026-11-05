        zeta=(delta-n*gyro)/(kpar*thermal)
        response=spec(is)%alph_s-1. + &
             (spec(is)%alph_s*delta+(1.-spec(is)%alph_s)*n*gyro)*zet_in(kpar,zeta)/(kpar*thermal)
