
        set term post eps color enhanced 12
        set out "likelihood_projection_for_mass25.eps"
        set xlabel "{/=28 mass:25}"
        set title "{/=28 likelihood: projection for mass:25}" 
        a=125.0
        b=-2.0
        c=17.356749492
        f(x)=1/(2*b**2)*(x-a)**2+c
        fit f(x) 'likelihood_projection_for_mass25.dat' u 1:2:3 via a,b,c
        plot 'likelihood_projection_for_mass25.dat' u 1:2:3 with errorbar title "" , 1/(2*b**2)*(x-a)**2+c title "fit" 
