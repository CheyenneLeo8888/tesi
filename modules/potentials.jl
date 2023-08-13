module Potentials

export harmonic_oscillator, hydrogen_atom, particle_in_a_box, double_well, coulomb_potential, lennard_jones_potential,
morse_potential

function harmonic_oscillator(x)
    return x^2
end

function hydrogen_atom(x, k)
    r = abs(x)
    return -k / r
end

function particle_in_a_box(x)
    return 0  
end

function double_well(x, a, b, V0)
    if a <= x <= b
        return 0 
    else
        return V0  
    end
end

function coulomb_potential(x, q1, q2, k)
    r = abs(x)
    return -k * q1 * q2 / r
end

function lennard_jones_potential(x, epsilon, sigma)
    r = abs(x)
    return 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6)
end

function morse_potential(x, D, a, re)
    r = abs(x)
    return D * (1 - exp(-a * (r - re)))^2
end

end  # module