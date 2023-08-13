using LinearAlgebra;
using Plots;
include("modules/potentials.jl");
using .Potentials ;


# Get user input for potential selection and constants
println("Select a potential:")
println("1. Harmonic Oscillator")
println("2. Hydrogen Atom")
println("3. Particle in a Box")
println("4. Double-Well")
println("5. Coulomb Potential")
println("6. Lennard-Jones Potential")
println("7. Morse Potential")
potential_choice = parse(Int, readline())

# Create the potential function based on the user's choice
if potential_choice == 1
    V_function = Potentials.harmonic_oscillator
elseif potential_choice == 2
    println("Enter the constant 'k' for the Hydrogen Atom potential: ")
    k = parse(Float64, readline())
    V_function = x -> Potentials.hydrogen_atom(x, k)
elseif potential_choice == 3
    V_function = Potentials.particle_in_a_box
elseif potential_choice == 4
    println("Enter the left boundary 'a' of the wells: ")
    a = parse(Float64, readline())
    
    println("Enter the right boundary 'b' of the wells: ")
    b = parse(Float64, readline())
    
    println("Enter the potential energy 'V0' outside the wells: ")
    V0 = parse(Float64, readline())
    
    V_function = x -> Potentials.double_well(x, a, b, V0)
elseif potential_choice == 5
    println("Enter the charge 'q1': ")
    q1 = parse(Float64, readline())
    
    println("Enter the charge 'q2': ")
    q2 = parse(Float64, readline())
    
    println("Enter the constant 'k' for the Coulomb Potential: ")
    k = parse(Float64, readline())
    
    V_function = x -> Potentials.coulomb_potential(x, q1, q2, k)
elseif potential_choice == 6
    println("Enter the constant 'epsilon' for the Lennard-Jones Potential: ")
    epsilon = parse(Float64, readline())
    
    println("Enter the constant 'sigma' for the Lennard-Jones Potential: ")
    sigma = parse(Float64, readline())
    
    V_function = x -> Potentials.lennard_jones_potential(x, epsilon, sigma)
elseif potential_choice == 7
    println("Enter the constant 'D' for the Morse Potential: ")
    D = parse(Float64, readline())
    
    println("Enter the constant 'a' for the Morse Potential: ")
    a = parse(Float64, readline())
    
    println("Enter the bond length 're' for the Morse Potential: ")
    re = parse(Float64, readline())
    
    V_function = x -> Potentials.morse_potential(x, D, a, re)
end
# Let's enter the boundaries of our system
a = -6;
b = 6;
N = 5000;
# Let's create our grid of points: h will be the step of the range, i.e. our shortest distance 
x = range(a, stop=b, length=N);
h = x[2]-x[1];
# We create the matrix for the Kinetic energy
T = zeros(N-2, N-2)

for i in 1:N-2
    for j in 1:N-2
        if i == j
            T[i, j] = -2
        elseif abs(i - j) == 1
            T[i, j] = 1
        end
    end
end
# We create the matrix for the potential
V = zeros(N-2, N-2)
for i in 1:N-2
    for j in 1:N-2
        if i == j
            V[i,j]= V_function(x[i+1])
        else
            V[i,j]=0
        end
    end
end
# Let's create the Hamiltonian matrix
H = -T/(2*h^2) + V
# Find eigenvalues(in ascending order) and eigenfunctions
val,vec=eigen(H)
z = sortperm(val)
z = z[1:4]
energies = val[z] ./ val[z][1]
println(energies)

# Set the plot size
p = plot(figsize=(12, 10))

# Loop through the indices in z
for i in 1:length(z)
    y = vcat(vec[:, z[i]], 0)
    y = vcat(0, y)
    x1 = range(0, stop=1, length=length(y))
    plot!(x1, y, lw=3, label="Eigenfunction $i")
end

# Set labels and title
xlabel!("x")
ylabel!("ψ(x)")
title!("Normalized Wavefunctions using FDM")

# Display the plot
display(p)
const output = IOBuffer()
using REPL
const out_terminal = REPL.Terminals.TerminalBuffer(output)
const basic_repl = REPL.BasicREPL(out_terminal)
const basic_display = REPL.REPLDisplay(basic_repl)
Base.pushdisplay(basic_display)
readline()
