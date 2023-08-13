using LinearAlgebra
using Plots

function solve_schrodinger_1D(n, dx, xmax, potential_func)
    # Discretization
    N = Int(xmax / dx)
    x = range(-xmax, stop=xmax, length=N)

    # Potential energy function
    V = potential_func(x)

    # Hamiltonian matrix
    H = zeros(N, N)
    for i in 1:N
        H[i, i] = 2.0 / (dx^2) + V[i]
        if i > 1
            H[i, i-1] = -1.0 / (dx^2)
            H[i-1, i] = -1.0 / (dx^2)
        end
    end

    # Eigenvalue problem
    eigenvals, eigenvecs = eigen(H)
    energies = eigenvals[1:n]
    wavefunctions = eigenvecs[:, 1:n]

    # Normalization
    norm = dx * LinearAlgebra.norm(wavefunctions[:, 1])   # First wavefunction as reference
    wavefunctions ./= sqrt(norm)

    return x, energies, wavefunctions
end

# Potential function for 1D quantum harmonic oscillator
function harmonic_oscillator_potential(x)
    return 0.5 * x.^2
end

# Parameters
n = 5       # Number of energy levels to compute
dx = 0.1    # Spatial step size
xmax = 5.0  # Maximum position

# Solve Schrödinger equation for 1D quantum harmonic oscillator
x, energies, wavefunctions = solve_schrodinger_1D(n, dx, xmax, harmonic_oscillator_potential)

# Plot wavefunctions
p = plot(x, abs2.(wavefunctions[:, 1:n]), title="Quantum Harmonic Oscillator", label=["E$i = $(energies[i])" for i in 1:n])

display(p)
const output = IOBuffer()
using REPL
const out_terminal = REPL.Terminals.TerminalBuffer(output)
const basic_repl = REPL.BasicREPL(out_terminal)
const basic_display = REPL.REPLDisplay(basic_repl)
Base.pushdisplay(basic_display)
readline()
