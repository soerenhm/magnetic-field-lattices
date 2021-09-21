include("core.jl")
import PyPlot as plt; plt.pygui(true)


λ = 1800e-9
θ = deg2rad(10)
ψ = deg2rad(5)

wfun = PlaneWave(λ, [ψ, π/2], [1.42, -1im+1])
wsh1 = PlaneWave(λ/2, [θ/2, 0], [1, 1im])
wsh2 = PlaneWave(λ/2, [-θ/2, 0], [1, -1im])

x = -10λ:0.1λ:10λ
y = copy(x)
z = 0.0
Efun = sample(wfun, x, y, z)
Esh1 = sample(wsh1, x, y, z)
Esh2 = sample(wsh2, x, y, z)
J = currentdensity(Efun, Esh1 .+ Esh2)[:,:,1]
Jx, Jy, Jz = map(n -> getindex.(J, n), 1:3)