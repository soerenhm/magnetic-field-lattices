include("core.jl")
import PyPlot as plt; plt.pygui(true)


λ = 1480e-9
θ = deg2rad(10)
ψ = deg2rad(5)

wfun = PlaneWave(λ, [ψ, π/2], [1.42, -1im+1])
wsh1 = PlaneWave(λ/2, [θ/2, 0], [1, 1im])
wsh2 = PlaneWave(λ/2, [θ/2, π], [1, -1im])

x = -5λ:0.1λ:5λ
y = copy(x)
z = 0.0
Efun = sample(wfun, x, y, z)
Esh1 = sample(wsh1, x, y, z)
Esh2 = sample(wsh2, x, y, z)
J = currentdensity(Efun, Esh1 .+ Esh2)[:,:,1]
Jx, Jy, Jz = map(n -> getindex.(J, n), 1:3)
B = magnetic_field(real.(J), x, y)
Bz = getindex.(B, 3)

plt.pcolormesh(x,y,Bz')